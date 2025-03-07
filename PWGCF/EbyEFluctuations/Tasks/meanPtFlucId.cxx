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
  Configurable<bool> cfgCorrection{"cfgCorrection", true, "Efficiency Correction"};
  Configurable<bool> cfgCorrectionPID{"cfgCorrectionPID", true, "ID particles Efficiency Correction"};
  Configurable<bool> cfgCorrectionPtRap{"cfgCorrectionPtRap", false, "Efficiency Correction for pT and eta"};
  Configurable<bool> cfgCorrectionPtRapPID{"cfgCorrectionPtRapPID", false, "ID particles Efficiency Correction for pT and eta"};
  Configurable<bool> cfgPidCut{"cfgPidCut", false, ""};
  Configurable<bool> cfgPDGCodeOnly{"cfgPDGCodeOnly", true, ""};
  Configurable<bool> cfgMCReco{"cfgMCReco", false, ""};
  Configurable<bool> cfgMCTruth{"cfgMCTruth", false, ""};
  Configurable<bool> cfgPosZ{"cfgPosZ", true, "Position Z"};
  Configurable<bool> cfgSel8{"cfgSel8", true, "Sel8 trigger"};
  Configurable<bool> cfgNoSameBunchPileup{"cfgNoSameBunchPileup", true, "kNoSameBunchPileup"};
  Configurable<bool> cfgIsVertexITSTPC{"cfgIsVertexITSTPC", true, "kIsVertexITSTPC"};
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
  Configurable<std::vector<float>> effPtRapValuesPi{"effPtRapValuesPi", {0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.5515, 0.5339, 0.5204, 0.5122, 0.5071, 0.5179, 0.4312, 0.2341, 0.2669, 0.2606, 0.2483, 0.2551, 0.2536, 0.2605, 0.2621, 0.2646, 0.2654, 0.2783, 0.2767, 0.2868, 0.2924, 0.2871, 0.3049, 0.3052, 0.3075, 0.3031, 0.3160, 0.2992, 0.3131, 0.3090, 0.3215, 0.3294, 0.3211, 0.3354, 0.3489, 0.3483, 0.0000, 0.5485, 0.5303, 0.5227, 0.5120, 0.5096, 0.5167, 0.5028, 0.2424, 0.2605, 0.2525, 0.2403, 0.2457, 0.2476, 0.2480, 0.2518, 0.2536, 0.2566, 0.2577, 0.2618, 0.2598, 0.2632, 0.2792, 0.2768, 0.2865, 0.2797, 0.2885, 0.2976, 0.3001, 0.2950, 0.2927, 0.2926, 0.3179, 0.3063, 0.3083, 0.3224, 0.3266, 0.0000, 0.5348, 0.5284, 0.5193, 0.5115, 0.5132, 0.5205, 0.5248, 0.2876, 0.2646, 0.2592, 0.2533, 0.2574, 0.2596, 0.2653, 0.2700, 0.2693, 0.2750, 0.2897, 0.2815, 0.2855, 0.2882, 0.2946, 0.3091, 0.3090, 0.3212, 0.3181, 0.3212, 0.3242, 0.3234, 0.3339, 0.3299, 0.3373, 0.3367, 0.3349, 0.3436, 0.3658, 0.0000, 0.4804, 0.5192, 0.5200, 0.5127, 0.5142, 0.5205, 0.5251, 0.3668, 0.2952, 0.2831, 0.2745, 0.2764, 0.2823, 0.2872, 0.2898, 0.2871, 0.2919, 0.2955, 0.2993, 0.3029, 0.3082, 0.3162, 0.3135, 0.3298, 0.3314, 0.3362, 0.3339, 0.3412, 0.3307, 0.3402, 0.3282, 0.3472, 0.3743, 0.3436, 0.3606, 0.3546, 0.0000, 0.4305, 0.5050, 0.5156, 0.5141, 0.5118, 0.5180, 0.5228, 0.4094, 0.2868, 0.2791, 0.2660, 0.2683, 0.2718, 0.2772, 0.2758, 0.2782, 0.2824, 0.2818, 0.2874, 0.2870, 0.2836, 0.2885, 0.3005, 0.3003, 0.3145, 0.3158, 0.3055, 0.3081, 0.3207, 0.3314, 0.3149, 0.3277, 0.3208, 0.3137, 0.3477, 0.3461, 0.0000, 0.4123, 0.4795, 0.5125, 0.5081, 0.5123, 0.5177, 0.5244, 0.4448, 0.2584, 0.2504, 0.2404, 0.2400, 0.2446, 0.2508, 0.2585, 0.2568, 0.2574, 0.2601, 0.2638, 0.2690, 0.2664, 0.2763, 0.2775, 0.2806, 0.2903, 0.2899, 0.2879, 0.2971, 0.3041, 0.2923, 0.3003, 0.3067, 0.3065, 0.3033, 0.2982, 0.3242, 0.0000, 0.4040, 0.4701, 0.5105, 0.5075, 0.5090, 0.5166, 0.5209, 0.4723, 0.2242, 0.2254, 0.2157, 0.2246, 0.2245, 0.2373, 0.2323, 0.2372, 0.2460, 0.2432, 0.2445, 0.2504, 0.2522, 0.2524, 0.2655, 0.2600, 0.2828, 0.2843, 0.2877, 0.2889, 0.2957, 0.2945, 0.3053, 0.3077, 0.2982, 0.3068, 0.3066, 0.3109, 0.0000, 0.4032, 0.4926, 0.5108, 0.4975, 0.4988, 0.5093, 0.5165, 0.4915, 0.2091, 0.2115, 0.1980, 0.2030, 0.2104, 0.2097, 0.2184, 0.2138, 0.2194, 0.2172, 0.2176, 0.2270, 0.2257, 0.2313, 0.2318, 0.2412, 0.2447, 0.2391, 0.2605, 0.2514, 0.2627, 0.2575, 0.2598, 0.2555, 0.2646, 0.2699, 0.2721, 0.2679, 0.0000, 0.3451, 0.4420, 0.4600, 0.4475, 0.4581, 0.4718, 0.4745, 0.4725, 0.1812, 0.1765, 0.1678, 0.1677, 0.1730, 0.1761, 0.1762, 0.1789, 0.1776, 0.1796, 0.1835, 0.1795, 0.1874, 0.1881, 0.1870, 0.1920, 0.1917, 0.1999, 0.2101, 0.2099, 0.2083, 0.2141, 0.2171, 0.2097, 0.2108, 0.2074, 0.2140, 0.2060, 0.0000, 0.2712, 0.3687, 0.4008, 0.3954, 0.4089, 0.4173, 0.4276, 0.4311, 0.1215, 0.1171, 0.1121, 0.1130, 0.1175, 0.1203, 0.1243, 0.1227, 0.1249, 0.1275, 0.1210, 0.1237, 0.1244, 0.1309, 0.1283, 0.1346, 0.1415, 0.1446, 0.1448, 0.1364, 0.1478, 0.1465, 0.1435, 0.1551, 0.1388, 0.1463, 0.1488, 0.1411, 0.0000, 0.2700, 0.3726, 0.4012, 0.3999, 0.4128, 0.4228, 0.4307, 0.4341, 0.1343, 0.1344, 0.1293, 0.1254, 0.1290, 0.1307, 0.1342, 0.1328, 0.1394, 0.1400, 0.1359, 0.1423, 0.1383, 0.1482, 0.1471, 0.1517, 0.1492, 0.1599, 0.1516, 0.1568, 0.1603, 0.1635, 0.1599, 0.1606, 0.1603, 0.1598, 0.1729, 0.1660, 0.0000, 0.3252, 0.4427, 0.4683, 0.4525, 0.4645, 0.4751, 0.4804, 0.4789, 0.1826, 0.1817, 0.1737, 0.1719, 0.1753, 0.1784, 0.1813, 0.1846, 0.1799, 0.1903, 0.1841, 0.1909, 0.1912, 0.1969, 0.1885, 0.2016, 0.2069, 0.2072, 0.1970, 0.2125, 0.2100, 0.2110, 0.2168, 0.2154, 0.2173, 0.2194, 0.2248, 0.2107, 0.0000, 0.3561, 0.4868, 0.5036, 0.4920, 0.4958, 0.5034, 0.5112, 0.4861, 0.1933, 0.1911, 0.1848, 0.1879, 0.1946, 0.1951, 0.1994, 0.2005, 0.2016, 0.2030, 0.2010, 0.2057, 0.2074, 0.2175, 0.2217, 0.2154, 0.2363, 0.2318, 0.2403, 0.2396, 0.2342, 0.2434, 0.2385, 0.2441, 0.2380, 0.2486, 0.2488, 0.2550, 0.0000, 0.3630, 0.4654, 0.5047, 0.4980, 0.5007, 0.5093, 0.5194, 0.4649, 0.2338, 0.2319, 0.2199, 0.2245, 0.2286, 0.2332, 0.2378, 0.2368, 0.2405, 0.2390, 0.2432, 0.2463, 0.2499, 0.2562, 0.2581, 0.2576, 0.2713, 0.2778, 0.2863, 0.2819, 0.2895, 0.2963, 0.2897, 0.2830, 0.2869, 0.3009, 0.2923, 0.3181, 0.0000, 0.3910, 0.4739, 0.5051, 0.5008, 0.5014, 0.5108, 0.5198, 0.4457, 0.2904, 0.2829, 0.2752, 0.2773, 0.2820, 0.2859, 0.2957, 0.2936, 0.2957, 0.2982, 0.2976, 0.3007, 0.3072, 0.3058, 0.3202, 0.3237, 0.3326, 0.3293, 0.3305, 0.3464, 0.3396, 0.3269, 0.3417, 0.3478, 0.3430, 0.3472, 0.3433, 0.3875, 0.0000, 0.4216, 0.4936, 0.5101, 0.5054, 0.5005, 0.5077, 0.5101, 0.4259, 0.3187, 0.3101, 0.2967, 0.3010, 0.3035, 0.3102, 0.3110, 0.3129, 0.3144, 0.3174, 0.3242, 0.3323, 0.3322, 0.3300, 0.3343, 0.3475, 0.3521, 0.3429, 0.3491, 0.3592, 0.3700, 0.3532, 0.3806, 0.3486, 0.3671, 0.3798, 0.3808, 0.3984, 0.0000, 0.4694, 0.5071, 0.5102, 0.5025, 0.5013, 0.5107, 0.5136, 0.3797, 0.3257, 0.3214, 0.3061, 0.3098, 0.3152, 0.3253, 0.3237, 0.3291, 0.3315, 0.3303, 0.3439, 0.3409, 0.3513, 0.3532, 0.3582, 0.3686, 0.3700, 0.3796, 0.3664, 0.3786, 0.3741, 0.3801, 0.3896, 0.3939, 0.4071, 0.3942, 0.4073, 0.4455, 0.0000, 0.5246, 0.5206, 0.5129, 0.5007, 0.5011, 0.5094, 0.5149, 0.3079, 0.2997, 0.2948, 0.2844, 0.2866, 0.2993, 0.3000, 0.3032, 0.3064, 0.3106, 0.3177, 0.3176, 0.3225, 0.3319, 0.3317, 0.3378, 0.3507, 0.3551, 0.3456, 0.3625, 0.3698, 0.3570, 0.3629, 0.3669, 0.3698, 0.3725, 0.3997, 0.3872, 0.4163, 0.0000, 0.5407, 0.5209, 0.5114, 0.5037, 0.4997, 0.5077, 0.4915, 0.2569, 0.2815, 0.2731, 0.2669, 0.2652, 0.2720, 0.2759, 0.2745, 0.2815, 0.2831, 0.2866, 0.2906, 0.2964, 0.2956, 0.3009, 0.3031, 0.3165, 0.3142, 0.3134, 0.3305, 0.3323, 0.3246, 0.3320, 0.3324, 0.3317, 0.3358, 0.3507, 0.3684, 0.3533, 0.0000, 0.5404, 0.5250, 0.5098, 0.4997, 0.4998, 0.5077, 0.4251, 0.2495, 0.2793, 0.2713, 0.2659, 0.2699, 0.2746, 0.2803, 0.2778, 0.2893, 0.2850, 0.2897, 0.2942, 0.3007, 0.3040, 0.3041, 0.3162, 0.3166, 0.3187, 0.3330, 0.3217, 0.3301, 0.3321, 0.3279, 0.3316, 0.3490, 0.3509, 0.3387, 0.3445, 0.3634, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000}, "pT rap effeciency values for Pions"};
  Configurable<std::vector<float>> effPtValuesKa{"effPtValuesKa", {0, 0, 0, 0.328845, 0.379771, 0.390088, 0.403074, 0.35504, 0.256438, 0.131726, 0.13796, 0.140295, 0.147229, 0.156968, 0.162245, 0.171312, 0.175851, 0.185823, 0.188763, 0.193965, 0.192999, 0.191121, 0.195547, 0.210082, 0.217502, 0.232456, 0.245035, 0.254051, 0.268206, 0.274664, 0.290428, 0.294979, 0.304817, 0.324206, 0.342578, 0.36466, 0.394134}, "effeciency values for Kaons"};
  Configurable<std::vector<float>> effPtRapValuesKa{"effPtRapValuesKa", {0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.1708, 0.4335, 0.4545, 0.4420, 0.1393, 0.1644, 0.1837, 0.1891, 0.1730, 0.1792, 0.1780, 0.1921, 0.2063, 0.2119, 0.2194, 0.2169, 0.2148, 0.2454, 0.2667, 0.2636, 0.2991, 0.2966, 0.3119, 0.3211, 0.3646, 0.3758, 0.3564, 0.4048, 0.3969, 0.4409, 0.5256, 0.5629, 0.4932, 0.5719, 0.0000, 0.0000, 0.0000, 0.4053, 0.4337, 0.4542, 0.4876, 0.2247, 0.1282, 0.1516, 0.1638, 0.1677, 0.1906, 0.2053, 0.2125, 0.2182, 0.2249, 0.2646, 0.2576, 0.2416, 0.2359, 0.2702, 0.2548, 0.2632, 0.2773, 0.3115, 0.3361, 0.3433, 0.3307, 0.3315, 0.3551, 0.3973, 0.4025, 0.4464, 0.4720, 0.5050, 0.5728, 0.0000, 0.0000, 0.0000, 0.3984, 0.4181, 0.4463, 0.4955, 0.3615, 0.1372, 0.1724, 0.1809, 0.1746, 0.1818, 0.1983, 0.2034, 0.2147, 0.2183, 0.2393, 0.2340, 0.2475, 0.2407, 0.2593, 0.2535, 0.2881, 0.2937, 0.2890, 0.3213, 0.3390, 0.3874, 0.4076, 0.4169, 0.3844, 0.4156, 0.4187, 0.5110, 0.5493, 0.5637, 0.0000, 0.0000, 0.0000, 0.4007, 0.4279, 0.4688, 0.4743, 0.4937, 0.1448, 0.1657, 0.1785, 0.1745, 0.1932, 0.2011, 0.2124, 0.2352, 0.2525, 0.2666, 0.2588, 0.2606, 0.2572, 0.2815, 0.2798, 0.2723, 0.2941, 0.3335, 0.3499, 0.3584, 0.3652, 0.4144, 0.3913, 0.4103, 0.4382, 0.4465, 0.4678, 0.5728, 0.5711, 0.0000, 0.0000, 0.0000, 0.4034, 0.4151, 0.4528, 0.4750, 0.5040, 0.2426, 0.1757, 0.1878, 0.1925, 0.2072, 0.2163, 0.2324, 0.2341, 0.2454, 0.2518, 0.2536, 0.2594, 0.2582, 0.2539, 0.2748, 0.2852, 0.2943, 0.3057, 0.3257, 0.3241, 0.3396, 0.3600, 0.3829, 0.4030, 0.3686, 0.3908, 0.4253, 0.4945, 0.5103, 0.0000, 0.0000, 0.0000, 0.3771, 0.4078, 0.4427, 0.4754, 0.5002, 0.3424, 0.1813, 0.1800, 0.1752, 0.1901, 0.1969, 0.1953, 0.2109, 0.2141, 0.2267, 0.2454, 0.2449, 0.2289, 0.2175, 0.2293, 0.2378, 0.2566, 0.2751, 0.2839, 0.2974, 0.3152, 0.3452, 0.3383, 0.3693, 0.3381, 0.4057, 0.3908, 0.4271, 0.4690, 0.0000, 0.0000, 0.0000, 0.3691, 0.4237, 0.4247, 0.4536, 0.4860, 0.4040, 0.1578, 0.1594, 0.1746, 0.1767, 0.1910, 0.2006, 0.2116, 0.1968, 0.2359, 0.2172, 0.2323, 0.2505, 0.2133, 0.2216, 0.2375, 0.2630, 0.2904, 0.2968, 0.3123, 0.3354, 0.3172, 0.3454, 0.3736, 0.3438, 0.3865, 0.3849, 0.4442, 0.4483, 0.0000, 0.0000, 0.0000, 0.3577, 0.4094, 0.4182, 0.4579, 0.4832, 0.4525, 0.1238, 0.1379, 0.1448, 0.1520, 0.1617, 0.1719, 0.1874, 0.1899, 0.2010, 0.2061, 0.2279, 0.2171, 0.1931, 0.2008, 0.2138, 0.2298, 0.2520, 0.2500, 0.2984, 0.3012, 0.2999, 0.3226, 0.2992, 0.3212, 0.3450, 0.3548, 0.3536, 0.3941, 0.0000, 0.0000, 0.0000, 0.3504, 0.3854, 0.3931, 0.4188, 0.4440, 0.4497, 0.1152, 0.1316, 0.1285, 0.1372, 0.1455, 0.1513, 0.1668, 0.1602, 0.1792, 0.1723, 0.1859, 0.1913, 0.1742, 0.1678, 0.1789, 0.1925, 0.2049, 0.2297, 0.2303, 0.2310, 0.2380, 0.2535, 0.2577, 0.2840, 0.3064, 0.2902, 0.3095, 0.3119, 0.0000, 0.0000, 0.0000, 0.2881, 0.3204, 0.3320, 0.3706, 0.3965, 0.4150, 0.0748, 0.0839, 0.0799, 0.0830, 0.0975, 0.0984, 0.1077, 0.1049, 0.1078, 0.1285, 0.1327, 0.1203, 0.1229, 0.1222, 0.1263, 0.1501, 0.1488, 0.1514, 0.1534, 0.1639, 0.1761, 0.1908, 0.1917, 0.1922, 0.1874, 0.2393, 0.2171, 0.2308, 0.0000, 0.0000, 0.0000, 0.2963, 0.3213, 0.3363, 0.3690, 0.3945, 0.4268, 0.0839, 0.0886, 0.0924, 0.1016, 0.1054, 0.1020, 0.1157, 0.1182, 0.1193, 0.1422, 0.1379, 0.1593, 0.1383, 0.1456, 0.1455, 0.1355, 0.1689, 0.1871, 0.1652, 0.1819, 0.1904, 0.2112, 0.2036, 0.2195, 0.2244, 0.2416, 0.2359, 0.2342, 0.0000, 0.0000, 0.0000, 0.3358, 0.3810, 0.3932, 0.4273, 0.4669, 0.4764, 0.1178, 0.1273, 0.1318, 0.1347, 0.1543, 0.1507, 0.1552, 0.1507, 0.1615, 0.1678, 0.1778, 0.1870, 0.1735, 0.1818, 0.1945, 0.2084, 0.2147, 0.2266, 0.2278, 0.2387, 0.2596, 0.2341, 0.2737, 0.2855, 0.2857, 0.2753, 0.2622, 0.3174, 0.0000, 0.0000, 0.0000, 0.3498, 0.3972, 0.4096, 0.4354, 0.4724, 0.4403, 0.1168, 0.1295, 0.1346, 0.1403, 0.1585, 0.1595, 0.1632, 0.1760, 0.1791, 0.1768, 0.2044, 0.2106, 0.1884, 0.1915, 0.2051, 0.2017, 0.2367, 0.2455, 0.2523, 0.2763, 0.2790, 0.2973, 0.2756, 0.2926, 0.3040, 0.3103, 0.3124, 0.3574, 0.0000, 0.0000, 0.0000, 0.3490, 0.4007, 0.4140, 0.4482, 0.4721, 0.3885, 0.1655, 0.1738, 0.1867, 0.1909, 0.2132, 0.2115, 0.2112, 0.2303, 0.2310, 0.2331, 0.2474, 0.2457, 0.2362, 0.2460, 0.2674, 0.2495, 0.2806, 0.3042, 0.3346, 0.3025, 0.3450, 0.3676, 0.3715, 0.3641, 0.3555, 0.4242, 0.4336, 0.4845, 0.0000, 0.0000, 0.0000, 0.3793, 0.4286, 0.4340, 0.4558, 0.4801, 0.3475, 0.1974, 0.2060, 0.2062, 0.2204, 0.2219, 0.2316, 0.2286, 0.2435, 0.2539, 0.2496, 0.2619, 0.2560, 0.2527, 0.2691, 0.2855, 0.3096, 0.2912, 0.3290, 0.3231, 0.3723, 0.3811, 0.4168, 0.4368, 0.4074, 0.4640, 0.4214, 0.4412, 0.6047, 0.0000, 0.0000, 0.0000, 0.3928, 0.4173, 0.4334, 0.4638, 0.5108, 0.2496, 0.1975, 0.2030, 0.2253, 0.2299, 0.2482, 0.2629, 0.2640, 0.2635, 0.2807, 0.3057, 0.3093, 0.2841, 0.2729, 0.2939, 0.3034, 0.3330, 0.3442, 0.3741, 0.3786, 0.3784, 0.4199, 0.4168, 0.4131, 0.4654, 0.4308, 0.4555, 0.5481, 0.5900, 0.0000, 0.0000, 0.0000, 0.4041, 0.4263, 0.4422, 0.4589, 0.4866, 0.1647, 0.1885, 0.1861, 0.1997, 0.2271, 0.2389, 0.2348, 0.2591, 0.2764, 0.2904, 0.3016, 0.3062, 0.2738, 0.2962, 0.3160, 0.3489, 0.3599, 0.3718, 0.3870, 0.4055, 0.4160, 0.4345, 0.4826, 0.4519, 0.4791, 0.4762, 0.5269, 0.6199, 0.6775, 0.0000, 0.0000, 0.0000, 0.4072, 0.4232, 0.4613, 0.4841, 0.3580, 0.1494, 0.1841, 0.1980, 0.2009, 0.1910, 0.2019, 0.2317, 0.2325, 0.2386, 0.2463, 0.2629, 0.2617, 0.2696, 0.2792, 0.2693, 0.3335, 0.2920, 0.3551, 0.3521, 0.4019, 0.4176, 0.3858, 0.4530, 0.4643, 0.4276, 0.5075, 0.5590, 0.5567, 0.6103, 0.0000, 0.0000, 0.0000, 0.4051, 0.4326, 0.4462, 0.4733, 0.2214, 0.1328, 0.1605, 0.1703, 0.1822, 0.1938, 0.2114, 0.2328, 0.2489, 0.2553, 0.2664, 0.2602, 0.2677, 0.2590, 0.2611, 0.2584, 0.3021, 0.3118, 0.3087, 0.3297, 0.3414, 0.3893, 0.3550, 0.3937, 0.4043, 0.4075, 0.4396, 0.5144, 0.5738, 0.5785, 0.0000, 0.0000, 0.0000, 0.1661, 0.4374, 0.4470, 0.4326, 0.1511, 0.1718, 0.2039, 0.1941, 0.1819, 0.1787, 0.1915, 0.1857, 0.2109, 0.2293, 0.2368, 0.2242, 0.2308, 0.2467, 0.2568, 0.2726, 0.2919, 0.3205, 0.3332, 0.3621, 0.3788, 0.3659, 0.4094, 0.4147, 0.4206, 0.4576, 0.5174, 0.5668, 0.5717, 0.6030, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000}, "pT rap effeciency values for Kaons"};
  Configurable<std::vector<float>> effPtValuesPr{"effPtValuesPr", {0, 0, 0, 0, 0, 0, 0, 0.413799, 0.443597, 0.478144, 0.505512, 0.514127, 0.523279, 0.506121, 0.481129, 0.436858, 0.365426, 0.244158, 0.246106, 0.249133, 0.251019, 0.256516, 0.263027, 0.258241, 0.260814, 0.270519, 0.272534, 0.271853, 0.274523, 0.279029, 0.279756, 0.285479, 0.292531, 0.292348, 0.294704, 0.295867, 0.289818}, "effeciency values for Kaons"};
  Configurable<std::vector<float>> effPtRapValuesPr{"effPtRapValuesPr", {0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0350, 0.3227, 0.6253, 0.6963, 0.6794, 0.3871, 0.3568, 0.3739, 0.3741, 0.3839, 0.3920, 0.3733, 0.3853, 0.3715, 0.3401, 0.3239, 0.3477, 0.3271, 0.3238, 0.3161, 0.3699, 0.3458, 0.3625, 0.3263, 0.3729, 0.3329, 0.4162, 0.3769, 0.3676, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0083, 0.2461, 0.5721, 0.6271, 0.6475, 0.6854, 0.7333, 0.6264, 0.3574, 0.3284, 0.3125, 0.3068, 0.2907, 0.2901, 0.3157, 0.3290, 0.3605, 0.3312, 0.3321, 0.3806, 0.3699, 0.3328, 0.3790, 0.4218, 0.3560, 0.3815, 0.4224, 0.4084, 0.3966, 0.4242, 0.3494, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.1002, 0.4701, 0.5444, 0.6183, 0.6271, 0.6466, 0.6825, 0.7341, 0.7063, 0.4644, 0.3044, 0.3042, 0.3036, 0.3436, 0.3403, 0.3202, 0.3508, 0.3570, 0.3379, 0.3234, 0.3496, 0.3952, 0.3592, 0.3558, 0.3744, 0.3872, 0.3675, 0.4067, 0.3756, 0.3148, 0.3804, 0.3394, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.5268, 0.5608, 0.5678, 0.6051, 0.6555, 0.6539, 0.6547, 0.6774, 0.7031, 0.7077, 0.3825, 0.3330, 0.3176, 0.3167, 0.3257, 0.3346, 0.3343, 0.3542, 0.3552, 0.3335, 0.3553, 0.3783, 0.3622, 0.3745, 0.3829, 0.4017, 0.4052, 0.4155, 0.4053, 0.4183, 0.4008, 0.4192, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.4958, 0.5522, 0.5850, 0.6070, 0.6356, 0.6566, 0.6527, 0.6986, 0.7230, 0.7258, 0.5841, 0.3389, 0.3065, 0.3507, 0.3404, 0.3470, 0.3700, 0.3559, 0.4174, 0.3907, 0.3927, 0.4110, 0.3821, 0.4087, 0.3856, 0.4179, 0.3846, 0.4074, 0.3998, 0.4417, 0.4030, 0.3701, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.5001, 0.5387, 0.5687, 0.5979, 0.6325, 0.6667, 0.6343, 0.6637, 0.6918, 0.7160, 0.7163, 0.4032, 0.3426, 0.3320, 0.3235, 0.3942, 0.3691, 0.3655, 0.3713, 0.3375, 0.3440, 0.3526, 0.3389, 0.3869, 0.3817, 0.3117, 0.3904, 0.3625, 0.3513, 0.4033, 0.3470, 0.3548, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.4804, 0.5030, 0.5471, 0.6193, 0.6334, 0.6618, 0.6741, 0.6791, 0.7283, 0.7427, 0.7476, 0.5547, 0.3485, 0.3388, 0.3286, 0.3421, 0.3439, 0.3235, 0.3386, 0.3548, 0.3266, 0.3326, 0.3649, 0.3570, 0.3457, 0.3576, 0.3624, 0.3652, 0.4041, 0.3640, 0.3371, 0.3712, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.4419, 0.5043, 0.5513, 0.5806, 0.6186, 0.6703, 0.6720, 0.6886, 0.7436, 0.7489, 0.7714, 0.6523, 0.2967, 0.2818, 0.2935, 0.2928, 0.3057, 0.3124, 0.3170, 0.3149, 0.3148, 0.3017, 0.3070, 0.3154, 0.3342, 0.3592, 0.3378, 0.3462, 0.3262, 0.3199, 0.3758, 0.3140, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.4038, 0.4921, 0.5492, 0.5645, 0.5991, 0.6307, 0.6477, 0.6584, 0.6855, 0.7222, 0.7616, 0.7143, 0.2660, 0.2750, 0.2842, 0.2531, 0.2539, 0.2648, 0.2608, 0.2451, 0.3028, 0.2926, 0.2565, 0.2665, 0.2860, 0.2966, 0.2979, 0.2980, 0.3018, 0.3010, 0.3075, 0.3206, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.3102, 0.4179, 0.4648, 0.5130, 0.5140, 0.5673, 0.6299, 0.5979, 0.6114, 0.6674, 0.7252, 0.6795, 0.1875, 0.1700, 0.1824, 0.1786, 0.1608, 0.1882, 0.1811, 0.1999, 0.2013, 0.1920, 0.1976, 0.1903, 0.1912, 0.2079, 0.2191, 0.2053, 0.1942, 0.1965, 0.2280, 0.1778, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.3163, 0.4231, 0.4802, 0.4893, 0.5333, 0.5497, 0.6038, 0.5997, 0.6288, 0.6778, 0.6772, 0.6903, 0.1937, 0.1960, 0.1922, 0.1820, 0.1821, 0.2076, 0.1988, 0.1971, 0.2013, 0.1875, 0.2262, 0.2070, 0.2248, 0.1937, 0.2380, 0.2450, 0.2245, 0.2393, 0.2302, 0.2090, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.4023, 0.4795, 0.5410, 0.5795, 0.5822, 0.6051, 0.6759, 0.6646, 0.6803, 0.6998, 0.7394, 0.7203, 0.2703, 0.2484, 0.2319, 0.2506, 0.2807, 0.2881, 0.2433, 0.2440, 0.2752, 0.2875, 0.2810, 0.2802, 0.2895, 0.2778, 0.2870, 0.2928, 0.2577, 0.2739, 0.2786, 0.2742, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.4392, 0.5116, 0.5504, 0.5738, 0.6041, 0.6422, 0.6585, 0.6790, 0.7010, 0.7327, 0.7515, 0.6178, 0.2717, 0.2596, 0.2611, 0.2460, 0.2794, 0.2702, 0.2864, 0.2747, 0.2795, 0.2828, 0.2950, 0.2705, 0.2973, 0.2855, 0.2977, 0.2903, 0.2951, 0.3171, 0.3202, 0.2941, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.4439, 0.4995, 0.5408, 0.5552, 0.6285, 0.6429, 0.6630, 0.6641, 0.6801, 0.7234, 0.7430, 0.5295, 0.3723, 0.3488, 0.3658, 0.3605, 0.3756, 0.3643, 0.3665, 0.3644, 0.3818, 0.3922, 0.3707, 0.3767, 0.3660, 0.3812, 0.3982, 0.3818, 0.4121, 0.4057, 0.3821, 0.4087, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.4802, 0.5405, 0.5537, 0.5879, 0.6195, 0.6429, 0.6356, 0.6538, 0.7105, 0.7092, 0.6975, 0.4333, 0.3931, 0.4090, 0.3939, 0.4047, 0.3970, 0.4104, 0.4154, 0.4295, 0.3995, 0.4013, 0.3880, 0.3569, 0.4243, 0.3714, 0.3962, 0.4198, 0.4361, 0.4200, 0.4040, 0.4374, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.5116, 0.5483, 0.5715, 0.5935, 0.6488, 0.6416, 0.6203, 0.6435, 0.6814, 0.6855, 0.5243, 0.3995, 0.3796, 0.3779, 0.3840, 0.4065, 0.4174, 0.4472, 0.4030, 0.4281, 0.4814, 0.4217, 0.4567, 0.4689, 0.4117, 0.4286, 0.3991, 0.4668, 0.4786, 0.4416, 0.4254, 0.5134, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.5333, 0.5601, 0.5965, 0.6047, 0.6303, 0.6673, 0.6368, 0.6791, 0.7238, 0.7025, 0.3686, 0.3411, 0.3513, 0.3422, 0.3791, 0.3428, 0.3462, 0.3873, 0.3683, 0.3951, 0.4290, 0.4181, 0.4171, 0.3960, 0.4049, 0.4350, 0.4397, 0.4206, 0.4598, 0.4500, 0.4367, 0.4547, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.1062, 0.4560, 0.5766, 0.6264, 0.6540, 0.6394, 0.6296, 0.6780, 0.6934, 0.4543, 0.3003, 0.3125, 0.3187, 0.3339, 0.3645, 0.3697, 0.3436, 0.3927, 0.3516, 0.3489, 0.3887, 0.3933, 0.3978, 0.3986, 0.3764, 0.3722, 0.4232, 0.4170, 0.3918, 0.4188, 0.4411, 0.3808, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0064, 0.2256, 0.5893, 0.6313, 0.6598, 0.6690, 0.7104, 0.6301, 0.3874, 0.3643, 0.3177, 0.2875, 0.2959, 0.2995, 0.2917, 0.3201, 0.3266, 0.3320, 0.3801, 0.3751, 0.3990, 0.3944, 0.3750, 0.3895, 0.3890, 0.3957, 0.3862, 0.4067, 0.4631, 0.3982, 0.4257, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0352, 0.3187, 0.6095, 0.6791, 0.6808, 0.4209, 0.3994, 0.3917, 0.4046, 0.3751, 0.3845, 0.3943, 0.3548, 0.3847, 0.3555, 0.3471, 0.3369, 0.3262, 0.3203, 0.3671, 0.3227, 0.3087, 0.3774, 0.3175, 0.3776, 0.3860, 0.4008, 0.4351, 0.4322, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000}, "pT rap effeciency values for Protons"};

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
    hist.add("QA/after/h3_nft0m_pt_nch", "Reco", kTHnSparseD, {axisMult, axisPt, axisMultFT0M});
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
    hist.add("QA/Pion/h2_pt_nch", "Reco", kTH2D, {{axisMult}, {axisPt}});
    hist.add("QA/Pion/h3_nft0m_pt_nch", "Reco", kTHnSparseD, {axisMult, axisPt, axisMultFT0M});
    hist.add("QA/Pion/h2_pt_nch_prof", "Reco", kTProfile, {axisMult});

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
    hist.add("Analysis/Charged/h_Mult_weighted", "Multiplicity", kTH1D, {axisMult});
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

    hist.add("Gen/Charged/h2_Nid_NidSim", "reco vs truth multiplicity", kTH2D, {{axisMultTPC}, {axisMultTPC}});
    hist.add("Gen/Charged/h_EtaTruth", "#eta ", kTH1D, {axisEta});
    hist.add("Gen/Charged/h_PhiTruth", "#phi ", kTH1D, {axisPhi});
    hist.add("Gen/Charged/h_PtTruth", "p_{T} ", kTH1D, {axisPt});
    hist.add("Gen/Charged/h2_PtTruth_Eta", "p_{T} vs #eta", kTH2D, {{axisEta}, {axisPt}});
    hist.add("Gen/Charged/h2_PtTruth_NFT0M", "p_{T} in Multiplicity Classes", kTH2D, {{axisPt}, {axisMultFT0M}});

    hist.add("Gen/Charged/h_Mult", "Multiplicity", kTH1D, {axisMult});
    hist.add("Gen/Charged/h_Mult_weighted", "Multiplicity", kTH1D, {axisMult});
    hist.add("Gen/Charged/h2_pt_nch", "Truth", kTH2D, {{axisMult}, {axisPt}});
    hist.add("Gen/Charged/h3_nft0m_pt_nch", "Truth", kTHnSparseD, {axisMult, axisPt, axisMultFT0M});
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
  void fillIdParticleQAHistos(T const& track, const std::vector<double>& ptBins, const std::vector<double>& rapBins, const std::vector<float>& effPtValues, const std::vector<float>& effPtRapValues, double rap, double nSigmaTPC, double nSigmaTOF, int nFT0M, int& N, int& NW, double& Q1, double& Q2, double& Q3, double& Q4)
  {
    double pt = track.pt();
    double weight = getCorrectedWeight(ptBins, rapBins, effPtValues, effPtRapValues, pt, rap, cfgCorrectionPtRapPID, cfgCorrectionPID);

    NW += weight;
    N++;
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
  void fillAnalysisHistos(int nTPC, int nFT0M, int N, int NW, double Q1, double Q2, double Q3, double Q4)
  {
    if (N == 0) {
      return;
    }
    double twopart1 = ((Q1 * Q1) - Q2);
    double threepart1 = ((Q1 * Q1 * Q1) - (3 * Q2 * Q1) + 2 * Q3);
    double fourpart1 = ((Q1 * Q1 * Q1 * Q1) - (6 * Q2 * Q1 * Q1) + (3 * Q2 * Q2) + (8 * Q3 * Q1) - 6 * Q4);

    hist.fill(HIST(Dire[Mode]) + HIST("h_Mult"), N);
    hist.fill(HIST(Dire[Mode]) + HIST("h_Mult_weighted"), NW);
    hist.fill(HIST(Dire[Mode]) + HIST("h_Q1"), nTPC, Q1, nFT0M);
    hist.fill(HIST(Dire[Mode]) + HIST("h_Q2"), nTPC, Q2, nFT0M);
    hist.fill(HIST(Dire[Mode]) + HIST("h_Q3"), nTPC, Q3, nFT0M);
    hist.fill(HIST(Dire[Mode]) + HIST("h_Q4"), nTPC, Q4, nFT0M);

    if (N > 1) {
      double meanPt = Q1 / static_cast<double>(NW);
      double nPair = (static_cast<double>(NW) * (static_cast<double>(NW) - 1));
      double twopart = twopart1 / nPair;
      double checkNDenoVar = (1 / std::sqrt(1 - (1 / static_cast<double>(NW))));
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
    int nChW = 0;

    int nPi = 0, nKa = 0, nPr = 0;
    int nPiW = 0, nKaW = 0, nPrW = 0;
    double ptCh = 0, q1Ch = 0, q2Ch = 0, q3Ch = 0, q4Ch = 0;
    double ptPi = 0, q1Pi = 0, q2Pi = 0, q3Pi = 0, q4Pi = 0;
    double ptPr = 0, q1Pr = 0, q2Pr = 0, q3Pr = 0, q4Pr = 0;
    double ptKa = 0, q1Ka = 0, q2Ka = 0, q3Ka = 0, q4Ka = 0;

    int nChSim = 0, nSim = 0, nFT0CSim = 0;
    int nPiSim = 0, nKaSim = 0, nPrSim = 0;
    double eta = 0, etaSim = -999, rapSim = -999;
    double ptChSim = 0, q1ChSim = 0, q2ChSim = 0, q3ChSim = 0, q4ChSim = 0;
    double ptPiSim = 0, q1PiSim = 0, q2PiSim = 0, q3PiSim = 0, q4PiSim = 0;
    double ptPrSim = 0, q1PrSim = 0, q2PrSim = 0, q3PrSim = 0, q4PrSim = 0;
    double ptKaSim = 0, q1KaSim = 0, q2KaSim = 0, q3KaSim = 0, q4KaSim = 0;

    double wghtCh = 1.0, wghtPi = 1.0, wghtKa = 1.0, wghtPr = 1.0;

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
          nChW += weight;
          nCh++;
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
            fillIdParticleQAHistos<QA_Pion>(track, ptBins, rapBins, effPtValuesPi, effPtRapValuesPi, rapPi, nSigmaTPCPi, nSigmaTOFPi, nFT0M, nPi, nPiW, q1Pi, q2Pi, q3Pi, q4Pi);
          }
        } else if (cfgSelOR == false && cfgSelAND == true) {
          if (selLowPi(track) == cfgSelLow && selHighPi(track) == cfgSelHigh) {
            fillIdParticleQAHistos<QA_Pion>(track, ptBins, rapBins, effPtValuesPi, effPtRapValuesPi, rapPi, nSigmaTPCPi, nSigmaTOFPi, nFT0M, nPi, nPiW, q1Pi, q2Pi, q3Pi, q4Pi);
          }
        }

        if (cfgSelOR == true && cfgSelAND == false) {
          if (selLowKa(track) == cfgSelLow || selHighKa(track) == cfgSelHigh) {
            fillIdParticleQAHistos<QA_Kaon>(track, ptBins, rapBins, effPtValuesKa, effPtRapValuesKa, rapKa, nSigmaTPCKa, nSigmaTOFKa, nFT0M, nKa, nKaW, q1Ka, q2Ka, q3Ka, q4Ka);
          }
        } else if (cfgSelOR == false && cfgSelAND == true) {
          if (selLowKa(track) == cfgSelLow && selHighKa(track) == cfgSelHigh) {
            fillIdParticleQAHistos<QA_Kaon>(track, ptBins, rapBins, effPtValuesKa, effPtRapValuesKa, rapKa, nSigmaTPCKa, nSigmaTOFKa, nFT0M, nKa, nKaW, q1Ka, q2Ka, q3Ka, q4Ka);
          }
        }

        if (cfgSelOR == true && cfgSelAND == false) {
          if (selLowPr(track) == cfgSelLow || selHighPr(track) == cfgSelHigh) {
            fillIdParticleQAHistos<QA_Proton>(track, ptBins, rapBins, effPtValuesPr, effPtRapValuesPr, rapPr, nSigmaTPCPr, nSigmaTOFPr, nFT0M, nPr, nPrW, q1Pr, q2Pr, q3Pr, q4Pr);
          }
        } else if (cfgSelOR == false && cfgSelAND == true) {
          if (selLowPr(track) == cfgSelLow && selHighPr(track) == cfgSelHigh) {
            fillIdParticleQAHistos<QA_Proton>(track, ptBins, rapBins, effPtValuesPr, effPtRapValuesPr, rapPr, nSigmaTPCPr, nSigmaTOFPr, nFT0M, nPr, nPrW, q1Pr, q2Pr, q3Pr, q4Pr);
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
            nChW += weight;
            nCh++;
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
              fillIdParticleQAHistos<QA_Pion>(track, ptBins, rapBins, effPtValuesPi, effPtRapValuesPi, rapPi, nSigmaTPCPi, nSigmaTOFPi, nFT0M, nPi, nPiW, q1Pi, q2Pi, q3Pi, q4Pi);
              fillPtMCHist<QA_Pion>(ptPi, eta, rapPi, nFT0M, pid, kPiPlus, kPiMinus);
            }

            if (std::abs(pid) == kKPlus && std::abs(rapKa) < 0.5 && track.pt() >= cfgCutKaPtMin) {
              ptKa = track.pt();
              fillIdParticleQAHistos<QA_Kaon>(track, ptBins, rapBins, effPtValuesKa, effPtRapValuesKa, rapKa, nSigmaTPCKa, nSigmaTOFKa, nFT0M, nKa, nKaW, q1Ka, q2Ka, q3Ka, q4Ka);
              fillPtMCHist<QA_Kaon>(ptKa, eta, rapKa, nFT0M, pid, kKPlus, kKMinus);
            }

            if (std::abs(pid) == kProton && std::abs(rapPr) < 0.5 && track.pt() >= cfgCutPrPtMin) {
              ptPr = track.pt();
              fillIdParticleQAHistos<QA_Proton>(track, ptBins, rapBins, effPtValuesPr, effPtRapValuesPr, rapPr, nSigmaTPCPr, nSigmaTOFPr, nFT0M, nPr, nPrW, q1Pr, q2Pr, q3Pr, q4Pr);
              fillPtMCHist<QA_Proton>(ptPr, eta, rapPr, nFT0M, pid, kProton, kProtonBar);
            }
          }

          if (cfgPidCut == true) {
            if (cfgSelOR == true && cfgSelAND == false) {
              if (selLowPi(track) == cfgSelLow || selHighPi(track) == cfgSelHigh) {
                ptPi = track.pt();
                fillIdParticleQAHistos<QA_Pion>(track, ptBins, rapBins, effPtValuesPi, effPtRapValuesPi, rapPi, nSigmaTPCPi, nSigmaTOFPi, nFT0M, nPi, nPiW, q1Pi, q2Pi, q3Pi, q4Pi);
                if (std::abs(pid) == kPiPlus) {
                  fillPtMCHist<QA_Pion>(ptPi, eta, rapPi, nFT0M, pid, kPiPlus, kPiMinus);
                }
              }
            } else if (cfgSelOR == false && cfgSelAND == true) {
              if (selLowPi(track) == cfgSelLow && selHighPi(track) == cfgSelHigh) {
                ptPi = track.pt();
                fillIdParticleQAHistos<QA_Pion>(track, ptBins, rapBins, effPtValuesPi, effPtRapValuesPi, rapPi, nSigmaTPCPi, nSigmaTOFPi, nFT0M, nPi, nPiW, q1Pi, q2Pi, q3Pi, q4Pi);
                if (std::abs(pid) == kPiPlus) {
                  fillPtMCHist<QA_Pion>(ptPi, eta, rapPi, nFT0M, pid, kPiPlus, kPiMinus);
                }
              }
            }

            if (cfgSelOR == true && cfgSelAND == false) {
              if (selLowKa(track) == cfgSelLow || selHighKa(track) == cfgSelHigh) {
                ptKa = track.pt();
                fillIdParticleQAHistos<QA_Kaon>(track, ptBins, rapBins, effPtValuesKa, effPtRapValuesKa, rapKa, nSigmaTPCKa, nSigmaTOFKa, nFT0M, nKa, nKaW, q1Ka, q2Ka, q3Ka, q4Ka);
                if (std::abs(pid) == kKPlus) {
                  fillPtMCHist<QA_Kaon>(ptKa, eta, rapKa, nFT0M, pid, kKPlus, kKMinus);
                }
              }
            } else if (cfgSelOR == false && cfgSelAND == true) {
              if (selLowKa(track) == cfgSelLow && selHighKa(track) == cfgSelHigh) {
                ptKa = track.pt();
                fillIdParticleQAHistos<QA_Kaon>(track, ptBins, rapBins, effPtValuesKa, effPtRapValuesKa, rapKa, nSigmaTPCKa, nSigmaTOFKa, nFT0M, nKa, nKaW, q1Ka, q2Ka, q3Ka, q4Ka);
                if (std::abs(pid) == kKPlus) {
                  fillPtMCHist<QA_Kaon>(ptKa, eta, rapKa, nFT0M, pid, kKPlus, kKMinus);
                }
              }
            }

            if (cfgSelOR == true && cfgSelAND == false) {
              if (selLowPr(track) == cfgSelLow || selHighPr(track) == cfgSelHigh) {
                ptPr = track.pt();
                fillIdParticleQAHistos<QA_Proton>(track, ptBins, rapBins, effPtValuesPr, effPtRapValuesPr, rapPr, nSigmaTPCPr, nSigmaTOFPr, nFT0M, nPr, nPrW, q1Pr, q2Pr, q3Pr, q4Pr);
                if (std::abs(pid) == kProton) {
                  fillPtMCHist<QA_Proton>(ptPr, eta, rapPr, nFT0M, pid, kProton, kProtonBar);
                }
              }
            } else if (cfgSelOR == false && cfgSelAND == true) {
              if (selLowPr(track) == cfgSelLow && selHighPr(track) == cfgSelHigh) {
                ptPr = track.pt();
                fillIdParticleQAHistos<QA_Proton>(track, ptBins, rapBins, effPtValuesPr, effPtRapValuesPr, rapPr, nSigmaTPCPr, nSigmaTOFPr, nFT0M, nPr, nPrW, q1Pr, q2Pr, q3Pr, q4Pr);
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
            wghtCh = getCorrectedWeight(ptBins, {}, effValuesCh, {}, pt, 0.0, false, cfgCorrection);
            hist.fill(HIST("QA/after/h2_pt_nch"), nCh, pt, wghtCh);
            hist.fill(HIST("QA/after/h3_nft0m_pt_nch"), nCh, pt, nFT0M, wghtCh);
            hist.fill(HIST("QA/after/h2_pt_nch_prof"), nCh, pt, wghtCh);
          }

          if (selLowPi(track) == cfgSelLow || selHighPi(track) == cfgSelHigh) {
            ptPi = track.pt();
            wghtPi = getCorrectedWeight(ptBins, rapBins, effPtValuesPi, effPtRapValuesPi, ptPi, rapPi, cfgCorrectionPtRapPID, cfgCorrectionPID);
            hist.fill(HIST("QA/Pion/h2_pt_nch"), nPi, ptPi, wghtPi);
            hist.fill(HIST("QA/Pion/h3_nft0m_pt_nch"), nPi, ptPi, nFT0M, wghtPi);
            hist.fill(HIST("QA/Pion/h2_pt_nch_prof"), nPi, ptPi, wghtPi);
          }

          if (selLowKa(track) == cfgSelLow || selHighKa(track) == cfgSelHigh) {
            ptKa = track.pt();
            wghtKa = getCorrectedWeight(ptBins, rapBins, effPtValuesKa, effPtRapValuesKa, ptKa, rapKa, cfgCorrectionPtRapPID, cfgCorrectionPID);
            hist.fill(HIST("QA/Kaon/h2_pt_nch"), nKa, ptKa, wghtKa);
            hist.fill(HIST("QA/Kaon/h3_nft0m_pt_nch"), nKa, ptKa, nFT0M, wghtKa);
            hist.fill(HIST("QA/Kaon/h2_pt_nch_prof"), nKa, ptKa, wghtKa);
          }

          if (selLowPr(track) == cfgSelLow || selHighPr(track) == cfgSelHigh) {
            ptPr = track.pt();
            wghtPr = getCorrectedWeight(ptBins, rapBins, effPtValuesPr, effPtRapValuesPr, ptPr, rapPr, cfgCorrectionPtRapPID, cfgCorrectionPID);
            hist.fill(HIST("QA/Proton/h2_pt_nch"), nPr, ptPr, wghtPr);
            hist.fill(HIST("QA/Proton/h3_nft0m_pt_nch"), nPr, ptPr, nFT0M, wghtPr);
            hist.fill(HIST("QA/Proton/h2_pt_nch_prof"), nPr, ptPr, wghtPr);
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
            hist.fill(HIST("Gen/Charged/h3_nft0m_pt_nch"), nChSim, pt, nFT0M);
            hist.fill(HIST("Gen/Charged/h2_pt_nch_prof"), nChSim, pt);
          }

          if (std::abs(mcPart.y()) >= 0.5)
            continue;

          if (std::abs(pid) == kPiPlus && mcPart.pt() >= cfgCutPiPtMin) {
            hist.fill(HIST("Gen/Pion/h2_pt_nch"), nPiSim, mcPart.pt());
            hist.fill(HIST("Gen/Pion/h3_nft0m_pt_nch"), nPiSim, mcPart.pt(), nFT0M);
            hist.fill(HIST("Gen/Pion/h2_pt_nch_prof"), nPiSim, mcPart.pt());
          }
          if (std::abs(pid) == kKPlus && mcPart.pt() >= cfgCutKaPtMin) {
            hist.fill(HIST("Gen/Kaon/h2_pt_nch"), nKaSim, mcPart.pt());
            hist.fill(HIST("Gen/Kaon/h3_nft0m_pt_nch"), nKaSim, mcPart.pt(), nFT0M);
            hist.fill(HIST("Gen/Kaon/h2_pt_nch_prof"), nKaSim, mcPart.pt());
          }
          if (std::abs(pid) == kProton && mcPart.pt() >= cfgCutPrPtMin) {
            hist.fill(HIST("Gen/Proton/h2_pt_nch"), nPrSim, mcPart.pt());
            hist.fill(HIST("Gen/Proton/h3_nft0m_pt_nch"), nPrSim, mcPart.pt(), nFT0M);
            hist.fill(HIST("Gen/Proton/h2_pt_nch_prof"), nPrSim, mcPart.pt());
          }
        }
      }

      hist.fill(HIST("Gen/h_Counts"), 2);
      hist.fill(HIST("QA/after/h_VtxZReco"), col.posZ());
      hist.fill(HIST("Gen/h_VtxZ"), col.mcCollision().posZ());

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

      if (nChSim > 0 && nCh > 0)
        hist.fill(HIST("Gen/Charged/h2_Nid_NidSim"), nChSim, nCh, wghtCh);

      if (nPi > 0 && nPiSim > 0)
        hist.fill(HIST("Gen/Pion/h2_Nid_NidSim"), nPiSim, nPi, wghtPi);

      if (nKa > 0 && nKaSim > 0)
        hist.fill(HIST("Gen/Kaon/h2_Nid_NidSim"), nKaSim, nKa, wghtKa);

      if (nPr > 0 && nPrSim > 0)
        hist.fill(HIST("Gen/Proton/h2_Nid_NidSim"), nPrSim, nPr, wghtPr);

      hist.fill(HIST("Gen/h_NTPC"), nTPC);
      hist.fill(HIST("Gen/h_NFT0C"), nFT0CSim);
      hist.fill(HIST("Gen/h2_NTPC_NFT0C"), nFT0CSim, nTPC);
      hist.fill(HIST("Gen/h2_NTPC_NFT0M"), nFT0M, nTPC);

      if (nFT0C != 0 && nFT0CSim != 0)
        hist.fill(HIST("Gen/h2_NFT0C_NFT0CSim"), nFT0CSim, nFT0C);

      fillAnalysisHistos<Gen_Charged>(nTPC, nFT0M, nChSim, nChSim, q1ChSim, q2ChSim, q3ChSim, q4ChSim);
      fillAnalysisHistos<Gen_Pion>(nTPC, nFT0M, nPiSim, nPiSim, q1PiSim, q2PiSim, q3PiSim, q4PiSim);
      fillAnalysisHistos<Gen_Kaon>(nTPC, nFT0M, nKaSim, nKaSim, q1KaSim, q2KaSim, q3KaSim, q4KaSim);
      fillAnalysisHistos<Gen_Proton>(nTPC, nFT0M, nPrSim, nPrSim, q1PrSim, q2PrSim, q3PrSim, q4PrSim);
    }

    if (nTPC > 0 && nCh > 0)
      hist.fill(HIST("QA/after/h2_NTPC_NCh"), nTPC, nCh, wghtCh);

    if (nPi > 0 && nTPC > 0)
      hist.fill(HIST("QA/after/h2_NTPC_NPi"), nTPC, nPi, wghtPi);

    if (nKa > 0 && nTPC > 0)
      hist.fill(HIST("QA/after/h2_NTPC_NKa"), nTPC, nKa, wghtKa);

    if (nPr > 0 && nTPC > 0)
      hist.fill(HIST("QA/after/h2_NTPC_NPr"), nTPC, nPr, wghtPr);

    fillAnalysisHistos<Analysis_Charged>(nTPC, nFT0M, nCh, nChW, q1Ch, q2Ch, q3Ch, q4Ch);
    fillAnalysisHistos<Analysis_Pion>(nTPC, nFT0M, nPi, nPiW, q1Pi, q2Pi, q3Pi, q4Pi);
    fillAnalysisHistos<Analysis_Kaon>(nTPC, nFT0M, nKa, nKaW, q1Ka, q2Ka, q3Ka, q4Ka);
    fillAnalysisHistos<Analysis_Proton>(nTPC, nFT0M, nPr, nPrW, q1Pr, q2Pr, q3Pr, q4Pr);
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
