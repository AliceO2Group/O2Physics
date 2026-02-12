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
/// \file   spectraTOF.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
///
/// \brief Task for the analysis of the spectra with the TOF detector.
///        Depending on the configuration it can also run on tiny tables.
///

// O2 includes

#include "PWGLF/DataModel/spectraTOF.h"

#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "PWGLF/DataModel/mcCentrality.h"
#include "PWGLF/Utils/inelGt.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/McCollisionExtra.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/StaticFor.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "TPDGCode.h"

#include <string>
#include <vector>
using namespace o2;
using namespace o2::track;
using namespace o2::framework;
using namespace o2::framework::expressions;

// Histograms
std::array<std::shared_ptr<TH3>, NpCharge> hDcaXYZ;
std::array<std::shared_ptr<TH3>, NpCharge> hDcaXYZPrm;
std::array<std::shared_ptr<TH3>, NpCharge> hDcaXYZStr;
std::array<std::shared_ptr<TH3>, NpCharge> hDcaXYZMat;
std::array<std::shared_ptr<TH2>, NpCharge> hDcaXYWrongCollisionPrm;
std::array<std::shared_ptr<TH2>, NpCharge> hDcaXYWrongCollisionStr;
std::array<std::shared_ptr<TH2>, NpCharge> hDcaXYWrongCollisionMat;
std::array<std::shared_ptr<TH2>, NpCharge> hDcaXYMC;             // DCA xy in the MC
std::array<std::shared_ptr<TH2>, NpCharge> hDcaZMC;              // DCA z in the MC
std::array<std::shared_ptr<TH2>, NpCharge> hDcaXYMCD0;           // DCA xy in the MC for particles from D0
std::array<std::shared_ptr<TH2>, NpCharge> hDcaZMCD0;            // DCA z in the MC for particles from D0
std::array<std::shared_ptr<TH2>, NpCharge> hDcaXYMCCharm;        // DCA xy in the MC for particles from charm
std::array<std::shared_ptr<TH2>, NpCharge> hdcaZMCCharm;         // DCA z in the MC for particles from charm
std::array<std::shared_ptr<TH2>, NpCharge> hDcaXYMCBeauty;       // DCA xy in the MC for particles from beauty
std::array<std::shared_ptr<TH2>, NpCharge> hDcaZMCBeauty;        // DCA z in the MC for particles from beauty
std::array<std::shared_ptr<TH2>, NpCharge> hDcaXYMCNotHF;        // DCA xy in the MC for particles from not a HF
std::array<std::shared_ptr<TH2>, NpCharge> hDcaZMCNotHF;         // DCA z in the MC for particles from not a HF
std::array<std::shared_ptr<TH2>, NpCharge> hDecayLengthStr;      // Decay Length for particles from Strange
std::array<std::shared_ptr<TH2>, NpCharge> hDecayLengthMCD0;     // Decay Length in the MC for particles from D0
std::array<std::shared_ptr<TH2>, NpCharge> hDecayLengthMCCharm;  // Decay Length in the MC for particles from charm
std::array<std::shared_ptr<TH2>, NpCharge> hDecayLengthMCBeauty; // Decay Length in the MC for particles from charm
std::array<std::shared_ptr<TH2>, NpCharge> hDecayLengthMCNotHF;  // Decay Length in the MC for particles from not a HF

std::array<std::shared_ptr<TH2>, NpCharge> hPtNumTOFMatchWithPIDSignalPrm; // Pt distribution of particles with a hit in the TOF and a compatible signal

std::array<std::array<std::shared_ptr<TH3>, NpCharge>, 3> hMCpdg_nsigmaTPC; // 2D array of nsigmaTPC histograms [Selection: pi,K,p][True PDG: 18 species]

// Spectra task
struct tofSpectra {
  struct : ConfigurableGroup {
    Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
    Configurable<int> cfgINELCut{"cfgINELCut", 0, "INEL event selection: 0 no sel, 1 INEL>0, 2 INEL>1"};
    Configurable<bool> askForCustomTVX{"askForCustomTVX", false, "Ask fro custom TVX rather than sel8"};
    Configurable<bool> removeITSROFrameBorder{"removeITSROFrameBorder", false, "Remove TF border"};
    Configurable<bool> removeNoSameBunchPileup{"removeNoSameBunchPileup", false, "Remove TF border"};
    Configurable<bool> requireIsGoodZvtxFT0vsPV{"requireIsGoodZvtxFT0vsPV", false, "Remove TF border"};
    Configurable<bool> requireIsVertexITSTPC{"requireIsVertexITSTPC", false, "Remove TF border"};
    Configurable<bool> removeNoTimeFrameBorder{"removeNoTimeFrameBorder", false, "Remove TF border"};
    Configurable<bool> requirekIsVertexTOFmatched{"requirekIsVertexTOFmatched", false, "Require requirekIsVertexTOFmatched"};
  } evselOptions;

  struct : ConfigurableGroup {
    Configurable<float> cfgCutEtaMax{"cfgCutEtaMax", 0.8f, "Max eta range for tracks"};
    Configurable<float> cfgCutNsigma{"cfgCutNsigma", 100.0f, "nsigma cut range for tracks"};
    Configurable<float> cfgCutEtaMin{"cfgCutEtaMin", -0.8f, "Min eta range for tracks"};
    Configurable<float> cfgCutY{"cfgCutY", 0.5f, "Y range for tracks"};
    Configurable<int> lastRequiredTrdCluster{"lastRequiredTrdCluster", 5, "Last cluster to require in TRD for track selection. -1 does not require any TRD cluster"};
    Configurable<bool> requireTrdOnly{"requireTrdOnly", false, "Require only tracks from TRD"};
    Configurable<bool> requireNoTrd{"requireNoTrd", false, "Require tracks without TRD"};
    Configurable<int> selectEvTime{"selectEvTime", 0, "Select event time flags; 0: any event time, 1: isEvTimeDefined, 2: IsEvTimeTOF, 3: IsEvTimeT0AC, 4: IsEvTimeTOFT0AV, 5: NOT isEvTimeDefined"};
  } trkselOptions;

  Configurable<bool> enableDcaGoodEvents{"enableDcaGoodEvents", true, "Enables the MC plots with the correct match between data and MC"};
  Configurable<bool> enablePureDCAHistogram{"enablePureDCAHistogram", false, "Enables the pure DCA histograms"};
  Configurable<bool> enableTrackCutHistograms{"enableTrackCutHistograms", true, "Enables track cut histograms, before and after the cut"};
  Configurable<bool> enableDeltaHistograms{"enableDeltaHistograms", true, "Enables the delta TPC and TOF histograms"};
  Configurable<bool> enableTPCTOFHistograms{"enableTPCTOFHistograms", true, "Enables TPC TOF histograms"};
  Configurable<bool> enableDCAxyzHistograms{"enableDCAxyzHistograms", false, "Enables DCAxyz correlation histograms"};
  Configurable<bool> enableDCAxyphiHistograms{"enableDCAxyphiHistograms", false, "Enables DCAxyphi correlation histograms"};
  Configurable<bool> enableDCAvsmotherHistograms{"enableDCAvsmotherHistograms", false, "Enables DCA vs mother histograms"};

  struct : ConfigurableGroup {
    ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0}, "Binning of the pT axis"};
    ConfigurableAxis binsEta{"binsEta", {100, -1, 1}, "Binning of the eta axis"};
    ConfigurableAxis binsnsigmaTPC{"binsnsigmaTPC", {200, -10, 10}, "Binning of the nsigmaTPC axis"};
    ConfigurableAxis binsnsigmaTOF{"binsnsigmaTOF", {200, -10, 10}, "Binning of the nsigmaTOF axis"};
    ConfigurableAxis binsdeltaTPC{"binsdeltaTPC", {500, -1000, 1000}, "Binning of the nsigmaTPC axis"};
    ConfigurableAxis binsdeltaTOF{"binsdeltaTOF", {500, -1000, 1000}, "Binning of the nsigmaTOF axis"};
    ConfigurableAxis binsDca{"binsDca", {VARIABLE_WIDTH, -3.0, -2.95, -2.9, -2.85, -2.8, -2.75, -2.7, -2.65, -2.6, -2.55, -2.5, -2.45, -2.4, -2.35, -2.3, -2.25, -2.2, -2.15, -2.1, -2.05, -2.0, -1.975, -1.95, -1.925, -1.9, -1.875, -1.85, -1.825, -1.8, -1.775, -1.75, -1.725, -1.7, -1.675, -1.65, -1.625, -1.6, -1.575, -1.55, -1.525, -1.5, -1.475, -1.45, -1.425, -1.4, -1.375, -1.35, -1.325, -1.3, -1.275, -1.25, -1.225, -1.2, -1.175, -1.15, -1.125, -1.1, -1.075, -1.05, -1.025, -1.0, -0.99, -0.98, -0.97, -0.96, -0.95, -0.94, -0.93, -0.92, -0.91, -0.9, -0.89, -0.88, -0.87, -0.86, -0.85, -0.84, -0.83, -0.82, -0.81, -0.8, -0.79, -0.78, -0.77, -0.76, -0.75, -0.74, -0.73, -0.72, -0.71, -0.7, -0.69, -0.68, -0.67, -0.66, -0.65, -0.64, -0.63, -0.62, -0.61, -0.6, -0.59, -0.58, -0.57, -0.56, -0.55, -0.54, -0.53, -0.52, -0.51, -0.5, -0.49, -0.48, -0.47, -0.46, -0.45, -0.44, -0.43, -0.42, -0.41, -0.4, -0.396, -0.392, -0.388, -0.384, -0.38, -0.376, -0.372, -0.368, -0.364, -0.36, -0.356, -0.352, -0.348, -0.344, -0.34, -0.336, -0.332, -0.328, -0.324, -0.32, -0.316, -0.312, -0.308, -0.304, -0.3, -0.296, -0.292, -0.288, -0.284, -0.28, -0.276, -0.272, -0.268, -0.264, -0.26, -0.256, -0.252, -0.248, -0.244, -0.24, -0.236, -0.232, -0.228, -0.224, -0.22, -0.216, -0.212, -0.208, -0.204, -0.2, -0.198, -0.196, -0.194, -0.192, -0.19, -0.188, -0.186, -0.184, -0.182, -0.18, -0.178, -0.176, -0.174, -0.172, -0.17, -0.168, -0.166, -0.164, -0.162, -0.16, -0.158, -0.156, -0.154, -0.152, -0.15, -0.148, -0.146, -0.144, -0.142, -0.14, -0.138, -0.136, -0.134, -0.132, -0.13, -0.128, -0.126, -0.124, -0.122, -0.12, -0.118, -0.116, -0.114, -0.112, -0.11, -0.108, -0.106, -0.104, -0.102, -0.1, -0.099, -0.098, -0.097, -0.096, -0.095, -0.094, -0.093, -0.092, -0.091, -0.09, -0.089, -0.088, -0.087, -0.086, -0.085, -0.084, -0.083, -0.082, -0.081, -0.08, -0.079, -0.078, -0.077, -0.076, -0.075, -0.074, -0.073, -0.072, -0.071, -0.07, -0.069, -0.068, -0.067, -0.066, -0.065, -0.064, -0.063, -0.062, -0.061, -0.06, -0.059, -0.058, -0.057, -0.056, -0.055, -0.054, -0.053, -0.052, -0.051, -0.05, -0.049, -0.048, -0.047, -0.046, -0.045, -0.044, -0.043, -0.042, -0.041, -0.04, -0.039, -0.038, -0.037, -0.036, -0.035, -0.034, -0.033, -0.032, -0.031, -0.03, -0.029, -0.028, -0.027, -0.026, -0.025, -0.024, -0.023, -0.022, -0.021, -0.02, -0.019, -0.018, -0.017, -0.016, -0.015, -0.014, -0.013, -0.012, -0.011, -0.01, -0.009, -0.008, -0.007, -0.006, -0.005, -0.004, -0.003, -0.002, -0.001, -0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.017, 0.018, 0.019, 0.02, 0.021, 0.022, 0.023, 0.024, 0.025, 0.026, 0.027, 0.028, 0.029, 0.03, 0.031, 0.032, 0.033, 0.034, 0.035, 0.036, 0.037, 0.038, 0.039, 0.04, 0.041, 0.042, 0.043, 0.044, 0.045, 0.046, 0.047, 0.048, 0.049, 0.05, 0.051, 0.052, 0.053, 0.054, 0.055, 0.056, 0.057, 0.058, 0.059, 0.06, 0.061, 0.062, 0.063, 0.064, 0.065, 0.066, 0.067, 0.068, 0.069, 0.07, 0.071, 0.072, 0.073, 0.074, 0.075, 0.076, 0.077, 0.078, 0.079, 0.08, 0.081, 0.082, 0.083, 0.084, 0.085, 0.086, 0.087, 0.088, 0.089, 0.09, 0.091, 0.092, 0.093, 0.094, 0.095, 0.096, 0.097, 0.098, 0.099, 0.1, 0.102, 0.104, 0.106, 0.108, 0.11, 0.112, 0.114, 0.116, 0.118, 0.12, 0.122, 0.124, 0.126, 0.128, 0.13, 0.132, 0.134, 0.136, 0.138, 0.14, 0.142, 0.144, 0.146, 0.148, 0.15, 0.152, 0.154, 0.156, 0.158, 0.16, 0.162, 0.164, 0.166, 0.168, 0.17, 0.172, 0.174, 0.176, 0.178, 0.18, 0.182, 0.184, 0.186, 0.188, 0.19, 0.192, 0.194, 0.196, 0.198, 0.2, 0.204, 0.208, 0.212, 0.216, 0.22, 0.224, 0.228, 0.232, 0.236, 0.24, 0.244, 0.248, 0.252, 0.256, 0.26, 0.264, 0.268, 0.272, 0.276, 0.28, 0.284, 0.288, 0.292, 0.296, 0.3, 0.304, 0.308, 0.312, 0.316, 0.32, 0.324, 0.328, 0.332, 0.336, 0.34, 0.344, 0.348, 0.352, 0.356, 0.36, 0.364, 0.368, 0.372, 0.376, 0.38, 0.384, 0.388, 0.392, 0.396, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.0, 1.025, 1.05, 1.075, 1.1, 1.125, 1.15, 1.175, 1.2, 1.225, 1.25, 1.275, 1.3, 1.325, 1.35, 1.375, 1.4, 1.425, 1.45, 1.475, 1.5, 1.525, 1.55, 1.575, 1.6, 1.625, 1.65, 1.675, 1.7, 1.725, 1.75, 1.775, 1.8, 1.825, 1.85, 1.875, 1.9, 1.925, 1.95, 1.975, 2.0, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.85, 2.9, 2.95, 3.0}, "Binning of DCA xy and z axis"};
    ConfigurableAxis binsMultiplicity{"binsMultiplicity", {100, 0, 100}, "Binning for multiplicity"};
    ConfigurableAxis binsOccupancy{"binsOccupancy", {3000, 0, 15000}, "Binning for occupancy"};
    ConfigurableAxis binsPercentile{"binsPercentile", {100, 0, 100}, "Binning for percentiles"};
    ConfigurableAxis binsImpactParam{"binsImpactParam", {2500, 0, 25}, "Binning for impact parameter"};
  } binsOptions;

  Configurable<int> multiplicityEstimator{"multiplicityEstimator", 0, "Flag to use a multiplicity estimator: 0 no multiplicity, 1 MultFV0M, 2 MultFT0M, 3 MultFDDM, 4 MultTracklets, 5 MultTPC, 6 MultNTracksPV, 7 MultNTracksPVeta1, 8 CentralityFT0C, 9 CentralityFT0M, 10 CentralityFV0A"};

  // Custom track cuts for the cut variation study
  TrackSelection customTrackCuts;
  Configurable<bool> kaonIsPvContrib{"kaonIsPvContrib", false, "Flag to ckeck if kaon tracks are from pv"};
  Configurable<bool> useCustomTrackCuts{"useCustomTrackCuts", false, "Flag to use custom track cuts"};
  Configurable<int> itsPattern{"itsPattern", 0, "0 = Run3ITSibAny, 1 = Run3ITSallAny, 2 = Run3ITSall7Layers, 3 = Run3ITSibTwo"};
  Configurable<bool> requireITS{"requireITS", true, "Additional cut on the ITS requirement"};
  Configurable<bool> requireTPC{"requireTPC", true, "Additional cut on the TPC requirement"};
  Configurable<bool> requireGoldenChi2{"requireGoldenChi2", true, "Additional cut on the GoldenChi2"};
  Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 70.f, "Additional cut on the minimum number of crossed rows in the TPC"};
  Configurable<float> minNCrossedRowsOverFindableClustersTPC{"minNCrossedRowsOverFindableClustersTPC", 0.8f, "Additional cut on the minimum value of the ratio between crossed rows and findable clusters in the TPC"};
  Configurable<float> maxChi2PerClusterTPC{"maxChi2PerClusterTPC", 4.f, "Additional cut on the maximum value of the chi2 per cluster in the TPC"};
  Configurable<float> minChi2PerClusterTPC{"minChi2PerClusterTPC", 0.5f, "Additional cut on the minimum value of the chi2 per cluster in the TPC"};
  Configurable<float> maxChi2PerClusterITS{"maxChi2PerClusterITS", 36.f, "Additional cut on the maximum value of the chi2 per cluster in the ITS"};
  Configurable<float> maxDcaXYFactor{"maxDcaXYFactor", 1.f, "Additional cut on the maximum value of the DCA xy (multiplicative factor)"};
  Configurable<float> maxDcaZ{"maxDcaZ", 2.f, "Additional cut on the maximum value of the DCA z"};
  Configurable<float> minTPCNClsFound{"minTPCNClsFound", 100.f, "Additional cut on the minimum value of the number of found clusters in the TPC"};
  Configurable<bool> makeTHnSparseChoice{"makeTHnSparseChoice", false, "choose if produce thnsparse"}; // RD
  Configurable<bool> enableTPCTOFvsEtaHistograms{"enableTPCTOFvsEtaHistograms", false, "choose if produce TPC tof vs Eta"};
  Configurable<bool> includeCentralityMC{"includeCentralityMC", false, "choose if include Centrality to MC"};
  Configurable<bool> isImpactParam{"isImpactParam", false, "choose if include impactparam to MC"};
  Configurable<bool> usePDGcode{"usePDGcode", false, "choose if include PDG code for MC closure test"};
  Configurable<bool> enableTPCTOFVsMult{"enableTPCTOFVsMult", false, "Produce TPC-TOF plots vs multiplicity"};
  Configurable<bool> includeCentralityToTracks{"includeCentralityToTracks", false, "choose if include Centrality to tracks"};
  Configurable<int> min_ITS_nClusters{"min_ITS_nClusters", 5, "minimum number of found ITS clusters"};

  // Histograms
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    // Standard process functions
    // Full
    if (doprocessFullEl) {
      LOG(info) << "Enabling process function processFullEl";
    }
    if (doprocessFullMu) {
      LOG(info) << "Enabling process function processFullMu";
    }
    if (doprocessFullPi) {
      LOG(info) << "Enabling process function processFullPi";
    }
    if (doprocessFullKa) {
      LOG(info) << "Enabling process function processFullKa";
    }
    if (doprocessFullPr) {
      LOG(info) << "Enabling process function processFullPr";
    }
    if (doprocessFullDe) {
      LOG(info) << "Enabling process function processFullDe";
    }
    if (doprocessFullTr) {
      LOG(info) << "Enabling process function processFullTr";
    }
    if (doprocessFullHe) {
      LOG(info) << "Enabling process function processFullHe";
    }
    if (doprocessFullAl) {
      LOG(info) << "Enabling process function processFullAl";
    }
    // LF Full
    if (doprocessLfFullEl) {
      LOG(info) << "Enabling process function processLfFullEl";
    }
    if (doprocessLfFullMu) {
      LOG(info) << "Enabling process function processLfFullMu";
    }
    if (doprocessLfFullPi) {
      LOG(info) << "Enabling process function processLfFullPi";
    }
    if (doprocessLfFullKa) {
      LOG(info) << "Enabling process function processLfFullKa";
    }
    if (doprocessLfFullPr) {
      LOG(info) << "Enabling process function processLfFullPr";
    }
    if (doprocessLfFullDe) {
      LOG(info) << "Enabling process function processLfFullDe";
    }
    if (doprocessLfFullTr) {
      LOG(info) << "Enabling process function processLfFullTr";
    }
    if (doprocessLfFullHe) {
      LOG(info) << "Enabling process function processLfFullHe";
    }
    if (doprocessLfFullAl) {
      LOG(info) << "Enabling process function processLfFullAl";
    }

    LOG(info) << "\tKaonIsPvContrib=" << kaonIsPvContrib.value;

    // Custom track cuts
    if (useCustomTrackCuts.value) {
      LOG(info) << "Using custom track cuts from values:";
      LOG(info) << "\trequireITS=" << requireITS.value;
      LOG(info) << "\trequireTPC=" << requireTPC.value;
      LOG(info) << "\trequireGoldenChi2=" << requireGoldenChi2.value;
      LOG(info) << "\tmaxChi2PerClusterTPC=" << maxChi2PerClusterTPC.value;
      LOG(info) << "\tminChi2PerClusterTPC=" << minChi2PerClusterTPC.value;
      LOG(info) << "\tminNCrossedRowsTPC=" << minNCrossedRowsTPC.value;
      LOG(info) << "\tmin_ITS_nClusters=" << min_ITS_nClusters.value;
      LOG(info) << "\tminTPCNClsFound=" << minTPCNClsFound.value;
      LOG(info) << "\tmaxChi2PerClusterITS=" << maxChi2PerClusterITS.value;
      LOG(info) << "\tmaxDcaZ=" << maxDcaZ.value;
      LOG(info) << "\tmakeTHnSparseChoice=" << makeTHnSparseChoice.value;

      customTrackCuts = getGlobalTrackSelectionRun3ITSMatch(itsPattern.value);
      LOG(info) << "Customizing track cuts:";
      customTrackCuts.SetRequireITSRefit(requireITS.value);
      customTrackCuts.SetRequireTPCRefit(requireTPC.value);
      customTrackCuts.SetMinNClustersITS(min_ITS_nClusters.value);
      customTrackCuts.SetRequireGoldenChi2(requireGoldenChi2.value);
      customTrackCuts.SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC.value);
      customTrackCuts.SetMaxChi2PerClusterITS(maxChi2PerClusterITS.value);
      customTrackCuts.SetMinNCrossedRowsTPC(minNCrossedRowsTPC.value);
      customTrackCuts.SetMinNClustersTPC(minTPCNClsFound.value);
      customTrackCuts.SetMinNCrossedRowsOverFindableClustersTPC(minNCrossedRowsOverFindableClustersTPC.value);
      customTrackCuts.SetMaxDcaXYPtDep([](float /*pt*/) { return 10000.f; }); // No DCAxy cut will be used, this is done via the member function of the task
      customTrackCuts.SetMaxDcaZ(maxDcaZ.value);
      customTrackCuts.print();
    }
    // Histograms
    const AxisSpec vtxZAxis{100, -20, 20, "Vtx_{z} (cm)"};
    const AxisSpec pAxis{binsOptions.binsPt, "#it{p} (GeV/#it{c})"};
    const AxisSpec ptAxis{binsOptions.binsPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec etaAxis{binsOptions.binsEta, "#eta"};
    const AxisSpec impParamAxis{binsOptions.binsImpactParam, "Impact parameter"};
    const AxisSpec occupancyAxis{binsOptions.binsOccupancy, "Occupancy Axis"};
    AxisSpec multAxis{binsOptions.binsMultiplicity, "Undefined multiplicity estimator"};
    switch (multiplicityEstimator) {
      case MultCodes::kNoMultiplicity: // No multiplicity
        break;
      case MultCodes::kMultFV0M: // MultFV0M
        multAxis.name = "MultFV0M";
        break;
      case MultCodes::kMultFT0M: // MultFT0M
        multAxis.name = "MultFT0M";
        break;
      case MultCodes::kMultFDDM: // MultFDDM
        multAxis.name = "MultFDDM";
        break;
      case MultCodes::kMultTracklets: // MultTracklets
        multAxis.name = "MultTracklets";
        break;
      case MultCodes::kMultTPC: // MultTPC
        multAxis.name = "MultTPC";
        break;
      case MultCodes::kMultNTracksPV: // MultNTracksPV
        multAxis.name = "MultNTracksPV";
        break;
      case MultCodes::kMultNTracksPVeta1: // MultNTracksPVeta1
        multAxis.name = "MultNTracksPVeta1";
        break;
      case MultCodes::kCentralityFT0C: // Centrality FT0C
        multAxis = {binsOptions.binsPercentile, "Centrality FT0C"};
        break;
      case MultCodes::kCentralityFT0M: // Centrality FT0M
        multAxis = {binsOptions.binsPercentile, "Centrality FT0M"};
        break;
      case MultCodes::kCentralityFV0A: // Centrality FV0A
        multAxis = {binsOptions.binsPercentile, "Centrality FV0A"};
        break;
      default:
        LOG(fatal) << "Unrecognized option for multiplicity " << multiplicityEstimator;
    }

    histos.add("event/vertexz", "", HistType::kTH1D, {vtxZAxis});
    histos.add("test_occupancy/event/vertexz", "", HistType::kTH1D, {vtxZAxis});
    auto h = histos.add<TH1>("evsel", "evsel", HistType::kTH1D, {{20, 0.5, 20.5}});
    h->GetXaxis()->SetBinLabel(1, "Events read");
    h->GetXaxis()->SetBinLabel(2, "INEL>0 (fraction)");
    h->GetXaxis()->SetBinLabel(3, "INEL>1 (fraction)");
    h->GetXaxis()->SetBinLabel(4, "Ev. sel. passed");
    h->GetXaxis()->SetBinLabel(5, "NoITSROFrameBorder");
    h->GetXaxis()->SetBinLabel(6, "NoITSROFrameBorder");
    h->GetXaxis()->SetBinLabel(7, "NoSameBunchPileup");
    h->GetXaxis()->SetBinLabel(8, "IsGoodZvtxFT0vsPV");
    h->GetXaxis()->SetBinLabel(9, "IsVertexITSTPC");
    h->GetXaxis()->SetBinLabel(10, "NoTimeFrameBorder");
    h->GetXaxis()->SetBinLabel(11, "INEL>0 (fraction)");
    h->GetXaxis()->SetBinLabel(12, "INEL>1 (fraction)");
    h->GetXaxis()->SetBinLabel(13, "posZ passed");
    h->GetXaxis()->SetBinLabel(14, evselOptions.cfgINELCut.value == 1 ? "INEL>0" : "INEL>0 (fraction)");
    h->GetXaxis()->SetBinLabel(15, evselOptions.cfgINELCut.value == 2 ? "INEL>1" : "INEL>1 (fraction)");

    h = histos.add<TH1>("tracksel", "tracksel", HistType::kTH1D, {{10, 0.5, 10.5}});
    h->GetXaxis()->SetBinLabel(1, "Tracks read");
    h->GetXaxis()->SetBinLabel(2, Form(" %.2f < #eta < %.2f ", trkselOptions.cfgCutEtaMin.value, trkselOptions.cfgCutEtaMax.value));
    h->GetXaxis()->SetBinLabel(3, "Quality passed");
    h->GetXaxis()->SetBinLabel(4, "TOF passed (partial)");

    h = histos.add<TH1>("evtime_tof", "event time selections from pidEvTimeFlags", kTH1D, {{10, -0.5, 9.5}});
    h->GetXaxis()->SetBinLabel(1, "AnyEvTime");
    h->GetXaxis()->SetBinLabel(2, "EvTimeDefined");
    h->GetXaxis()->SetBinLabel(3, "EvTimeTOF");
    h->GetXaxis()->SetBinLabel(4, "EvTimeT0AC");
    h->GetXaxis()->SetBinLabel(5, "EvTimeTOFT0AC");
    h->GetXaxis()->SetBinLabel(6, "AnyEvTime (selected)");
    h->GetXaxis()->SetBinLabel(7, "EvTimeDefined (selected)");
    h->GetXaxis()->SetBinLabel(8, "EvTimeTOF (selected)");
    h->GetXaxis()->SetBinLabel(9, "EvTimeT0AC (selected)");
    h->GetXaxis()->SetBinLabel(10, "EvTimeTOFT0AC (selected)");

    histos.add("Centrality/FV0A", "FV0A", HistType::kTH1D, {{binsOptions.binsPercentile, "Centrality FV0A"}});
    histos.add("Centrality/FT0M", "FT0M", HistType::kTH1D, {{binsOptions.binsPercentile, "Centrality FT0M"}});
    histos.add("Centrality/FT0A", "FT0A", HistType::kTH1D, {{binsOptions.binsPercentile, "Centrality FT0A"}});
    histos.add("Centrality/FT0C", "FT0C", HistType::kTH1D, {{binsOptions.binsPercentile, "Centrality FT0C"}});
    histos.add("Centrality/FDDM", "FDDM", HistType::kTH1D, {{binsOptions.binsPercentile, "Centrality FDDM"}});
    histos.add("Centrality/NTPV", "NTPV", HistType::kTH1D, {{binsOptions.binsPercentile, "Centrality NTPV"}});

    histos.add("Mult/FV0M", "MultFV0M", HistType::kTH1D, {{binsOptions.binsMultiplicity, "MultFV0M"}});
    histos.add("Mult/FT0M", "MultFT0M", HistType::kTH1D, {{binsOptions.binsMultiplicity, "MultFT0M"}});
    histos.add("Mult/FDDM", "MultFDDM", HistType::kTH1D, {{binsOptions.binsMultiplicity, "MultFDDM"}});
    if (doprocessBC) {
      histos.add("Mult/PerBC/FT0M", "FT0M", HistType::kTH1D, {{binsOptions.binsMultiplicity, "Multiplicity FT0M"}});
      histos.add("Mult/PerBC/FT0A", "FT0A", HistType::kTH1D, {{binsOptions.binsMultiplicity, "Multiplicity FT0A"}});
      histos.add("Mult/PerBC/FT0C", "FT0C", HistType::kTH1D, {{binsOptions.binsMultiplicity, "Multiplicity FT0C"}});
      histos.add("Mult/PerBC/sel8/FT0M", "FT0M", HistType::kTH1D, {{binsOptions.binsMultiplicity, "Multiplicity FT0M"}});
      histos.add("Mult/PerBC/sel8/FT0A", "FT0A", HistType::kTH1D, {{binsOptions.binsMultiplicity, "Multiplicity FT0A"}});
      histos.add("Mult/PerBC/sel8/FT0C", "FT0C", HistType::kTH1D, {{binsOptions.binsMultiplicity, "Multiplicity FT0C"}});
      histos.add("Mult/PerBC/sel8/FT0AvsFT0C", "FT0C", HistType::kTH2D, {{binsOptions.binsMultiplicity, "Multiplicity FT0A"}, {binsOptions.binsMultiplicity, "Multiplicity FT0C"}});
    }
    // histos.add("Mult/Tracklets", "MultTracklets", HistType::kTH1D, {{binsOptions.binsMultiplicity, "MultTracklets"}});
    histos.add("Mult/TPC", "MultTPC", HistType::kTH1D, {{binsOptions.binsMultiplicity, "MultTPC"}});
    histos.add("Mult/NTracksPV", "MultNTracksPV", HistType::kTH1D, {{binsOptions.binsMultiplicity, "MultNTracksPV"}});
    histos.add("Mult/NTracksPVeta1", "MultNTracksPVeta1", HistType::kTH1D, {{binsOptions.binsMultiplicity, "MultNTracksPVeta1"}});
    histos.add("test_occupancy/tpcCount", "TPC Track Count", HistType::kTH1F, {{100, 0, 10000, "Number of TPC Tracks"}});
    histos.add("test_occupancy/tofCount", "TOF Track Count", HistType::kTH1F, {{100, 0, 10000, "Number of TOF Tracks"}});

    const AxisSpec dcaXyAxis{binsOptions.binsDca, "DCA_{xy} (cm)"};
    const AxisSpec phiAxis{200, 0, 7, "#it{#varphi} (rad)"};
    const AxisSpec dcaZAxis{binsOptions.binsDca, "DCA_{z} (cm)"};
    const AxisSpec lengthAxis{100, 0, 600, "Track length (cm)"};
    const AxisSpec decayLengthAxis{100, 0, 1.0, "Decay Length (cm)"};

    if (enableTrackCutHistograms) {
      const AxisSpec chargeAxis{2, -2.f, 2.f, "Charge"};
      histos.add("track/pos/Eta", "Eta Positive tracks", HistType::kTH1D, {{binsOptions.binsEta, "#eta tracks"}});
      histos.add("track/neg/Eta", "Eta Negative tracks", HistType::kTH1D, {{binsOptions.binsEta, "#eta tracks"}});
      // its histograms
      histos.add("track/ITS/itsNCls", "number of found ITS clusters;# clusters ITS", kTH2D, {{8, -0.5, 7.5}, chargeAxis});
      histos.add("track/ITS/itsChi2NCl", "chi2 per ITS cluster;chi2 / cluster ITS", kTH2D, {{100, 0, 40}, chargeAxis});

      // tpc histograms
      histos.add("track/TPC/tpcNClsFindable", "number of findable TPC clusters;# findable clusters TPC", kTH2D, {{165, -0.5, 164.5}, chargeAxis});
      histos.add("track/TPC/tpcNClsFound", "number of found TPC clusters;# clusters TPC", kTH2D, {{165, -0.5, 164.5}, chargeAxis});
      histos.add("track/TPC/tpcNClsShared", "number of shared TPC clusters;# shared clusters TPC", kTH2D, {{165, -0.5, 164.5}, chargeAxis});
      histos.add("track/TPC/tpcCrossedRows", "number of crossed TPC rows;# crossed rows TPC", kTH2D, {{165, -0.5, 164.5}, chargeAxis});
      histos.add("track/TPC/tpcFractionSharedCls", "fraction of shared TPC clusters;fraction shared clusters TPC", kTH2D, {{100, 0., 1.}, chargeAxis});
      histos.add("track/TPC/tpcCrossedRowsOverFindableCls", "crossed TPC rows over findable clusters;crossed rows / findable clusters TPC", kTH2D, {{60, 0.7, 1.3}, chargeAxis});
      histos.add("track/TPC/tpcChi2NCl", "chi2 per cluster in TPC;chi2 / cluster TPC", kTH2D, {{100, 0, 10}, chargeAxis});
      histos.add("Vertex/histGenVtxMC", "MC generated vertex z position", HistType::kTH1F, {{400, -40., +40., "z position (cm)"}});
      histos.add("Vertex/RecoEvs/histGenVtxMC", "MC generated vertex z position", HistType::kTH1F, {{400, -40., +40., "z position (cm)"}});
      histos.add("Centrality/ImpParm", "Centrality", HistType::kTH1F, {impParamAxis});
      histos.add("Centrality/RecoEvs/ImpParm", "Centrality", HistType::kTH1F, {impParamAxis});

      histos.addClone("track/ITS/itsNCls", "track/selected/ITS/itsNCls");
      histos.addClone("track/ITS/itsChi2NCl", "track/selected/ITS/itsChi2NCl");
      histos.addClone("track/TPC/tpcNClsFindable", "track/selected/TPC/tpcNClsFindable");
      histos.addClone("track/TPC/tpcNClsFound", "track/selected/TPC/tpcNClsFound");
      histos.addClone("track/TPC/tpcNClsShared", "track/selected/TPC/tpcNClsShared");
      histos.addClone("track/TPC/tpcCrossedRows", "track/selected/TPC/tpcCrossedRows");
      histos.addClone("track/TPC/tpcFractionSharedCls", "track/selected/TPC/tpcFractionSharedCls");
      histos.addClone("track/TPC/tpcCrossedRowsOverFindableCls", "track/selected/TPC/tpcCrossedRowsOverFindableCls");
      histos.addClone("track/TPC/tpcChi2NCl", "track/selected/TPC/tpcChi2NCl");

      // trd histograms
      histos.add("track/TRD/trdSignal", "", HistType::kTH2D, {pAxis, {1000, 0, 1000, "TRD signal (a.u.)"}});
      histos.add("track/TRD/length", "", HistType::kTH1D, {lengthAxis});
      histos.add("track/TRD/lengthnotrd", "", HistType::kTH1D, {lengthAxis});
    }

    // 4 detectors
    histos.add("Data/pos/pt/its_tpc_trd_tof", "pos ITS-TPC-TRD-TOF", kTH1D, {ptAxis});
    histos.add("Data/neg/pt/its_tpc_trd_tof", "neg ITS-TPC-TRD-TOF", kTH1D, {ptAxis});

    // 3 detectors
    histos.add("Data/pos/pt/its_tpc_trd", "pos ITS-TPC-TRD", kTH1D, {ptAxis});
    histos.add("Data/neg/pt/its_tpc_trd", "neg ITS-TPC-TRD", kTH1D, {ptAxis});

    histos.add("Data/pos/pt/its_tpc_tof", "pos ITS-TPC-TOF", kTH1D, {ptAxis});
    histos.add("Data/neg/pt/its_tpc_tof", "neg ITS-TPC-TOF", kTH1D, {ptAxis});

    histos.add("Data/pos/pt/tpc_trd_tof", "pos TPC-TRD-TOF", kTH1D, {ptAxis});
    histos.add("Data/neg/pt/tpc_trd_tof", "neg TPC-TRD-TOF", kTH1D, {ptAxis});

    histos.add("Data/pos/pt/its_trd_tof", "pos ITS-TRD-TOF", kTH1D, {ptAxis});
    histos.add("Data/neg/pt/its_trd_tof", "neg ITS-TRD-TOF", kTH1D, {ptAxis});

    // 2 detectors
    histos.add("Data/pos/pt/its_tpc", "pos ITS-TPC", kTH1D, {ptAxis});
    histos.add("Data/neg/pt/its_tpc", "neg ITS-TPC", kTH1D, {ptAxis});

    histos.add("Data/pos/pt/trd_tof", "pos TRD-TOF", kTH1D, {ptAxis});
    histos.add("Data/neg/pt/trd_tof", "neg TRD-TOF", kTH1D, {ptAxis});

    histos.add("Data/pos/pt/tpc_tof", "pos TPC-TOF", kTH1D, {ptAxis});
    histos.add("Data/neg/pt/tpc_tof", "neg TPC-TOF", kTH1D, {ptAxis});

    histos.add("Data/pos/pt/its_tof", "pos ITS-TOF", kTH1D, {ptAxis});
    histos.add("Data/neg/pt/its_tof", "neg ITS-TOF", kTH1D, {ptAxis});

    histos.add("Data/pos/pt/tpc_trd", "pos TPC-TRD", kTH1D, {ptAxis});
    histos.add("Data/neg/pt/tpc_trd", "neg TPC-TRD", kTH1D, {ptAxis});

    histos.add("Data/pos/pt/its_trd", "pos ITS-TRD", kTH1D, {ptAxis});
    histos.add("Data/neg/pt/its_trd", "neg ITS-TRD", kTH1D, {ptAxis});

    // 1 detectors
    histos.add("Data/pos/pt/its", "pos ITS", kTH1D, {ptAxis});
    histos.add("Data/neg/pt/its", "neg ITS", kTH1D, {ptAxis});
    histos.add("Data/pos/pt/tpc", "pos TPC", kTH1D, {ptAxis});
    histos.add("Data/neg/pt/tpc", "neg TPC", kTH1D, {ptAxis});

    if (includeCentralityToTracks) {
      histos.add("Data/cent/pos/pt/its_tpc_tof", "pos ITS-TPC-TOF", kTH3D, {ptAxis, multAxis, occupancyAxis});
      histos.add("Data/cent/neg/pt/its_tpc_tof", "neg ITS-TPC-TOF", kTH3D, {ptAxis, multAxis, occupancyAxis});
      histos.add("Data/cent/pos/pt/its_tpc", "pos ITS-TPC", kTH3D, {ptAxis, multAxis, occupancyAxis});
      histos.add("Data/cent/neg/pt/its_tpc", "neg ITS-TPC", kTH3D, {ptAxis, multAxis, occupancyAxis});
      histos.add("Data/cent/pos/pt/its_tof", "pos ITS-TOF", kTH3D, {ptAxis, multAxis, occupancyAxis});
      histos.add("Data/cent/neg/pt/its_tof", "neg ITS-TOF", kTH3D, {ptAxis, multAxis, occupancyAxis});
    }
    const AxisSpec nsigmaTPCAxisOccupancy{binsOptions.binsnsigmaTPC, "nsigmaTPC"};
    if (doprocessMCclosure) {
      histos.add("nsigmatpc/mc_closure/pos/pi", "mc_closure dependent pion", kTHnSparseD, {ptAxis, nsigmaTPCAxisOccupancy, multAxis});
      histos.add("nsigmatpc/mc_closure/neg/pi", "mc_closure dependent pion", kTHnSparseD, {ptAxis, nsigmaTPCAxisOccupancy, multAxis});
      histos.add("nsigmatof/mc_closure/pos/pi", "mc_closure dependent pion", kTHnSparseD, {ptAxis, nsigmaTPCAxisOccupancy, multAxis});
      histos.add("nsigmatof/mc_closure/neg/pi", "mc_closure dependent pion", kTHnSparseD, {ptAxis, nsigmaTPCAxisOccupancy, multAxis});

      histos.add("nsigmatpc/mc_closure/pos/ka", "mc_closure dependent kaon", kTHnSparseD, {ptAxis, nsigmaTPCAxisOccupancy, multAxis});
      histos.add("nsigmatpc/mc_closure/neg/ka", "mc_closure dependent kaon", kTHnSparseD, {ptAxis, nsigmaTPCAxisOccupancy, multAxis});
      histos.add("nsigmatof/mc_closure/pos/ka", "mc_closure dependent kaon", kTHnSparseD, {ptAxis, nsigmaTPCAxisOccupancy, multAxis});
      histos.add("nsigmatof/mc_closure/neg/ka", "mc_closure dependent kaon", kTHnSparseD, {ptAxis, nsigmaTPCAxisOccupancy, multAxis});

      histos.add("nsigmatpc/mc_closure/pos/pr", "mc_closure dependent proton", kTHnSparseD, {ptAxis, nsigmaTPCAxisOccupancy, multAxis});
      histos.add("nsigmatpc/mc_closure/neg/pr", "mc_closure dependent proton", kTHnSparseD, {ptAxis, nsigmaTPCAxisOccupancy, multAxis});
      histos.add("nsigmatof/mc_closure/pos/pr", "mc_closure dependent proton", kTHnSparseD, {ptAxis, nsigmaTPCAxisOccupancy, multAxis});
      histos.add("nsigmatof/mc_closure/neg/pr", "mc_closure dependent proton", kTHnSparseD, {ptAxis, nsigmaTPCAxisOccupancy, multAxis});
    }

    if (doprocessOccupancy) {
      histos.add("nsigmatpc/test_occupancy/Mult_vs_Occupancy", "occuppancy vs Multiplicity", kTHnSparseD, {multAxis, occupancyAxis});
      histos.add("nsigmatpc/test_occupancy/pos/pi", "occuppancy dependent pion", kTHnSparseD, {ptAxis, nsigmaTPCAxisOccupancy, multAxis, occupancyAxis});
      histos.add("nsigmatpc/test_occupancy/neg/pi", "occuppancy dependent pion", kTHnSparseD, {ptAxis, nsigmaTPCAxisOccupancy, multAxis, occupancyAxis});
      histos.add("nsigmatpc/test_occupancy/pos/ka", "occuppancy dependent kaon", kTHnSparseD, {ptAxis, nsigmaTPCAxisOccupancy, multAxis, occupancyAxis});
      histos.add("nsigmatpc/test_occupancy/neg/ka", "occuppancy dependent kaon", kTHnSparseD, {ptAxis, nsigmaTPCAxisOccupancy, multAxis, occupancyAxis});
      histos.add("nsigmatpc/test_occupancy/pos/pr", "occuppancy dependent proton", kTHnSparseD, {ptAxis, nsigmaTPCAxisOccupancy, multAxis, occupancyAxis});
      histos.add("nsigmatpc/test_occupancy/neg/pr", "occuppancy dependent proton", kTHnSparseD, {ptAxis, nsigmaTPCAxisOccupancy, multAxis, occupancyAxis});

      histos.add("nsigmatof/test_occupancy/pos/pi", "occuppancy dependent pion", kTHnSparseD, {ptAxis, nsigmaTPCAxisOccupancy, multAxis, occupancyAxis});
      histos.add("nsigmatof/test_occupancy/neg/pi", "occuppancy dependent pion", kTHnSparseD, {ptAxis, nsigmaTPCAxisOccupancy, multAxis, occupancyAxis});
      histos.add("nsigmatof/test_occupancy/pos/ka", "occuppancy dependent kaon", kTHnSparseD, {ptAxis, nsigmaTPCAxisOccupancy, multAxis, occupancyAxis});
      histos.add("nsigmatof/test_occupancy/neg/ka", "occuppancy dependent kaon", kTHnSparseD, {ptAxis, nsigmaTPCAxisOccupancy, multAxis, occupancyAxis});
      histos.add("nsigmatof/test_occupancy/pos/pr", "occuppancy dependent proton", kTHnSparseD, {ptAxis, nsigmaTPCAxisOccupancy, multAxis, occupancyAxis});
      histos.add("nsigmatof/test_occupancy/neg/pr", "occuppancy dependent proton", kTHnSparseD, {ptAxis, nsigmaTPCAxisOccupancy, multAxis, occupancyAxis});
    }

    if (doprocessMC) {
      histos.add("MC/fake/pos", "Fake positive tracks", kTH1D, {ptAxis});
      histos.add("MC/fake/neg", "Fake negative tracks", kTH1D, {ptAxis});
      histos.add("MC/no_collision/pos", "No collision pos track", kTH1D, {ptAxis});
      histos.add("MC/no_collision/neg", "No collision neg track", kTH1D, {ptAxis});
      if (isImpactParam) {
        histos.add("MC/withPID/pi/pos/prm/pt/num", "recons. MC #pi^{+}", kTHnSparseD, {ptAxis, impParamAxis});
        histos.add("MC/withPID/pi/neg/prm/pt/num", "recons. MC #pi^{-}", kTHnSparseD, {ptAxis, impParamAxis});
        histos.add("MC/withPID/ka/pos/prm/pt/num", "recons. MC K^{+}", kTHnSparseD, {ptAxis, impParamAxis});
        histos.add("MC/withPID/ka/neg/prm/pt/num", "recons. MC K^{-}", kTHnSparseD, {ptAxis, impParamAxis});
        histos.add("MC/withPID/pr/pos/prm/pt/num", "recons. MC p", kTHnSparseD, {ptAxis, impParamAxis});
        histos.add("MC/withPID/pr/neg/prm/pt/num", "recons. MC #bar{p}", kTHnSparseD, {ptAxis, impParamAxis});
        histos.add("MC/withPID/pi/pos/prm/pt/num_str", "recons. MC #pi^{+}", kTHnSparseD, {ptAxis, impParamAxis});
        histos.add("MC/withPID/pi/neg/prm/pt/num_str", "recons. MC #pi^{-}", kTHnSparseD, {ptAxis, impParamAxis});
        histos.add("MC/withPID/ka/pos/prm/pt/num_str", "recons. MC K^{+}", kTHnSparseD, {ptAxis, impParamAxis});
        histos.add("MC/withPID/ka/neg/prm/pt/num_str", "recons. MC K^{-}", kTHnSparseD, {ptAxis, impParamAxis});
        histos.add("MC/withPID/pr/pos/prm/pt/num_str", "recons. MC p", kTHnSparseD, {ptAxis, impParamAxis});
        histos.add("MC/withPID/pr/neg/prm/pt/num_str", "recons. MC #bar{p}", kTHnSparseD, {ptAxis, impParamAxis});
        histos.add("MC/withPID/pi/pos/prm/pt/num_mat", "recons. MC #pi^{+}", kTHnSparseD, {ptAxis, impParamAxis});
        histos.add("MC/withPID/pi/neg/prm/pt/num_mat", "recons. MC #pi^{-}", kTHnSparseD, {ptAxis, impParamAxis});
        histos.add("MC/withPID/ka/pos/prm/pt/num_mat", "recons. MC K^{+}", kTHnSparseD, {ptAxis, impParamAxis});
        histos.add("MC/withPID/ka/neg/prm/pt/num_mat", "recons. MC K^{-}", kTHnSparseD, {ptAxis, impParamAxis});
        histos.add("MC/withPID/pr/pos/prm/pt/num_mat", "recons. MC p", kTHnSparseD, {ptAxis, impParamAxis});
        histos.add("MC/withPID/pr/neg/prm/pt/num_mat", "recons. MC #bar{p}", kTHnSparseD, {ptAxis, impParamAxis});
        histos.add("MC/withPID/pi/pos/prm/pt/numtof", "recons. MC #pi^{+}", kTHnSparseD, {ptAxis, impParamAxis});
        histos.add("MC/withPID/pi/neg/prm/pt/numtof", "recons. MC #pi^{-}", kTHnSparseD, {ptAxis, impParamAxis});
        histos.add("MC/withPID/ka/pos/prm/pt/numtof", "recons. MC K^{+}", kTHnSparseD, {ptAxis, impParamAxis});
        histos.add("MC/withPID/ka/neg/prm/pt/numtof", "recons. MC K^{-}", kTHnSparseD, {ptAxis, impParamAxis});
        histos.add("MC/withPID/pr/pos/prm/pt/numtof", "recons. MC p", kTHnSparseD, {ptAxis, impParamAxis});
        histos.add("MC/withPID/pr/neg/prm/pt/numtof", "recons. MC #bar{p}", kTHnSparseD, {ptAxis, impParamAxis});
        histos.add("MC/withPID/pi/pos/prm/pt/numtof_matched", "recons. MC #pi^{+}", kTHnSparseD, {ptAxis, impParamAxis});
        histos.add("MC/withPID/pi/neg/prm/pt/numtof_matched", "recons. MC #pi^{-}", kTHnSparseD, {ptAxis, impParamAxis});
        histos.add("MC/withPID/ka/pos/prm/pt/numtof_matched", "recons. MC K^{+}", kTHnSparseD, {ptAxis, impParamAxis});
        histos.add("MC/withPID/ka/neg/prm/pt/numtof_matched", "recons. MC K^{-}", kTHnSparseD, {ptAxis, impParamAxis});
        histos.add("MC/withPID/pr/pos/prm/pt/numtof_matched", "recons. MC p", kTHnSparseD, {ptAxis, impParamAxis});
        histos.add("MC/withPID/pr/neg/prm/pt/numtof_matched", "recons. MC #bar{p}", kTHnSparseD, {ptAxis, impParamAxis});
      } else {
        histos.add("MC/withPID/pi/pos/prm/pt/num", "recons. MC #pi^{+}", kTHnSparseD, {ptAxis, multAxis});
        histos.add("MC/withPID/pi/neg/prm/pt/num", "recons. MC #pi^{-}", kTHnSparseD, {ptAxis, multAxis});
        histos.add("MC/withPID/ka/pos/prm/pt/num", "recons. MC K^{+}", kTHnSparseD, {ptAxis, multAxis});
        histos.add("MC/withPID/ka/neg/prm/pt/num", "recons. MC K^{-}", kTHnSparseD, {ptAxis, multAxis});
        histos.add("MC/withPID/pr/pos/prm/pt/num", "recons. MC p", kTHnSparseD, {ptAxis, multAxis});
        histos.add("MC/withPID/pr/neg/prm/pt/num", "recons. MC #bar{p}", kTHnSparseD, {ptAxis, multAxis});
        histos.add("MC/withPID/pi/pos/prm/pt/num_str", "recons. MC #pi^{+}", kTHnSparseD, {ptAxis, multAxis});
        histos.add("MC/withPID/pi/neg/prm/pt/num_str", "recons. MC #pi^{-}", kTHnSparseD, {ptAxis, multAxis});
        histos.add("MC/withPID/ka/pos/prm/pt/num_str", "recons. MC K^{+}", kTHnSparseD, {ptAxis, multAxis});
        histos.add("MC/withPID/ka/neg/prm/pt/num_str", "recons. MC K^{-}", kTHnSparseD, {ptAxis, multAxis});
        histos.add("MC/withPID/pr/pos/prm/pt/num_str", "recons. MC p", kTHnSparseD, {ptAxis, multAxis});
        histos.add("MC/withPID/pr/neg/prm/pt/num_str", "recons. MC #bar{p}", kTHnSparseD, {ptAxis, multAxis});
        histos.add("MC/withPID/pi/pos/prm/pt/num_mat", "recons. MC #pi^{+}", kTHnSparseD, {ptAxis, multAxis});
        histos.add("MC/withPID/pi/neg/prm/pt/num_mat", "recons. MC #pi^{-}", kTHnSparseD, {ptAxis, multAxis});
        histos.add("MC/withPID/ka/pos/prm/pt/num_mat", "recons. MC K^{+}", kTHnSparseD, {ptAxis, multAxis});
        histos.add("MC/withPID/ka/neg/prm/pt/num_mat", "recons. MC K^{-}", kTHnSparseD, {ptAxis, multAxis});
        histos.add("MC/withPID/pr/pos/prm/pt/num_mat", "recons. MC p", kTHnSparseD, {ptAxis, multAxis});
        histos.add("MC/withPID/pr/neg/prm/pt/num_mat", "recons. MC #bar{p}", kTHnSparseD, {ptAxis, multAxis});
        histos.add("MC/withPID/pi/pos/prm/pt/numtof", "recons. MC #pi^{+}", kTHnSparseD, {ptAxis, multAxis});
        histos.add("MC/withPID/pi/neg/prm/pt/numtof", "recons. MC #pi^{-}", kTHnSparseD, {ptAxis, multAxis});
        histos.add("MC/withPID/ka/pos/prm/pt/numtof", "recons. MC K^{+}", kTHnSparseD, {ptAxis, multAxis});
        histos.add("MC/withPID/ka/neg/prm/pt/numtof", "recons. MC K^{-}", kTHnSparseD, {ptAxis, multAxis});
        histos.add("MC/withPID/pr/pos/prm/pt/numtof", "recons. MC p", kTHnSparseD, {ptAxis, multAxis});
        histos.add("MC/withPID/pr/neg/prm/pt/numtof", "recons. MC #bar{p}", kTHnSparseD, {ptAxis, multAxis});
        histos.add("MC/withPID/pi/pos/prm/pt/numtof_matched", "recons. MC #pi^{+}", kTHnSparseD, {ptAxis, multAxis});
        histos.add("MC/withPID/pi/neg/prm/pt/numtof_matched", "recons. MC #pi^{-}", kTHnSparseD, {ptAxis, multAxis});
        histos.add("MC/withPID/ka/pos/prm/pt/numtof_matched", "recons. MC K^{+}", kTHnSparseD, {ptAxis, multAxis});
        histos.add("MC/withPID/ka/neg/prm/pt/numtof_matched", "recons. MC K^{-}", kTHnSparseD, {ptAxis, multAxis});
        histos.add("MC/withPID/pr/pos/prm/pt/numtof_matched", "recons. MC p", kTHnSparseD, {ptAxis, multAxis});
        histos.add("MC/withPID/pr/neg/prm/pt/numtof_matched", "recons. MC #bar{p}", kTHnSparseD, {ptAxis, multAxis});
      }
      if (doprocessMCgen) {
        histos.add("MC/test/pi/pos/prm/pt/den", "generated MC #pi^{+}", kTHnSparseD, {ptAxis, impParamAxis});
        histos.add("MC/test/pi/neg/prm/pt/den", "generated MC #pi^{-}", kTHnSparseD, {ptAxis, impParamAxis});
        histos.add("MC/test/ka/pos/prm/pt/den", "generated MC K^{+}", kTHnSparseD, {ptAxis, impParamAxis});
        histos.add("MC/test/ka/neg/prm/pt/den", "generated MC K^{-}", kTHnSparseD, {ptAxis, impParamAxis});
        histos.add("MC/test/pr/pos/prm/pt/den", "generated MC p", kTHnSparseD, {ptAxis, impParamAxis});
        histos.add("MC/test/pr/neg/prm/pt/den", "generated MC #bar{p}", kTHnSparseD, {ptAxis, impParamAxis});
        if (doprocessMCgen_RecoEvs) {
          histos.add("MC/test/RecoEvs/pi/pos/prm/pt/num", "generated MC #pi^{+} from recons. events", kTHnSparseD, {ptAxis, impParamAxis});
          histos.add("MC/test/RecoEvs/pi/neg/prm/pt/num", "generated MC #pi^{-} from recons. events", kTHnSparseD, {ptAxis, impParamAxis});
          histos.add("MC/test/RecoEvs/ka/pos/prm/pt/num", "generated MC K^{+} from recons. events", kTHnSparseD, {ptAxis, impParamAxis});
          histos.add("MC/test/RecoEvs/ka/neg/prm/pt/num", "generated MC K^{-} from recons. events", kTHnSparseD, {ptAxis, impParamAxis});
          histos.add("MC/test/RecoEvs/pr/pos/prm/pt/num", "generated MC p from recons. events", kTHnSparseD, {ptAxis, impParamAxis});
          histos.add("MC/test/RecoEvs/pr/neg/prm/pt/num", "generated MC #bar{p} from recons. events", kTHnSparseD, {ptAxis, impParamAxis});
          histos.add("MC/test/RecoEvs/pi/pos/prm/pt/numtof", "generated MC #pi^{+} from recons. events", kTHnSparseD, {ptAxis, impParamAxis});
          histos.add("MC/test/RecoEvs/pi/neg/prm/pt/numtof", "generated MC #pi^{-} from recons. events", kTHnSparseD, {ptAxis, impParamAxis});
          histos.add("MC/test/RecoEvs/ka/pos/prm/pt/numtof", "generated MC K^{+} from recons. events", kTHnSparseD, {ptAxis, impParamAxis});
          histos.add("MC/test/RecoEvs/ka/neg/prm/pt/numtof", "generated MC K^{-} from recons. events", kTHnSparseD, {ptAxis, impParamAxis});
          histos.add("MC/test/RecoEvs/pr/pos/prm/pt/numtof", "generated MC p from recons. events", kTHnSparseD, {ptAxis, impParamAxis});
          histos.add("MC/test/RecoEvs/pr/neg/prm/pt/numtof", "generated MC #bar{p} from recons. events", kTHnSparseD, {ptAxis, impParamAxis});
        }
      }
      auto hh = histos.add<TH1>("MC/GenRecoCollisions", "Generated and Reconstructed MC Collisions", kTH1D, {{10, 0.5, 10.5}});
      hh->GetXaxis()->SetBinLabel(1, "Collisions generated");
      hh->GetXaxis()->SetBinLabel(2, "Collisions reconstructed");
      hh->GetXaxis()->SetBinLabel(3, "INEL>0");
      hh->GetXaxis()->SetBinLabel(4, "INEL>1");
      hh->GetXaxis()->SetBinLabel(5, "hasParticleInFT0C && hasParticleInFT0A");
      histos.add("MC/MultiplicityRecoEv", "MC multiplicity", kTH1D, {multAxis});
      histos.add("MC/Multiplicity", "MC multiplicity", kTH1D, {multAxis});
      histos.add("MC/MultiplicityMCINELgt0", "MC multiplicity", kTH1D, {multAxis});
      histos.add("MC/MultiplicityMCINELgt1", "MC multiplicity", kTH1D, {multAxis});
    }
    if (doprocessTrackMCLabels) {
      for (int par = 2; par <= 4; par++) {
        for (int i = 0; i < NpCharge; i++) {
          hMCpdg_nsigmaTPC[par - 2][i] = histos.add<TH3>(Form("test_mclabels/nsigmatpc/%s/%s/pdg_%i", (i < Np) ? "pos" : "neg", pN[par], PDGs[i % Np]), Form("True %s (%i) in %s selection", pTCharge[i], PDGs[i], (i < Np) ? pTCharge[par] : pTCharge[par + Np]), kTH3D, {ptAxis, nsigmaTPCAxisOccupancy, multAxis});
          hMCpdg_nsigmaTOF[par - 2][i] = histos.add<TH3>(Form("test_mclabels/nsigmatof/%s/%s/pdg_%i", (i < Np) ? "pos" : "neg", pN[par], PDGs[i % Np]), Form("True %s (%i) in %s selection", pTCharge[i], PDGs[i], (i < Np) ? pTCharge[par] : pTCharge[par + Np]), kTH3D, {ptAxis, nsigmaTOFAxisOccupancy, multAxis});
        }
      }
    }

    hMultiplicityvsPercentile = histos.add<TH2>("Mult/vsPercentile", "Multiplicity vs percentile", HistType::kTH2D, {{150, 0, 150}, {100, 0, 100, "Track multiplicity"}});

    for (int i = 0; i < NpCharge; i++) {
      switch (i) {
        case 0:
        case Np:
          if (doprocessFullEl == false && doprocessLfFullEl == false) {
            continue;
          }
          break;
        case 1:
        case Np + 1:
          if (doprocessFullMu == false && doprocessLfFullMu == false) {
            continue;
          }
          break;
        case 2:
        case Np + 2:
          if (doprocessFullPi == false && doprocessLfFullPi == false && doprocessDerived == false) {
            continue;
          }
          break;
        case 3:
        case Np + 3:
          if (doprocessFullKa == false && doprocessLfFullKa == false && doprocessDerived == false) {
            continue;
          }
          break;
        case 4:
        case Np + 4:
          if (doprocessFullPr == false && doprocessLfFullPr == false && doprocessDerived == false) {
            continue;
          }
          break;
        case 5:
        case Np + 5:
          if (doprocessFullDe == false && doprocessLfFullDe == false) {
            continue;
          }
          break;
        case 6:
        case Np + 6:
          if (doprocessFullTr == false && doprocessLfFullTr == false) {
            continue;
          }
          break;
        case 7:
        case Np + 7:
          if (doprocessFullHe == false && doprocessLfFullHe == false) {
            continue;
          }
          break;
        case 8:
        case Np + 8:
          if (doprocessFullAl == false && doprocessLfFullAl == false) {
            continue;
          }
          break;
      }

      const AxisSpec nsigmaTPCAxis{binsOptions.binsnsigmaTPC, Form("N_{#sigma}^{TPC}(%s)", pTCharge[i])};
      const AxisSpec nsigmaTOFAxis{binsOptions.binsnsigmaTOF, Form("N_{#sigma}^{TOF}(%s)", pTCharge[i])};
      const AxisSpec deltaTPCAxis{binsOptions.binsdeltaTPC, Form("#Delta^{TPC}(%s)", pTCharge[i])};
      const AxisSpec deltaTOFAxis{binsOptions.binsdeltaTOF, Form("#Delta^{TOF}(%s)", pTCharge[i])};
      if (multiplicityEstimator == MultCodes::kNoMultiplicity) {
        histos.add(hnsigmatof[i].data(), pTCharge[i], kTH2D, {ptAxis, nsigmaTOFAxis});
        histos.add(hnsigmatpc[i].data(), pTCharge[i], kTH2D, {ptAxis, nsigmaTPCAxis});
        if (enableDeltaHistograms) {
          histos.add(hdeltatof[i].data(), pTCharge[i], kTH2D, {ptAxis, deltaTOFAxis});
          histos.add(hdeltatpc[i].data(), pTCharge[i], kTH2D, {ptAxis, deltaTPCAxis});
        }
        if (enableTPCTOFHistograms) {
          if (enableTPCTOFvsEtaHistograms) {
            histos.add(hnsigmatpctof[i].data(), pTCharge[i], kTHnSparseD, {ptAxis, etaAxis, nsigmaTPCAxis, nsigmaTOFAxis});
          } else {
            histos.add(hnsigmatpctof[i].data(), pTCharge[i], kTH3D, {ptAxis, nsigmaTPCAxis, nsigmaTOFAxis});
          }
        }
      } else {
        if (makeTHnSparseChoice) {                                                                                  // RD
          histos.add(hnsigmatof[i].data(), pTCharge[i], kTHnSparseD, {ptAxis, nsigmaTOFAxis, multAxis, dcaXyAxis}); // RD
          histos.add(hnsigmatpc[i].data(), pTCharge[i], kTHnSparseD, {ptAxis, nsigmaTPCAxis, multAxis, dcaXyAxis}); // RD
        } else {
          histos.add(hnsigmatof[i].data(), pTCharge[i], kTH3D, {ptAxis, nsigmaTOFAxis, multAxis});
          histos.add(hnsigmatpc[i].data(), pTCharge[i], kTH3D, {ptAxis, nsigmaTPCAxis, multAxis});
        }
        if (enableDeltaHistograms) {
          histos.add(hdeltatof[i].data(), pTCharge[i], kTH3D, {ptAxis, deltaTOFAxis, multAxis});
          histos.add(hdeltatpc[i].data(), pTCharge[i], kTH3D, {ptAxis, deltaTPCAxis, multAxis});
        }
        if (enableTPCTOFHistograms) {
          if (enableTPCTOFVsMult) {
            if (enableTPCTOFvsEtaHistograms) {
              histos.add(hnsigmatpctof[i].data(), pTCharge[i], kTHnSparseD, {ptAxis, etaAxis, nsigmaTPCAxis, nsigmaTOFAxis, multAxis});
            } else {
              histos.add(hnsigmatpctof[i].data(), pTCharge[i], kTHnSparseD, {ptAxis, nsigmaTPCAxis, nsigmaTOFAxis, multAxis});
            }
          } else {
            if (enableTPCTOFvsEtaHistograms) {
              histos.add(hnsigmatpctof[i].data(), pTCharge[i], kTHnSparseD, {ptAxis, etaAxis, nsigmaTPCAxis, nsigmaTOFAxis});
            } else {
              histos.add(hnsigmatpctof[i].data(), pTCharge[i], kTH3D, {ptAxis, nsigmaTPCAxis, nsigmaTOFAxis});
            }
          }
        }
      }
      if (enableDCAxyphiHistograms) {
        histos.add(hdcaxyphi[i].data(), Form("%s -- 0.9 < #it{p}_{T} < 1.1 GeV/#it{c}", pTCharge[i]), kTH3D, {phiAxis, dcaXyAxis, dcaZAxis});
      }
      if (enableDCAxyzHistograms) {
        hDcaXYZ[i] = histos.add<TH3>(Form("dca/%s/%s", (i < Np) ? "pos" : "neg", pN[i % Np]), pTCharge[i], kTH3D, {ptAxis, dcaXyAxis, dcaZAxis});
      } else {
        histos.add(hdcaxy[i].data(), pTCharge[i], kTH2D, {ptAxis, dcaXyAxis});
        histos.add(hdcaz[i].data(), pTCharge[i], kTH2D, {ptAxis, dcaZAxis});
      }

      if (doprocessMC) {
        const std::string cpName = Form("/%s/%s", (i < Np) ? "pos" : "neg", pN[i % Np]);
        if (includeCentralityMC) {
          //*************************************RD**********************************************

          if (isImpactParam) {
            histos.add(hpt_num_prm[i].data(), pTCharge[i], kTHnSparseD, {ptAxis, impParamAxis, dcaXyAxis, occupancyAxis});
            histos.add(hpt_num_str[i].data(), pTCharge[i], kTHnSparseD, {ptAxis, impParamAxis, dcaXyAxis});
            histos.add(hpt_num_mat[i].data(), pTCharge[i], kTHnSparseD, {ptAxis, impParamAxis, dcaXyAxis});

            histos.add(hpt_numtof_prm[i].data(), pTCharge[i], kTHnSparseD, {ptAxis, impParamAxis, dcaXyAxis, occupancyAxis});
            histos.add(hpt_numtof_str[i].data(), pTCharge[i], kTHnSparseD, {ptAxis, impParamAxis, dcaXyAxis});
            histos.add(hpt_numtof_mat[i].data(), pTCharge[i], kTHnSparseD, {ptAxis, impParamAxis, dcaXyAxis});
          } else {
            histos.add(hpt_num_prm[i].data(), pTCharge[i], kTHnSparseD, {ptAxis, multAxis, dcaXyAxis, occupancyAxis});
            histos.add(hpt_num_str[i].data(), pTCharge[i], kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
            histos.add(hpt_num_mat[i].data(), pTCharge[i], kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});

            histos.add(hpt_numtof_prm[i].data(), pTCharge[i], kTHnSparseD, {ptAxis, multAxis, dcaXyAxis, occupancyAxis});
            histos.add(hpt_numtof_str[i].data(), pTCharge[i], kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
            histos.add(hpt_numtof_mat[i].data(), pTCharge[i], kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
          }

          histos.add(hpt_numtofgoodmatch_prm[i].data(), pTCharge[i], kTH3D, {ptAxis, multAxis, etaAxis});

          histos.add(hpt_den_prm[i].data(), pTCharge[i], kTHnSparseD, {ptAxis, multAxis});
          histos.add(hpt_den_str[i].data(), pTCharge[i], kTHnSparseD, {ptAxis, multAxis});
          histos.add(hpt_den_mat[i].data(), pTCharge[i], kTHnSparseD, {ptAxis, multAxis});

          histos.add(hpt_den_prm_recoev[i].data(), pTCharge[i], kTH3D, {ptAxis, multAxis, etaAxis});
          histos.add(hpt_den_prm_evsel[i].data(), pTCharge[i], kTH3D, {ptAxis, multAxis, etaAxis});
          histos.add(hpt_den_prm_goodev[i].data(), pTCharge[i], kTH3D, {ptAxis, multAxis, etaAxis});
          //***************************************************************************************
        } else {
          histos.add(hpt_num_prm[i].data(), pTCharge[i], kTH2D, {ptAxis, multAxis});
          histos.add(hpt_num_str[i].data(), pTCharge[i], kTH2D, {ptAxis, multAxis});
          histos.add(hpt_num_mat[i].data(), pTCharge[i], kTH2D, {ptAxis, multAxis});

          histos.add(hpt_numtof_prm[i].data(), pTCharge[i], kTH2D, {ptAxis, multAxis});
          histos.add(hpt_numtof_str[i].data(), pTCharge[i], kTH2D, {ptAxis, multAxis});
          histos.add(hpt_numtof_mat[i].data(), pTCharge[i], kTH2D, {ptAxis, multAxis});

          histos.add(hpt_numtofgoodmatch_prm[i].data(), pTCharge[i], kTH2D, {ptAxis, multAxis});
          hPtNumTOFMatchWithPIDSignalPrm[i] = histos.add<TH2>("MC" + cpName + "/prm/pt/numtofwithpid", pTCharge[i], kTH2D, {ptAxis, multAxis});

          histos.add(hpt_den_prm[i].data(), pTCharge[i], kTH2D, {ptAxis, multAxis});
          histos.add(hpt_den_str[i].data(), pTCharge[i], kTH2D, {ptAxis, multAxis});
          histos.add(hpt_den_mat[i].data(), pTCharge[i], kTH2D, {ptAxis, multAxis});

          histos.add(hpt_den_prm_recoev[i].data(), pTCharge[i], kTH2D, {ptAxis, multAxis});
          histos.add(hpt_den_prm_evsel[i].data(), pTCharge[i], kTH2D, {ptAxis, multAxis});
          histos.add(hpt_den_prm_goodev[i].data(), pTCharge[i], kTH2D, {ptAxis, multAxis});
        }
        histos.add(hpt_den_prm_mcgoodev[i].data(), pTCharge[i], kTH2D, {ptAxis, multAxis});
        histos.add(hpt_den_prm_mcbadev[i].data(), pTCharge[i], kTH2D, {ptAxis, multAxis});
        if (enableDCAxyzHistograms) {
          hDcaXYZPrm[i] = histos.add<TH3>("dcaprm" + cpName, pTCharge[i], kTH3D, {ptAxis, dcaXyAxis, dcaZAxis});
          hDcaXYZStr[i] = histos.add<TH3>("dcastr" + cpName, pTCharge[i], kTH3D, {ptAxis, dcaXyAxis, dcaZAxis});
          hDcaXYZMat[i] = histos.add<TH3>("dcamat" + cpName, pTCharge[i], kTH3D, {ptAxis, dcaXyAxis, dcaZAxis});
          if (enableDcaGoodEvents) {
            histos.add(hdcaxyprmgoodevs[i].data(), pTCharge[i], kTH3D, {ptAxis, dcaXyAxis, dcaZAxis});
          }
        } else {
          hDcaXYWrongCollisionPrm[i] = histos.add<TH2>("dcaxywrongcollprm" + cpName, pTCharge[i], kTH2D, {ptAxis, dcaXyAxis});
          hDcaXYWrongCollisionStr[i] = histos.add<TH2>("dcaxywrongcollstr" + cpName, pTCharge[i], kTH2D, {ptAxis, dcaXyAxis});
          hDcaXYWrongCollisionMat[i] = histos.add<TH2>("dcaxywrongcollmat" + cpName, pTCharge[i], kTH2D, {ptAxis, dcaXyAxis});

          histos.add(hdcaxyprm[i].data(), pTCharge[i], kTH2D, {ptAxis, dcaXyAxis});
          histos.add(hdcazprm[i].data(), pTCharge[i], kTH2D, {ptAxis, dcaZAxis});
          histos.add(hdcaxystr[i].data(), pTCharge[i], kTH2D, {ptAxis, dcaXyAxis});
          histos.add(hdcazstr[i].data(), pTCharge[i], kTH2D, {ptAxis, dcaZAxis});
          histos.add(hdcaxymat[i].data(), pTCharge[i], kTH2D, {ptAxis, dcaXyAxis});
          histos.add(hdcazmat[i].data(), pTCharge[i], kTH2D, {ptAxis, dcaZAxis});
          hDecayLengthStr[i] = histos.add<TH2>("decaylengthstr" + cpName, pTCharge[i], kTH2D, {ptAxis, decayLengthAxis});
          if (enableDcaGoodEvents) {
            histos.add(hdcaxyprmgoodevs[i].data(), pTCharge[i], kTH2D, {ptAxis, dcaXyAxis});
            histos.add(hdcazprmgoodevs[i].data(), pTCharge[i], kTH2D, {ptAxis, dcaZAxis});
          }
          if (enableDCAvsmotherHistograms) {
            hDcaXYMC[i] = histos.add<TH2>("dcaxymc" + cpName, pTCharge[i], kTH2D, {ptAxis, dcaXyAxis});
            hDcaZMC[i] = histos.add<TH2>("dcazmc" + cpName, pTCharge[i], kTH2D, {ptAxis, dcaZAxis});
            hDcaXYMCNotHF[i] = histos.add<TH2>("dcaxynothf" + cpName, pTCharge[i], kTH2D, {ptAxis, dcaXyAxis});
            hDcaZMCNotHF[i] = histos.add<TH2>("dcaznothf" + cpName, pTCharge[i], kTH2D, {ptAxis, dcaZAxis});
            hDcaXYMCD0[i] = histos.add<TH2>("dcaxyD0" + cpName, pTCharge[i], kTH2D, {ptAxis, dcaXyAxis});
            hDcaZMCD0[i] = histos.add<TH2>("dcazD0" + cpName, pTCharge[i], kTH2D, {ptAxis, dcaZAxis});
            hDcaXYMCCharm[i] = histos.add<TH2>("dcaxycharm" + cpName, pTCharge[i], kTH2D, {ptAxis, dcaXyAxis});
            hdcaZMCCharm[i] = histos.add<TH2>("dcazcharm" + cpName, pTCharge[i], kTH2D, {ptAxis, dcaZAxis});
            hDcaXYMCBeauty[i] = histos.add<TH2>("dcaxybeauty" + cpName, pTCharge[i], kTH2D, {ptAxis, dcaXyAxis});
            hDcaZMCBeauty[i] = histos.add<TH2>("dcazbeauty" + cpName, pTCharge[i], kTH2D, {ptAxis, dcaZAxis});
            hDecayLengthMCD0[i] = histos.add<TH2>("decaylengthD0" + cpName, pTCharge[i], kTH2D, {ptAxis, decayLengthAxis});
            hDecayLengthMCCharm[i] = histos.add<TH2>("decaylengthcharm" + cpName, pTCharge[i], kTH2D, {ptAxis, decayLengthAxis});
            hDecayLengthMCBeauty[i] = histos.add<TH2>("decaylengthbeauty" + cpName, pTCharge[i], kTH2D, {ptAxis, decayLengthAxis});
            hDecayLengthMCNotHF[i] = histos.add<TH2>("decaylengthnothf" + cpName, pTCharge[i], kTH2D, {ptAxis, decayLengthAxis});
          }
        }

        // Mismatched info
        histos.add(hpt_mism_its_prm[i].data(), pTCharge[i], kTH1D, {ptAxis});
        histos.add(hpt_mism_tpc_prm[i].data(), pTCharge[i], kTH1D, {ptAxis});
        histos.add(hpt_mism_trd_prm[i].data(), pTCharge[i], kTH1D, {ptAxis});
        histos.add(hpt_mism_tof_prm[i].data(), pTCharge[i], kTH1D, {ptAxis});
      }
    }
    // Print output histograms statistics
    LOG(info) << "Size of the histograms in spectraTOF";
    histos.print();
  }

  void processBC(soa::Join<aod::BCs, aod::BcSels, aod::Run3MatchedToBCSparse>::iterator const& bc,
                 aod::FT0s const&)
  {
    if (!bc.has_ft0()) {
      return;
    }
    const bool sel8 = bc.selection_bit(aod::evsel::kIsTriggerTVX) && bc.selection_bit(aod::evsel::kNoTimeFrameBorder) && bc.selection_bit(aod::evsel::kNoITSROFrameBorder);
    const auto& ft0 = bc.ft0();

    histos.fill(HIST("Mult/PerBC/FT0M"), ft0.sumAmpA() + ft0.sumAmpC());
    histos.fill(HIST("Mult/PerBC/FT0A"), ft0.sumAmpA());
    histos.fill(HIST("Mult/PerBC/FT0C"), ft0.sumAmpC());

    if (!sel8) {
      return;
    }

    histos.fill(HIST("Mult/PerBC/sel8/FT0M"), ft0.sumAmpA() + ft0.sumAmpC());
    histos.fill(HIST("Mult/PerBC/sel8/FT0A"), ft0.sumAmpA());
    histos.fill(HIST("Mult/PerBC/sel8/FT0C"), ft0.sumAmpC());

    histos.fill(HIST("Mult/PerBC/sel8/FT0AvsFT0C"), ft0.sumAmpA(), ft0.sumAmpC());

  } // end of the process function
  PROCESS_SWITCH(tofSpectra, processBC, "Processor of BCs for the FT0 calibration", true);

  template <bool fillFullInfo, PID::ID id, typename T, typename C>
  void fillParticleHistos(const T& track, const C& collision)
  {
    if (std::abs(track.rapidity(PID::getMass(id))) > trkselOptions.cfgCutY) {
      return;
    }
    if constexpr (id == PID::Kaon) {
      if (kaonIsPvContrib && !track.isPVContributor()) {
        return;
      }
    }
    const auto& nsigmaTOF = o2::aod::pidutils::tofNSigma<id>(track);
    const auto& nsigmaTPC = o2::aod::pidutils::tpcNSigma<id>(track);

    // const auto id = track.sign() > 0 ? id : id + Np;
    const float multiplicity = getMultiplicity(collision);
    if (multiplicityEstimator == MultCodes::kNoMultiplicity) {
      if (track.sign() > 0) {
        histos.fill(HIST(hnsigmatpc[id]), track.pt(), nsigmaTPC);
      } else {
        histos.fill(HIST(hnsigmatpc[id + Np]), track.pt(), nsigmaTPC);
      }
    } else if (makeTHnSparseChoice) {                                                               // RD
      if (track.sign() > 0) {                                                                       // RD
        histos.fill(HIST(hnsigmatpc[id]), track.pt(), nsigmaTPC, multiplicity, track.dcaXY());      // RD
      } else {                                                                                      // RD
        histos.fill(HIST(hnsigmatpc[id + Np]), track.pt(), nsigmaTPC, multiplicity, track.dcaXY()); // RD
      }
    } else {
      if (track.sign() > 0) {
        histos.fill(HIST(hnsigmatpc[id]), track.pt(), nsigmaTPC, multiplicity);
      } else {
        histos.fill(HIST(hnsigmatpc[id + Np]), track.pt(), nsigmaTPC, multiplicity);
      }
    }

    if constexpr (fillFullInfo) {
      if (enableDeltaHistograms) {
        const auto& deltaTPC = o2::aod::pidutils::tpcExpSignalDiff<id>(track);
        if (multiplicityEstimator == MultCodes::kNoMultiplicity) {
          if (track.sign() > 0) {
            histos.fill(HIST(hdeltatpc[id]), track.pt(), deltaTPC);
          } else {
            histos.fill(HIST(hdeltatpc[id + Np]), track.pt(), deltaTPC);
          }
        } else {
          if (track.sign() > 0) {
            histos.fill(HIST(hdeltatpc[id]), track.pt(), deltaTPC, multiplicity);
          } else {
            histos.fill(HIST(hdeltatpc[id + Np]), track.pt(), deltaTPC, multiplicity);
          }
        }
      }
    }

    // TOF part
    if (!track.hasTOF()) {
      return;
    }
    if (trkselOptions.requireTrdOnly == true && !track.hasTRD()) {
      return;
    }
    if (trkselOptions.requireNoTrd == true && track.hasTRD()) {
      return;
    }
    histos.fill(HIST("evtime_tof"), 0.f);
    if (track.isEvTimeDefined()) {
      histos.fill(HIST("evtime_tof"), 1.f);
    }
    if (track.isEvTimeTOF()) {
      histos.fill(HIST("evtime_tof"), 2.f);
    }
    if (track.isEvTimeT0AC()) {
      histos.fill(HIST("evtime_tof"), 3.f);
    }
    if (track.isEvTimeTOFT0AC()) {
      histos.fill(HIST("evtime_tof"), 4.f);
    }
    switch (trkselOptions.selectEvTime) {
      case 0:
        break;
      case 1:
        if (!track.isEvTimeDefined()) {
          return;
        }
        break;
      case 2:
        if (!track.isEvTimeTOF()) {
          return;
        }
        break;
      case 3:
        if (!track.isEvTimeT0AC()) {
          return;
        }
        break;
      case 4:
        if (!track.isEvTimeTOFT0AC()) {
          return;
        }
        break;
      case 5:
        if (track.isEvTimeDefined()) {
          return;
        }
        break;
      default:
        LOG(fatal) << "Fatal did not recognise value select event time" << trkselOptions.selectEvTime;
    }
    histos.fill(HIST("evtime_tof"), 5.f);
    if (track.isEvTimeDefined()) {
      histos.fill(HIST("evtime_tof"), 6.f);
    }
    if (track.isEvTimeTOF()) {
      histos.fill(HIST("evtime_tof"), 7.f);
    }
    if (track.isEvTimeT0AC()) {
      histos.fill(HIST("evtime_tof"), 8.f);
    }
    if (track.isEvTimeTOFT0AC()) {
      histos.fill(HIST("evtime_tof"), 9.f);
    }

    if (track.hasTRD() && (trkselOptions.lastRequiredTrdCluster > 0)) {
      int lastLayer = 0;
      for (int l = 7; l >= 0; l--) {
        if (track.trdPattern() & (1 << l)) {
          lastLayer = l;
          break;
        }
      }
      if (lastLayer < trkselOptions.lastRequiredTrdCluster) {
        return;
      }
    }

    if (multiplicityEstimator == MultCodes::kNoMultiplicity) {
      if (track.sign() > 0) {
        histos.fill(HIST(hnsigmatof[id]), track.pt(), nsigmaTOF);
      } else {
        histos.fill(HIST(hnsigmatof[id + Np]), track.pt(), nsigmaTOF);
      }
    } else {
      if (makeTHnSparseChoice) {                                                                      // RD
        if (track.sign() > 0) {                                                                       // RD
          histos.fill(HIST(hnsigmatof[id]), track.pt(), nsigmaTOF, multiplicity, track.dcaXY());      // RD
        } else {                                                                                      // RD
          histos.fill(HIST(hnsigmatof[id + Np]), track.pt(), nsigmaTOF, multiplicity, track.dcaXY()); // RD
        }
      } else {
        if (track.sign() > 0) {
          histos.fill(HIST(hnsigmatof[id]), track.pt(), nsigmaTOF, multiplicity);
        } else {
          histos.fill(HIST(hnsigmatof[id + Np]), track.pt(), nsigmaTOF, multiplicity);
        }
      }
    }

    if (enableTPCTOFHistograms) {
      if (enableTPCTOFVsMult) {
        if (enableTPCTOFvsEtaHistograms) {
          if (track.sign() > 0) {
            histos.fill(HIST(hnsigmatpctof[id]), track.pt(), track.eta(), nsigmaTPC, nsigmaTOF, multiplicity);
          } else {
            histos.fill(HIST(hnsigmatpctof[id + Np]), track.pt(), track.eta(), nsigmaTPC, nsigmaTOF, multiplicity);
          }
        } else {
          if (track.sign() > 0) {
            histos.fill(HIST(hnsigmatpctof[id]), track.pt(), nsigmaTPC, nsigmaTOF, multiplicity);
          } else {
            histos.fill(HIST(hnsigmatpctof[id + Np]), track.pt(), nsigmaTPC, nsigmaTOF, multiplicity);
          }
        }
      } else {
        if (enableTPCTOFvsEtaHistograms) {
          if (track.sign() > 0) {
            histos.fill(HIST(hnsigmatpctof[id]), track.pt(), track.eta(), nsigmaTPC, nsigmaTOF);
          } else {
            histos.fill(HIST(hnsigmatpctof[id + Np]), track.pt(), track.eta(), nsigmaTPC, nsigmaTOF);
          }
        } else {
          if (track.sign() > 0) {
            histos.fill(HIST(hnsigmatpctof[id]), track.pt(), nsigmaTPC, nsigmaTOF);
          } else {
            histos.fill(HIST(hnsigmatpctof[id + Np]), track.pt(), nsigmaTPC, nsigmaTOF);
          }
        }
      }
    }
    if constexpr (fillFullInfo) {
      if (enableDeltaHistograms) {
        const auto& deltaTOF = o2::aod::pidutils::tofExpSignalDiff<id>(track);
        if (multiplicityEstimator == MultCodes::kNoMultiplicity) {
          if (track.sign() > 0) {
            histos.fill(HIST(hdeltatof[id]), track.pt(), deltaTOF);
          } else {
            histos.fill(HIST(hdeltatof[id + Np]), track.pt(), deltaTOF);
          }
        } else {
          if (track.sign() > 0) {
            histos.fill(HIST(hdeltatof[id]), track.pt(), deltaTOF, multiplicity);
          } else {
            histos.fill(HIST(hdeltatof[id + Np]), track.pt(), deltaTOF, multiplicity);
          }
        }
      }
    }
    // Filling DCA info with the TPC+TOF PID
    bool isDCAPureSample = (std::sqrt(nsigmaTOF * nsigmaTOF + nsigmaTPC * nsigmaTPC) < 2.f);
    if (track.pt() <= 0.4) {
      isDCAPureSample = (nsigmaTPC < 1.f);
    }
    if (isDCAPureSample) {
      const bool isInPtRangeForPhi = track.pt() < 1.1f && track.pt() > 0.9f;
      if (enableDCAxyzHistograms) {
        if (track.sign() > 0) {
          hDcaXYZ[id]->Fill(track.pt(), track.dcaXY(), track.dcaZ());
          if (isInPtRangeForPhi) {
            if (enableDCAxyphiHistograms) {
              histos.fill(HIST(hdcaxyphi[id]), track.phi(), track.dcaXY(), track.dcaZ());
            }
          }
        } else {
          hDcaXYZ[id + Np]->Fill(track.pt(), track.dcaXY(), track.dcaZ());
          if (isInPtRangeForPhi) {
            if (enableDCAxyphiHistograms) {
              histos.fill(HIST(hdcaxyphi[id + Np]), track.phi(), track.dcaXY(), track.dcaZ());
            }
          }
        }
      } else {
        if (track.sign() > 0) {
          histos.fill(HIST(hdcaxy[id]), track.pt(), track.dcaXY());
          histos.fill(HIST(hdcaz[id]), track.pt(), track.dcaZ());
          if (isInPtRangeForPhi) {
            if (enableDCAxyphiHistograms) {
              histos.fill(HIST(hdcaxyphi[id]), track.phi(), track.dcaXY());
            }
          }
        } else {
          histos.fill(HIST(hdcaxy[id + Np]), track.pt(), track.dcaXY());
          histos.fill(HIST(hdcaz[id + Np]), track.pt(), track.dcaZ());
          if (isInPtRangeForPhi) {
            if (enableDCAxyphiHistograms) {
              histos.fill(HIST(hdcaxyphi[id + Np]), track.phi(), track.dcaXY());
            }
          }
        }
      }
    }
    if (!passesDCAxyCut(track)) {
      return;
    }

    if constexpr (fillFullInfo) {
    }
  }

  template <bool fillHistograms = false, bool fillMultiplicity = false, typename CollisionType>
  bool isEventSelected(CollisionType const& collision)
  {
    if constexpr (fillHistograms) {
      histos.fill(HIST("evsel"), 1.f);
    }
    if constexpr (fillHistograms) {
      if (collision.isInelGt0()) {
        histos.fill(HIST("evsel"), 2.f);
      }
      if (collision.isInelGt1()) {
        histos.fill(HIST("evsel"), 3.f);
      }
    }
    if (evselOptions.askForCustomTVX) {
      if (!collision.selection_bit(aod::evsel::kIsTriggerTVX)) {
        return false;
      }
    } else {
      if (!collision.sel8()) {
        return false;
      }
    }
    if (evselOptions.requirekIsVertexTOFmatched && !collision.selection_bit(aod::evsel::kIsVertexTOFmatched)) {
      return false;
    }
    if constexpr (fillHistograms) {
      histos.fill(HIST("evsel"), 4.f);
    }
    if (evselOptions.removeITSROFrameBorder && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    if constexpr (fillHistograms) {
      histos.fill(HIST("evsel"), 5.f);
    }
    if (evselOptions.removeNoSameBunchPileup && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if constexpr (fillHistograms) {
      histos.fill(HIST("evsel"), 6.f);
    }
    if (evselOptions.requireIsGoodZvtxFT0vsPV && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if constexpr (fillHistograms) {
      histos.fill(HIST("evsel"), 7.f);
    }
    if (evselOptions.requireIsVertexITSTPC && !collision.selection_bit(aod::evsel::kIsVertexITSTPC)) {
      return false;
    }
    if constexpr (fillHistograms) {
      histos.fill(HIST("evsel"), 8.f);
    }
    if (evselOptions.removeNoTimeFrameBorder && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    if constexpr (fillHistograms) {
      histos.fill(HIST("evsel"), 9.f);
    }
    if constexpr (fillHistograms) {
      histos.fill(HIST("evsel"), 10.f);
      if (collision.isInelGt0()) {
        histos.fill(HIST("evsel"), 11.f);
      }
      if (collision.isInelGt1()) {
        histos.fill(HIST("evsel"), 12.f);
      }
    }
    if (std::abs(collision.posZ()) > evselOptions.cfgCutVertex) {
      return false;
    }
    if constexpr (fillHistograms) {
      histos.fill(HIST("evsel"), 13.f);
      if (collision.isInelGt0()) {
        histos.fill(HIST("evsel"), 14.f);
      } else if (evselOptions.cfgINELCut == 1) {
        return false;
      }
      if (collision.isInelGt1()) {
        histos.fill(HIST("evsel"), 15.f);
      } else if (evselOptions.cfgINELCut == 2) {
        return false;
      }
      histos.fill(HIST("event/vertexz"), collision.posZ());

      if constexpr (fillMultiplicity) {
        histos.fill(HIST("Centrality/FV0A"), collision.centFV0A());
        histos.fill(HIST("Centrality/FT0M"), collision.centFT0M());
        histos.fill(HIST("Centrality/FT0A"), collision.centFT0A());
        histos.fill(HIST("Centrality/FT0C"), collision.centFT0C());
        // histos.fill(HIST("Centrality/FDDM"), collision.centFDDM());
        // histos.fill(HIST("Centrality/NTPV"), collision.centNTPV());

        histos.fill(HIST("Mult/FV0M"), collision.multZeqFV0A());
        histos.fill(HIST("Mult/FT0M"), collision.multZeqFT0A() + collision.multZeqFT0C());
        histos.fill(HIST("Mult/FDDM"), collision.multZeqFDDA() + collision.multZeqFDDC());

        // histos.fill(HIST("Mult/Tracklets"), collision.multTracklets());
        histos.fill(HIST("Mult/TPC"), collision.multTPC());
        histos.fill(HIST("Mult/NTracksPV"), collision.multZeqNTracksPV());
        histos.fill(HIST("Mult/NTracksPVeta1"), collision.multNTracksPVeta1());
      }
    }
    return true;
  }

  template <typename TrackType>
  bool passesDCAxyCut(TrackType const& track) const
  {
    if (useCustomTrackCuts.value) {
      for (int i = 0; i < static_cast<int>(TrackSelection::TrackCuts::kNCuts); i++) {
        if (i == static_cast<int>(TrackSelection::TrackCuts::kDCAxy)) {
          continue;
        }
        if (!customTrackCuts.IsSelected(track, static_cast<TrackSelection::TrackCuts>(i))) {
          return false;
        }
      }
      return (std::abs(track.dcaXY()) <= (maxDcaXYFactor.value * (0.0105f + 0.0350f / std::pow(track.pt(), 1.1f))));
    }
    return track.isGlobalTrack();
  }

  template <typename TrackType>
  bool passesCutWoDCA(TrackType const& track) const
  {
    if (useCustomTrackCuts.value) {
      for (int i = 0; i < static_cast<int>(TrackSelection::TrackCuts::kNCuts); i++) {
        if (i == static_cast<int>(TrackSelection::TrackCuts::kDCAxy)) {
          continue;
        }
        if (i == static_cast<int>(TrackSelection::TrackCuts::kDCAz)) {
          continue;
        }
        if (!customTrackCuts.IsSelected(track, static_cast<TrackSelection::TrackCuts>(i))) {
          return false;
        }
      }
      return true;
    }
    return track.isGlobalTrackWoDCA();
  }

  template <bool fillHistograms = false, typename TrackType, typename CollisionType>
  bool isTrackSelected(TrackType const& track, CollisionType const& /*collision*/)
  {
    if constexpr (fillHistograms) {
      histos.fill(HIST("tracksel"), 1);
    }
    if (track.eta() < trkselOptions.cfgCutEtaMin || track.eta() > trkselOptions.cfgCutEtaMax) {
      return false;
    }
    if constexpr (fillHistograms) {
      histos.fill(HIST("tracksel"), 2);
      if (enableTrackCutHistograms) {
        if (track.sign() > 0) {
          histos.fill(HIST("track/pos/Eta"), track.eta());
        } else {
          histos.fill(HIST("track/neg/Eta"), track.eta());
        }
        if (track.hasITS() && track.hasTPC()) {
          histos.fill(HIST("track/ITS/itsNCls"), track.itsNCls(), track.sign());
          histos.fill(HIST("track/ITS/itsChi2NCl"), track.itsChi2NCl(), track.sign());

          histos.fill(HIST("track/TPC/tpcNClsFindable"), track.tpcNClsFindable(), track.sign());
          histos.fill(HIST("track/TPC/tpcNClsFound"), track.tpcNClsFound(), track.sign());
          histos.fill(HIST("track/TPC/tpcNClsShared"), track.tpcNClsShared(), track.sign());
          histos.fill(HIST("track/TPC/tpcCrossedRows"), track.tpcNClsCrossedRows(), track.sign());
          histos.fill(HIST("track/TPC/tpcCrossedRowsOverFindableCls"), track.tpcCrossedRowsOverFindableCls(), track.sign());
          histos.fill(HIST("track/TPC/tpcFractionSharedCls"), track.tpcFractionSharedCls(), track.sign());
          histos.fill(HIST("track/TPC/tpcChi2NCl"), track.tpcChi2NCl(), track.sign());

          histos.fill(HIST("track/TRD/trdSignal"), track.p(), track.trdSignal(), track.sign());
          if (track.hasTRD()) {
            histos.fill(HIST("track/TRD/length"), track.length());
          } else {
            histos.fill(HIST("track/TRD/lengthnotrd"), track.length());
          }
        }
      }
      if (track.hasITS() && track.isQualityTrackITS() && track.isInAcceptanceTrack()) {
        if (track.sign() > 0) {
          histos.fill(HIST("Data/pos/pt/its"), track.pt());
        } else {
          histos.fill(HIST("Data/neg/pt/its"), track.pt());
        }
      }
      if (track.hasTPC() && track.isQualityTrackTPC() && track.isInAcceptanceTrack()) {
        if (track.sign() > 0) {
          histos.fill(HIST("Data/pos/pt/tpc"), track.pt());
        } else {
          histos.fill(HIST("Data/neg/pt/tpc"), track.pt());
        }
      }
    }

    if (track.tpcChi2NCl() < minChi2PerClusterTPC || track.tpcChi2NCl() > maxChi2PerClusterTPC) {
      return false;
    }

    if (!passesCutWoDCA(track)) {
      return false;
    }

    if constexpr (fillHistograms) {
      histos.fill(HIST("tracksel"), 3);
      if (track.hasTOF()) {
        histos.fill(HIST("tracksel"), 4);
      }
      if (enableTrackCutHistograms) {
        histos.fill(HIST("track/selected/ITS/itsNCls"), track.itsNCls(), track.sign());
        histos.fill(HIST("track/selected/ITS/itsChi2NCl"), track.itsChi2NCl(), track.sign());

        histos.fill(HIST("track/selected/TPC/tpcNClsFindable"), track.tpcNClsFindable(), track.sign());
        histos.fill(HIST("track/selected/TPC/tpcNClsFound"), track.tpcNClsFound(), track.sign());
        histos.fill(HIST("track/selected/TPC/tpcNClsShared"), track.tpcNClsShared(), track.sign());
        histos.fill(HIST("track/selected/TPC/tpcCrossedRows"), track.tpcNClsCrossedRows(), track.sign());
        histos.fill(HIST("track/selected/TPC/tpcCrossedRowsOverFindableCls"), track.tpcCrossedRowsOverFindableCls(), track.sign());
        histos.fill(HIST("track/selected/TPC/tpcFractionSharedCls"), track.tpcFractionSharedCls(), track.sign());
        histos.fill(HIST("track/selected/TPC/tpcChi2NCl"), track.tpcChi2NCl(), track.sign());
      }
    }
    if constexpr (fillHistograms) {
      if (track.hasITS() && track.hasTPC() && track.hasTRD() && track.hasTOF()) {
        if (track.sign() > 0) {
          histos.fill(HIST("Data/pos/pt/its_tpc_trd_tof"), track.pt());
        } else {
          histos.fill(HIST("Data/neg/pt/its_tpc_trd_tof"), track.pt());
        }
      }
      if (track.hasITS() && track.hasTPC() && track.hasTRD()) {
        if (track.sign() > 0) {
          histos.fill(HIST("Data/pos/pt/its_tpc_trd"), track.pt());
        } else {
          histos.fill(HIST("Data/neg/pt/its_tpc_trd"), track.pt());
        }
      }
      if (track.hasITS() && track.hasTPC() && track.hasTOF()) {
        if (track.sign() > 0) {
          histos.fill(HIST("Data/pos/pt/its_tpc_tof"), track.pt());
        } else {
          histos.fill(HIST("Data/neg/pt/its_tpc_tof"), track.pt());
        }
      }
      if (track.hasITS() && track.hasTRD() && track.hasTOF()) {
        if (track.sign() > 0) {
          histos.fill(HIST("Data/pos/pt/its_trd_tof"), track.pt());
        } else {
          histos.fill(HIST("Data/neg/pt/its_trd_tof"), track.pt());
        }
      }
      if (track.hasTPC() && track.hasTRD() && track.hasTOF()) {
        if (track.sign() > 0) {
          histos.fill(HIST("Data/pos/pt/tpc_trd_tof"), track.pt());
        } else {
          histos.fill(HIST("Data/neg/pt/tpc_trd_tof"), track.pt());
        }
      }
      if (track.hasITS() && track.hasTPC()) {
        if (track.sign() > 0) {
          histos.fill(HIST("Data/pos/pt/its_tpc"), track.pt());
        } else {
          histos.fill(HIST("Data/neg/pt/its_tpc"), track.pt());
        }
      }
      if (track.hasTRD() && track.hasTOF()) {
        if (track.sign() > 0) {
          histos.fill(HIST("Data/pos/pt/trd_tof"), track.pt());
        } else {
          histos.fill(HIST("Data/neg/pt/trd_tof"), track.pt());
        }
      }
      if (track.hasTPC() && track.hasTOF()) {
        if (track.sign() > 0) {
          histos.fill(HIST("Data/pos/pt/tpc_tof"), track.pt());
        } else {
          histos.fill(HIST("Data/neg/pt/tpc_tof"), track.pt());
        }
      }
      if (track.hasITS() && track.hasTOF()) {
        if (track.sign() > 0) {
          histos.fill(HIST("Data/pos/pt/its_tof"), track.pt());
        } else {
          histos.fill(HIST("Data/neg/pt/its_tof"), track.pt());
        }
      }
      if (track.hasTPC() && track.hasTRD()) {
        if (track.sign() > 0) {
          histos.fill(HIST("Data/pos/pt/tpc_trd"), track.pt());
        } else {
          histos.fill(HIST("Data/neg/pt/tpc_trd"), track.pt());
        }
      }
      if (track.hasITS() && track.hasTRD()) {
        if (track.sign() > 0) {
          histos.fill(HIST("Data/pos/pt/its_trd"), track.pt());
        } else {
          histos.fill(HIST("Data/neg/pt/its_trd"), track.pt());
        }
      }
    }
    return true;
  }

  template <typename ParticleType>
  bool isMismatchedTrack(const ParticleType& track, const int detector)
  {
    switch (detector) {
      case 0: // ITS
        for (int i = 0; i < 7; i++) {
          if (track.mcMask() & 1 << i) {
            return true;
          }
        }
        return false;
        break;
      case 1: // TPC
        for (int i = 7; i < 10; i++) {
          if (track.mcMask() & 1 << i) {
            return true;
          }
        }
        break;
      case 2: // TRD
        if (track.mcMask() & 1 << 10) {
          return true;
        }
        return false;
        break;
      case 3: // TOF
        if (track.mcMask() & 1 << 11) {
          return true;
        }
        return false;
        break;
      default: // All
        if (track.mcMask() & 1 << 15) {
          return true;
        }
        return false;
    }
    return false;
  }

  using CollisionCandidates = soa::Join<aod::Collisions, aod::EvSels,
                                        aod::TPCMults, aod::PVMults, aod::MultZeqs,
                                        aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>;
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,
                                    aod::pidEvTimeFlags, aod::TrackSelection, aod::TOFSignal>;
  void processMCclosure(CollisionCandidates::iterator const& collisions,
                        soa::Join<TrackCandidates,
                                  aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                                  aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr> const& tracks,
                        aod::McTrackLabels const& mcTrackLabels, aod::McParticles const& mcParticles)
  {
    const float multiplicity = getMultiplicity(collisions);
    // int trackwoCut = 0; int trackwCut = 0;

    for (const auto& track : tracks) {
      if (!track.has_collision()) {
        continue;
      }
      const auto& collision = track.collision_as<CollisionCandidates>();
      if (!isEventSelected<false, false>(collision)) {
        continue;
      }
      if (!isTrackSelected<true>(track, collision)) {
        continue;
      }
      // trackwoCut++;
      if (std::abs(track.dcaXY()) > 0.05) { // Skipping tracks that don't pass the standard cuts
        return;
      }

      // trackwCut++;
      const auto& mcLabel = mcTrackLabels.iteratorAt(track.globalIndex());
      const auto& mcParticle = mcParticles.iteratorAt(mcLabel.mcParticleId());
      int pdgCode = mcParticle.pdgCode();
      const auto& nsigmaTPCPi = o2::aod::pidutils::tpcNSigma<2>(track);
      const auto& nsigmaTPCKa = o2::aod::pidutils::tpcNSigma<3>(track);
      const auto& nsigmaTPCPr = o2::aod::pidutils::tpcNSigma<4>(track);

      bool isTPCPion = std::abs(nsigmaTPCPi) < trkselOptions.cfgCutNsigma;
      bool isTPCKaon = std::abs(nsigmaTPCKa) < trkselOptions.cfgCutNsigma;
      bool isTPCProton = std::abs(nsigmaTPCPr) < trkselOptions.cfgCutNsigma;

      const auto& nsigmaTOFPi = o2::aod::pidutils::tofNSigma<2>(track);
      const auto& nsigmaTOFKa = o2::aod::pidutils::tofNSigma<3>(track);
      const auto& nsigmaTOFPr = o2::aod::pidutils::tofNSigma<4>(track);

      bool isTOFPion = track.hasTOF() && std::abs(nsigmaTOFPi) < trkselOptions.cfgCutNsigma;
      bool isTOFKaon = track.hasTOF() && std::abs(nsigmaTOFKa) < trkselOptions.cfgCutNsigma;
      bool isTOFProton = track.hasTOF() && std::abs(nsigmaTOFPr) < trkselOptions.cfgCutNsigma;
      // Precompute rapidity values to avoid redundant calculations
      double rapidityPi = std::abs(track.rapidity(PID::getMass(2)));
      double rapidityKa = std::abs(track.rapidity(PID::getMass(3)));
      double rapidityPr = std::abs(track.rapidity(PID::getMass(4)));
      if (track.eta() < trkselOptions.cfgCutEtaMin || track.eta() > trkselOptions.cfgCutEtaMax) {
        return;
      }
      if (mcParticle.isPhysicalPrimary()) {
        if (isTPCPion && rapidityPi <= trkselOptions.cfgCutY) {
          if (usePDGcode) {
            if (pdgCode == 211) {
              histos.fill(HIST("nsigmatpc/mc_closure/pos/pi"), track.pt(), nsigmaTPCPi, multiplicity);
            } else if (pdgCode == -211) {
              histos.fill(HIST("nsigmatpc/mc_closure/neg/pi"), track.pt(), nsigmaTPCPi, multiplicity);
            }
          } else {
            histos.fill(HIST("nsigmatpc/mc_closure/pos/pi"), track.pt(), nsigmaTPCPi, multiplicity);
            histos.fill(HIST("nsigmatpc/mc_closure/neg/pi"), track.pt(), nsigmaTPCPi, multiplicity);
          }
        }
        if (isTPCKaon && rapidityKa <= trkselOptions.cfgCutY) {
          if (usePDGcode) {
            if (pdgCode == 321) {
              histos.fill(HIST("nsigmatpc/mc_closure/pos/ka"), track.pt(), nsigmaTPCKa, multiplicity);
            } else if (pdgCode == -321) {
              histos.fill(HIST("nsigmatpc/mc_closure/neg/ka"), track.pt(), nsigmaTPCKa, multiplicity);
            }
          } else {
            histos.fill(HIST("nsigmatpc/mc_closure/pos/ka"), track.pt(), nsigmaTPCKa, multiplicity);
            histos.fill(HIST("nsigmatpc/mc_closure/neg/ka"), track.pt(), nsigmaTPCKa, multiplicity);
          }
        }
        if (isTPCProton && rapidityPr <= trkselOptions.cfgCutY) {
          if (usePDGcode) {
            if (pdgCode == 2212) {
              histos.fill(HIST("nsigmatpc/mc_closure/pos/pr"), track.pt(), nsigmaTPCPr, multiplicity);
            } else if (pdgCode == -2212) {
              histos.fill(HIST("nsigmatpc/mc_closure/neg/pr"), track.pt(), nsigmaTPCPr, multiplicity);
            }
          } else {
            histos.fill(HIST("nsigmatpc/mc_closure/pos/pr"), track.pt(), nsigmaTPCPr, multiplicity);
            histos.fill(HIST("nsigmatpc/mc_closure/neg/pr"), track.pt(), nsigmaTPCPr, multiplicity);
          }
        }

        // TOF Selection and Histogram Filling
        if (isTOFPion && rapidityPi <= trkselOptions.cfgCutY) {
          if (usePDGcode) {
            if (pdgCode == 211) {
              histos.fill(HIST("nsigmatof/mc_closure/pos/pi"), track.pt(), nsigmaTOFPi, multiplicity);
            } else if (pdgCode == -211) {
              histos.fill(HIST("nsigmatof/mc_closure/neg/pi"), track.pt(), nsigmaTOFPi, multiplicity);
            }
          } else {
            histos.fill(HIST("nsigmatof/mc_closure/pos/pi"), track.pt(), nsigmaTOFPi, multiplicity);
            histos.fill(HIST("nsigmatof/mc_closure/neg/pi"), track.pt(), nsigmaTOFPi, multiplicity);
          }
        }
        if (isTOFKaon && rapidityKa <= trkselOptions.cfgCutY) {
          if (usePDGcode) {
            if (pdgCode == 321) {
              histos.fill(HIST("nsigmatof/mc_closure/pos/ka"), track.pt(), nsigmaTOFKa, multiplicity);
            } else if (pdgCode == -321) {
              histos.fill(HIST("nsigmatof/mc_closure/neg/ka"), track.pt(), nsigmaTOFKa, multiplicity);
            }
          } else {
            histos.fill(HIST("nsigmatof/mc_closure/pos/ka"), track.pt(), nsigmaTOFKa, multiplicity);
            histos.fill(HIST("nsigmatof/mc_closure/neg/ka"), track.pt(), nsigmaTOFKa, multiplicity);
          }
        }
        if (isTOFProton && rapidityPr <= trkselOptions.cfgCutY) {
          if (usePDGcode) {
            if (pdgCode == 2212) {
              histos.fill(HIST("nsigmatof/mc_closure/pos/pr"), track.pt(), nsigmaTOFPr, multiplicity);
            } else if (pdgCode == -2212) {
              histos.fill(HIST("nsigmatof/mc_closure/neg/pr"), track.pt(), nsigmaTOFPr, multiplicity);
            }
          } else {
            histos.fill(HIST("nsigmatof/mc_closure/pos/pr"), track.pt(), nsigmaTOFPr, multiplicity);
            histos.fill(HIST("nsigmatof/mc_closure/neg/pr"), track.pt(), nsigmaTOFPr, multiplicity);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(tofSpectra, processMCclosure, "MC closure test", false);

  void processOccupancy(CollisionCandidates::iterator const& collision,
                        soa::Join<TrackCandidates,
                                  aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                                  aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr> const& tracks)
  {
    // Event selection criteria
    /*if (!collision.sel8() || std::abs(collision.posZ()) > 10 || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder) ||
        !collision.selection_bit(aod::evsel::kNoTimeFrameBorder) ||
        !collision.selection_bit(aod::evsel::kNoSameBunchPileup) ||
        !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
        return;
    }*/
    if (!isEventSelected<true, true>(collision)) {
      return;
    }
    histos.fill(HIST("test_occupancy/event/vertexz"), collision.posZ());

    // Multiplicity and occupancy
    int occupancy = collision.trackOccupancyInTimeRange();
    const float multiplicity = getMultiplicity(collision);
    histos.fill(HIST("nsigmatpc/test_occupancy/Mult_vs_Occupancy"), multiplicity, occupancy);

    int tpcCount = 0, tofCount = 0;

    for (const auto& track : tracks) {
      // Track selection criteria
      /*   if (track.tpcNClsCrossedRows() < minNCrossedRowsTPC || track.tpcChi2NCl() > maxChi2PerClusterTPC || track.tpcChi2NCl() > maxChi2PerClusterTPC ||
     track.itsChi2NCl() > maxChi2PerClusterITS || std::abs(track.dcaXY()) > maxDcaXYFactor.value * (0.0105f + 0.0350f / pow(track.pt(), 1.1f)) || std::abs(track.dcaZ()) > maxDcaZ.value || track.eta() < trkselOptions.cfgCutEtaMin || track.eta() > trkselOptions.cfgCutEtaMax || track.tpcCrossedRowsOverFindableCls() < minNCrossedRowsOverFindableClustersTPC || track.tpcNClsFound() < minTPCNClsFound ||
     !(o2::aod::track::ITSrefit) || !(o2::aod::track::TPCrefit)) {
     continue;
 }*/
      if (!isTrackSelected<true>(track, collision)) {
        continue;
      }
      const auto& nsigmaTPCPi = o2::aod::pidutils::tpcNSigma<2>(track);
      const auto& nsigmaTPCKa = o2::aod::pidutils::tpcNSigma<3>(track);
      const auto& nsigmaTPCPr = o2::aod::pidutils::tpcNSigma<4>(track);

      const bool isTPCPion = std::abs(nsigmaTPCPi) < trkselOptions.cfgCutNsigma;
      const bool isTPCKaon = std::abs(nsigmaTPCKa) < trkselOptions.cfgCutNsigma;
      const bool isTPCProton = std::abs(nsigmaTPCPr) < trkselOptions.cfgCutNsigma;

      const auto& nsigmaTOFPi = o2::aod::pidutils::tofNSigma<2>(track);
      const auto& nsigmaTOFKa = o2::aod::pidutils::tofNSigma<3>(track);
      const auto& nsigmaTOFPr = o2::aod::pidutils::tofNSigma<4>(track);

      const bool isTOFPion = track.hasTOF() && std::abs(nsigmaTOFPi) < trkselOptions.cfgCutNsigma;
      const bool isTOFKaon = track.hasTOF() && std::abs(nsigmaTOFKa) < trkselOptions.cfgCutNsigma;
      const bool isTOFProton = track.hasTOF() && std::abs(nsigmaTOFPr) < trkselOptions.cfgCutNsigma;

      // Apply rapidity cut for identified particles
      if (isTPCPion && std::abs(track.rapidity(PID::getMass(2))) < trkselOptions.cfgCutY) {
        tpcCount++;
        if (track.sign() > 0) {
          histos.fill(HIST("nsigmatpc/test_occupancy/pos/pi"), track.pt(), nsigmaTPCPi, multiplicity, occupancy);
        } else {
          histos.fill(HIST("nsigmatpc/test_occupancy/neg/pi"), track.pt(), nsigmaTPCPi, multiplicity, occupancy);
        }
      } else if (isTPCKaon && std::abs(track.rapidity(PID::getMass(3))) < trkselOptions.cfgCutY) {
        tpcCount++;
        if (track.sign() > 0) {
          histos.fill(HIST("nsigmatpc/test_occupancy/pos/ka"), track.pt(), nsigmaTPCKa, multiplicity, occupancy);
        } else {
          histos.fill(HIST("nsigmatpc/test_occupancy/neg/ka"), track.pt(), nsigmaTPCKa, multiplicity, occupancy);
        }
      } else if (isTPCProton && std::abs(track.rapidity(PID::getMass(4))) < trkselOptions.cfgCutY) {
        tpcCount++;
        if (track.sign() > 0) {
          histos.fill(HIST("nsigmatpc/test_occupancy/pos/pr"), track.pt(), nsigmaTPCPr, multiplicity, occupancy);
        } else {
          histos.fill(HIST("nsigmatpc/test_occupancy/neg/pr"), track.pt(), nsigmaTPCPr, multiplicity, occupancy);
        }
      }

      // TOF PID histograms
      if (isTOFPion && std::abs(track.rapidity(PID::getMass(2))) < trkselOptions.cfgCutY) {
        tofCount++;
        if (track.sign() > 0) {
          histos.fill(HIST("nsigmatof/test_occupancy/pos/pi"), track.pt(), nsigmaTOFPi, multiplicity, occupancy);
        } else {
          histos.fill(HIST("nsigmatof/test_occupancy/neg/pi"), track.pt(), nsigmaTOFPi, multiplicity, occupancy);
        }
      } else if (isTOFKaon && std::abs(track.rapidity(PID::getMass(3))) < trkselOptions.cfgCutY) {
        tofCount++;
        if (track.sign() > 0) {
          histos.fill(HIST("nsigmatof/test_occupancy/pos/ka"), track.pt(), nsigmaTOFKa, multiplicity, occupancy);
        } else {
          histos.fill(HIST("nsigmatof/test_occupancy/neg/ka"), track.pt(), nsigmaTOFKa, multiplicity, occupancy);
        }
      } else if (isTOFProton && std::abs(track.rapidity(PID::getMass(4))) < trkselOptions.cfgCutY) {
        tofCount++;
        if (track.sign() > 0) {
          histos.fill(HIST("nsigmatof/test_occupancy/pos/pr"), track.pt(), nsigmaTOFPr, multiplicity, occupancy);
        } else {
          histos.fill(HIST("nsigmatof/test_occupancy/neg/pr"), track.pt(), nsigmaTOFPr, multiplicity, occupancy);
        }
      }
    }
    histos.fill(HIST("test_occupancy/tpcCount"), tpcCount);
    histos.fill(HIST("test_occupancy/tofCount"), tofCount);
  } // process function
  PROCESS_SWITCH(tofSpectra, processOccupancy, "check for occupancy plots", true);

  void processStandard(CollisionCandidates::iterator const& collision,
                       TrackCandidates const& tracks)
  {
    if (!isEventSelected<true, true>(collision)) {
      return;
    }
    hMultiplicityvsPercentile->Fill(getMultiplicity(collision), collision.multNTracksPV());
    for (const auto& track : tracks) {
      if (!isTrackSelected<true>(track, collision)) {
        continue;
      }
    }
  } // end of the process function
  PROCESS_SWITCH(tofSpectra, processStandard, "Standard processor from AO2D", true);

  Preslice<aod::SpTracks> spPerCol = aod::spectra::collisionId;
  SliceCache cacheTrk;
  void processDerived(aod::SpColls const& collisions,
                      aod::SpTracks const& tracks)
  {
    for (const auto& collision : collisions) {
      if (!isEventSelected<true, true>(collision)) {
        return;
      }
      const auto& tracksInCollision = tracks.sliceByCached(aod::spectra::collisionId, collision.globalIndex(), cacheTrk);
      for (const auto& track : tracksInCollision) {
        if (!isTrackSelected<true>(track, collision)) {
          continue;
        }
        fillParticleHistos<false, PID::Pion>(track, collision);
        fillParticleHistos<false, PID::Kaon>(track, collision);
        fillParticleHistos<false, PID::Proton>(track, collision);
      }
    }
  } // end of the process function
  PROCESS_SWITCH(tofSpectra, processDerived, "Derived data processor", false);

#define makeProcessFunction(processorName, inputPid, particleId, isFull, tofTable, tpcTable)   \
  void process##processorName##inputPid(CollisionCandidates::iterator const& collision,        \
                                        soa::Join<TrackCandidates,                             \
                                                  aod::pid##tofTable##inputPid,                \
                                                  aod::pid##tpcTable##inputPid> const& tracks) \
  {                                                                                            \
    if (!isEventSelected<false, false>(collision)) {                                           \
      return;                                                                                  \
    }                                                                                          \
    for (const auto& track : tracks) {                                                         \
      if (!isTrackSelected<false>(track, collision)) {                                         \
        continue;                                                                              \
      }                                                                                        \
      fillParticleHistos<isFull, PID::particleId>(track, collision);                           \
    }                                                                                          \
  }                                                                                            \
  PROCESS_SWITCH(tofSpectra, process##processorName##inputPid, Form("Process for the %s hypothesis from %s tables", #particleId, #processorName), false);

// Full tables
#define makeProcessFunctionFull(inputPid, particleId) makeProcessFunction(Full, inputPid, particleId, true, TOFFull, TPCFull)

  makeProcessFunctionFull(El, Electron);
  makeProcessFunctionFull(Mu, Muon);
  makeProcessFunctionFull(Pi, Pion);
  makeProcessFunctionFull(Ka, Kaon);
  makeProcessFunctionFull(Pr, Proton);
  makeProcessFunctionFull(De, Deuteron);
  makeProcessFunctionFull(Tr, Triton);
  makeProcessFunctionFull(He, Helium3);
  makeProcessFunctionFull(Al, Alpha);
#undef makeProcessFunctionFull

// Full LF tables
#define makeProcessFunctionFull(inputPid, particleId) makeProcessFunction(LfFull, inputPid, particleId, true, TOFFull, TPCLfFull)

  makeProcessFunctionFull(El, Electron);
  makeProcessFunctionFull(Mu, Muon);
  makeProcessFunctionFull(Pi, Pion);
  makeProcessFunctionFull(Ka, Kaon);
  makeProcessFunctionFull(Pr, Proton);
  makeProcessFunctionFull(De, Deuteron);
  makeProcessFunctionFull(Tr, Triton);
  makeProcessFunctionFull(He, Helium3);
  makeProcessFunctionFull(Al, Alpha);
#undef makeProcessFunctionFull

  template <typename CollisionType, bool isMC = false>
  float getMultiplicity(const CollisionType& collision)
  {

    switch (multiplicityEstimator) {
      case MultCodes::kNoMultiplicity: // No multiplicity
        return 50.f;                   // to check if its filled
        break;
      case MultCodes::kMultFV0M: // MultFV0M
        if constexpr (!isMC) {
          return collision.multZeqFV0A();
        } else {
          return 50.f; // Not implemented yet
        }
        break;
      case MultCodes::kMultFT0M:
        if constexpr (!isMC) {
          return collision.multZeqFT0A() + collision.multZeqFT0C();
        } else {
          return 50.f; // Not implemented yet
        }
        break;
      case MultCodes::kMultFDDM: // MultFDDM
        if constexpr (!isMC) {
          return collision.multZeqFDDA() + collision.multZeqFDDC();
        } else {
          return 50.f; // Not implemented yet
        }
        break;
      case MultCodes::kMultTracklets: // MultTracklets
        if constexpr (!isMC) {
          // return collision.multTracklets();
        } else {
          return 50.f; // Not implemented yet
        }
        return 0.f; // Undefined in Run3
        break;
      case MultCodes::kMultTPC: // MultTPC
        if constexpr (!isMC) {
          return collision.multTPC();
        } else {
          return 50.f; // Not implemented yet
        }
        break;
      case MultCodes::kMultNTracksPV: // MultNTracksPV
        if constexpr (!isMC) {
          // return collision.multNTracksPV();
          return collision.multZeqNTracksPV();
        } else {
          return 50.f; // Not implemented yet
        }
        break;
      case MultCodes::kMultNTracksPVeta1: // MultNTracksPVeta1
        if constexpr (!isMC) {
          return collision.multNTracksPVeta1();
        } else {
          return 50.f; // Not implemented yet
        }
        break;
      case MultCodes::kCentralityFT0C: // Centrality FT0C
        if constexpr (!isMC) {
          return collision.centFT0C();
        } else {
          return 50.f; // Not implemented yet
        }
        break;
      case MultCodes::kCentralityFT0M: // Centrality FT0M
        return collision.centFT0M();
        break;
      default:
        LOG(fatal) << "Unknown multiplicity estimator: " << multiplicityEstimator;
        return 0.f;
    }
  }

  using GenMCCollisions = soa::Join<aod::McCollisions, aod::McCentFT0Ms, aod::MultsExtraMC>;
  float getMultiplicityMC(const GenMCCollisions::iterator& collision) { return getMultiplicity<GenMCCollisions::iterator, true>(collision); }

  template <std::size_t id>
  bool isParticleEnabled()
  {
    if constexpr (id == 0 || id == Np) {
      if (doprocessFullEl == true || doprocessLfFullEl == true) {
        return true;
      }
    } else if constexpr (id == 1 || id == Np + 1) {
      if (doprocessFullMu == true || doprocessLfFullMu == true) {
        return true;
      }
    } else if constexpr (id == 2 || id == Np + 2) {
      if (doprocessFullPi == true || doprocessLfFullPi == true) {
        return true;
      }
    } else if constexpr (id == 3 || id == Np + 3) {
      if (doprocessFullKa == true || doprocessLfFullKa == true) {
        return true;
      }
    } else if constexpr (id == 4 || id == Np + 4) {
      if (doprocessFullPr == true || doprocessLfFullPr == true) {
        return true;
      }
    } else if constexpr (id == 5 || id == Np + 5) {
      if (doprocessFullDe == true || doprocessLfFullDe == true) {
        return true;
      }
    } else if constexpr (id == 6 || id == Np + 6) {
      if (doprocessFullTr == true || doprocessLfFullTr == true) {
        return true;
      }
    } else if constexpr (id == 7 || id == Np + 7) {
      if (doprocessFullHe == true || doprocessLfFullHe == true) {
        return true;
      }
    } else if constexpr (id == 8 || id == Np + 8) {
      if (doprocessFullAl == true || doprocessLfFullAl == true) {
        return true;
      }
    } else {
      LOG(fatal) << "Unknown particle id: " << id;
    }
    return false;
  }

  using RecoMCCollisions = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0As, aod::CentFT0Cs, aod::TPCMults, aod::PVMults, aod::MultZeqs, aod::CentFT0Ms>; // RD
  template <std::size_t i, typename TrackType, typename ParticleType>
  void fillTrackHistograms_MC(TrackType const& track,
                              ParticleType::iterator const& mcParticle,
                              RecoMCCollisions::iterator const& collision,
                              ParticleType const& mcParticles)
  {
    if (!isParticleEnabled<i>()) { // Check if the particle is enabled
      return;
    }
    if (!collision.has_mcCollision()) {
      return; // Skips processing if no corresponding MC collision is found (rare case!)
    }
    const auto& mcCollision = collision.mcCollision_as<GenMCCollisions>();
    const float multiplicity = getMultiplicity(collision);
    const int occupancy = collision.trackOccupancyInTimeRange();
    //************************************RD**************************************************
    const float impParam = mcCollision.impactParameter();
    //************************************RD**************************************************

    if (mcParticle.pdgCode() != PDGs[i]) {
      return;
    }

    if (track.eta() < trkselOptions.cfgCutEtaMin || track.eta() > trkselOptions.cfgCutEtaMax) {
      return;
    }

    if (std::abs(mcParticle.y()) > trkselOptions.cfgCutY) {
      return;
    }
    if (enablePureDCAHistogram) {
      const auto& nsigmaTPCKa = o2::aod::pidutils::tpcNSigma<3>(track);
      const auto& nsigmaTOFKa = o2::aod::pidutils::tofNSigma<3>(track);

      // Filling DCA info with the TPC+TOF PID
      bool isDCAPureSample = (std::sqrt(nsigmaTOFKa * nsigmaTOFKa + nsigmaTPCKa * nsigmaTPCKa) < 2.f);
      if (track.pt() <= 0.4) {
        isDCAPureSample = (nsigmaTPCKa < 1.f);
      }

      if (isDCAPureSample) {
        if (enableDCAvsmotherHistograms) {
          hDcaXYMC[i]->Fill(track.pt(), track.dcaXY());
          hDcaZMC[i]->Fill(track.pt(), track.dcaZ());
        }

        if (!mcParticle.isPhysicalPrimary()) { // Secondaries (weak decays and material)
          if (mcParticle.getProcess() == 4) {  // Particles from decay
            if (enableDCAxyzHistograms) {
              hDcaXYZStr[i]->Fill(track.pt(), track.dcaXY(), track.dcaZ());
            } else {
              histos.fill(HIST(hdcaxystr[i]), track.pt(), track.dcaXY());
              histos.fill(HIST(hdcazstr[i]), track.pt(), track.dcaZ());
            }

            if (mcParticle.has_mothers()) {
              for (const auto& mother : mcParticle.template mothers_as<aod::McParticles>()) {
                auto daughter0 = mother.template daughters_as<aod::McParticles>().begin();
                double vertexDau[3] = {daughter0.vx(), daughter0.vy(), daughter0.vz()};
                double vertexMoth[3] = {mother.vx(), mother.vy(), mother.vz()};
                auto decayLength = RecoDecay::distance(vertexMoth, vertexDau);
                hDecayLengthStr[i]->Fill(track.pt(), decayLength);
              }
            }
          } else { // Particles from the material
            if (enableDCAxyzHistograms) {
              hDcaXYZMat[i]->Fill(track.pt(), track.dcaXY(), track.dcaZ());
            } else {
              histos.fill(HIST(hdcaxymat[i]), track.pt(), track.dcaXY());
              histos.fill(HIST(hdcazmat[i]), track.pt(), track.dcaZ());
            }
          }
        } else { // Primaries
          if (enableDCAxyzHistograms) {
            hDcaXYZPrm[i]->Fill(track.pt(), track.dcaXY(), track.dcaZ());
            if (enableDcaGoodEvents.value && collision.has_mcCollision()) {
              histos.fill(HIST(hdcaxyprmgoodevs[i]), track.pt(), track.dcaXY(), track.dcaZ());
            }
          } else {
            // DCAxy for all primaries
            histos.fill(HIST(hdcaxyprm[i]), track.pt(), track.dcaXY());
            histos.fill(HIST(hdcazprm[i]), track.pt(), track.dcaZ());
          }
          if (enableDcaGoodEvents.value && collision.has_mcCollision()) {
            histos.fill(HIST(hdcaxyprmgoodevs[i]), track.pt(), track.dcaXY());
            histos.fill(HIST(hdcazprmgoodevs[i]), track.pt(), track.dcaZ());
          }

          if (enableDCAvsmotherHistograms) {
            bool IsD0Mother = false;
            bool IsCharmMother = false;
            bool IsBeautyMother = false;
            bool IsNotHFMother = false;
            if (mcParticle.has_mothers()) {
              const int charmOrigin = RecoDecay::getCharmHadronOrigin(mcParticles, mcParticle, false);
              for (const auto& mother : mcParticle.template mothers_as<aod::McParticles>()) {
                const int motherPdgCode = std::abs(mother.pdgCode());
                if (motherPdgCode == 421) {
                  IsD0Mother = true;
                }
                if (charmOrigin == RecoDecay::OriginType::NonPrompt) {
                  IsBeautyMother = true;
                }
                if (charmOrigin == RecoDecay::OriginType::Prompt) {
                  IsCharmMother = true;
                }
                if (charmOrigin == RecoDecay::OriginType::None) {
                  IsNotHFMother = true;
                }
              }
            }
            if (IsD0Mother) {
              hDcaXYMCD0[i]->Fill(track.pt(), track.dcaXY());
              hDcaZMCD0[i]->Fill(track.pt(), track.dcaZ());
            }
            if (IsCharmMother) {
              hDcaXYMCCharm[i]->Fill(track.pt(), track.dcaXY());
              hdcaZMCCharm[i]->Fill(track.pt(), track.dcaZ());
            }
            if (IsBeautyMother) {
              hDcaXYMCBeauty[i]->Fill(track.pt(), track.dcaXY());
              hDcaZMCBeauty[i]->Fill(track.pt(), track.dcaZ());
            }
            if (IsNotHFMother) {
              hDcaXYMCNotHF[i]->Fill(track.pt(), track.dcaXY());
              hDcaZMCNotHF[i]->Fill(track.pt(), track.dcaZ());
            }

            if (mcParticle.has_mothers()) {
              for (const auto& mother : mcParticle.template mothers_as<aod::McParticles>()) {
                auto daughter0 = mother.template daughters_as<aod::McParticles>().begin();
                double vertexDau[3] = {daughter0.vx(), daughter0.vy(), daughter0.vz()};
                double vertexMoth[3] = {mother.vx(), mother.vy(), mother.vz()};
                auto decayLength = RecoDecay::distance(vertexMoth, vertexDau);

                if (IsD0Mother) {
                  hDecayLengthMCD0[i]->Fill(track.pt(), decayLength);
                }
                if (IsCharmMother) {
                  hDecayLengthMCCharm[i]->Fill(track.pt(), decayLength);
                }
                if (IsBeautyMother) {
                  hDecayLengthMCBeauty[i]->Fill(track.pt(), decayLength);
                }
                if (IsNotHFMother) {
                  hDecayLengthMCNotHF[i]->Fill(track.pt(), decayLength);
                }
              }
            }
          }
        }
      }
    } else {
      if (!mcParticle.isPhysicalPrimary()) {
        if (mcParticle.getProcess() == 4) {
          histos.fill(HIST(hdcaxystr[i]), track.pt(), track.dcaXY());
        } else {
          histos.fill(HIST(hdcaxymat[i]), track.pt(), track.dcaXY());
        }
      } else {
        histos.fill(HIST(hdcaxyprm[i]), track.pt(), track.dcaXY());
      }
    }

    if ((collision.has_mcCollision() && (mcParticle.mcCollisionId() != collision.mcCollisionId())) || !collision.has_mcCollision()) {
      if (!mcParticle.isPhysicalPrimary()) {
        if (mcParticle.getProcess() == 4) {
          hDcaXYWrongCollisionStr[i]->Fill(track.pt(), track.dcaXY());
        } else {
          hDcaXYWrongCollisionMat[i]->Fill(track.pt(), track.dcaXY());
        }
      } else {
        hDcaXYWrongCollisionPrm[i]->Fill(track.pt(), track.dcaXY());
      }
    }

    if (!passesDCAxyCut(track)) { // Skipping tracks that don't pass the standard cuts
      return;
    }
    const int pdgCode = mcParticle.pdgCode();
    const auto& nsigmaTPCPi = o2::aod::pidutils::tpcNSigma<2>(track);
    const auto& nsigmaTPCKa = o2::aod::pidutils::tpcNSigma<3>(track);
    const auto& nsigmaTPCPr = o2::aod::pidutils::tpcNSigma<4>(track);

    const bool isPionTPC = std::abs(nsigmaTPCPi) < trkselOptions.cfgCutNsigma;
    const bool isKaonTPC = std::abs(nsigmaTPCKa) < trkselOptions.cfgCutNsigma;
    const bool isProtonTPC = std::abs(nsigmaTPCPr) < trkselOptions.cfgCutNsigma;

    const auto& nsigmaTOFPi = o2::aod::pidutils::tofNSigma<2>(track);
    const auto& nsigmaTOFKa = o2::aod::pidutils::tofNSigma<3>(track);
    const auto& nsigmaTOFPr = o2::aod::pidutils::tofNSigma<4>(track);

    const bool isPionTOF = std::abs(nsigmaTOFPi) < trkselOptions.cfgCutNsigma;
    const bool isKaonTOF = std::abs(nsigmaTOFKa) < trkselOptions.cfgCutNsigma;
    const bool isProtonTOF = std::abs(nsigmaTOFPr) < trkselOptions.cfgCutNsigma;

    if (!mcParticle.isPhysicalPrimary()) { // Is not physical primary
      if (mcParticle.getProcess() == 4) {  // Is from decay
        if (includeCentralityMC) {
          if (includeCentralityMC) {
            histos.fill(HIST(hpt_num_str[i]), track.pt(), multiplicity, track.dcaXY());
          }
        } else {
          histos.fill(HIST(hpt_num_str[i]), track.pt(), multiplicity);
        }
        if (track.hasTOF()) {
          if (includeCentralityMC) {
            histos.fill(HIST(hpt_numtof_str[i]), track.pt(), multiplicity, track.dcaXY());
          } else {
            histos.fill(HIST(hpt_numtof_str[i]), track.pt(), multiplicity);
          }
        }
      } else {
        if (includeCentralityMC) {
          histos.fill(HIST(hpt_num_mat[i]), track.pt(), multiplicity, track.dcaXY());
          if (track.hasTOF()) {
            histos.fill(HIST(hpt_numtof_mat[i]), track.pt(), multiplicity, track.dcaXY());
          }
        } else {
          histos.fill(HIST(hpt_num_mat[i]), track.pt(), multiplicity);
          if (track.hasTOF()) {
            histos.fill(HIST(hpt_numtof_mat[i]), track.pt(), multiplicity);
          }
        }
      }
    } else { // Is physical primary
      if (includeCentralityMC) {
        if (isImpactParam) {
          histos.fill(HIST(hpt_num_prm[i]), track.pt(), impParam, track.dcaXY(), occupancy);
        } else {
          histos.fill(HIST(hpt_num_prm[i]), track.pt(), multiplicity, track.dcaXY(), occupancy);
        }
      } else {
        histos.fill(HIST(hpt_num_prm[i]), track.pt(), multiplicity);
      }
      if (isPionTPC || isKaonTPC || isProtonTPC) {
        if (pdgCode == 2212) {
          if (isImpactParam) {
            histos.fill(HIST("MC/withPID/pr/pos/prm/pt/num"), track.pt(), impParam);
            if (!mcParticle.isPhysicalPrimary()) {
              if (mcParticle.getProcess() == 4) {
                histos.fill(HIST("MC/withPID/pr/pos/prm/pt/num_str"), track.pt(), impParam);
              } else {
                histos.fill(HIST("MC/withPID/pr/pos/prm/pt/num_mat"), track.pt(), impParam);
              }
            }
          } else {
            histos.fill(HIST("MC/withPID/pr/pos/prm/pt/num"), track.pt(), multiplicity);
            if (!mcParticle.isPhysicalPrimary()) {
              if (mcParticle.getProcess() == 4) {
                histos.fill(HIST("MC/withPID/pr/pos/prm/pt/num_str"), track.pt(), multiplicity);
              } else {
                histos.fill(HIST("MC/withPID/pr/pos/prm/pt/num_mat"), track.pt(), multiplicity);
              }
            }
          }
        } else if (pdgCode == -2212) {
          if (isImpactParam) {
            histos.fill(HIST("MC/withPID/pr/neg/prm/pt/num"), track.pt(), impParam);
            if (!mcParticle.isPhysicalPrimary()) {
              if (mcParticle.getProcess() == 4) {
                histos.fill(HIST("MC/withPID/pr/neg/prm/pt/num_str"), track.pt(), impParam);
              } else {
                histos.fill(HIST("MC/withPID/pr/neg/prm/pt/num_mat"), track.pt(), impParam);
              }
            }
          } else {
            histos.fill(HIST("MC/withPID/pr/neg/prm/pt/num"), track.pt(), multiplicity);
            if (!mcParticle.isPhysicalPrimary()) {
              if (mcParticle.getProcess() == 4) {
                histos.fill(HIST("MC/withPID/pr/neg/prm/pt/num_str"), track.pt(), multiplicity);
              } else {
                histos.fill(HIST("MC/withPID/pr/neg/prm/pt/num_mat"), track.pt(), multiplicity);
              }
            }
          }
        } else if (pdgCode == 211) {
          if (isImpactParam) {
            histos.fill(HIST("MC/withPID/pi/pos/prm/pt/num"), track.pt(), impParam);
            if (!mcParticle.isPhysicalPrimary()) {
              if (mcParticle.getProcess() == 4) {
                histos.fill(HIST("MC/withPID/pi/pos/prm/pt/num_str"), track.pt(), impParam);
              } else {
                histos.fill(HIST("MC/withPID/pi/pos/prm/pt/num_mat"), track.pt(), impParam);
              }
            }
          } else {
            histos.fill(HIST("MC/withPID/pi/pos/prm/pt/num"), track.pt(), multiplicity);
            if (!mcParticle.isPhysicalPrimary()) {
              if (mcParticle.getProcess() == 4) {
                histos.fill(HIST("MC/withPID/pi/pos/prm/pt/num_str"), track.pt(), multiplicity);
              } else {
                histos.fill(HIST("MC/withPID/pi/pos/prm/pt/num_mat"), track.pt(), multiplicity);
              }
            }
          }
        } else if (pdgCode == -211) {
          if (isImpactParam) {
            histos.fill(HIST("MC/withPID/pi/neg/prm/pt/num"), track.pt(), impParam);
            if (!mcParticle.isPhysicalPrimary()) {
              if (mcParticle.getProcess() == 4) {
                histos.fill(HIST("MC/withPID/pi/neg/prm/pt/num_str"), track.pt(), impParam);
              } else {
                histos.fill(HIST("MC/withPID/pi/neg/prm/pt/num_mat"), track.pt(), impParam);
              }
            }
          } else {
            histos.fill(HIST("MC/withPID/pi/neg/prm/pt/num"), track.pt(), multiplicity);
            if (!mcParticle.isPhysicalPrimary()) {
              if (mcParticle.getProcess() == 4) {
                histos.fill(HIST("MC/withPID/pi/neg/prm/pt/num_str"), track.pt(), multiplicity);
              } else {
                histos.fill(HIST("MC/withPID/pi/neg/prm/pt/num_mat"), track.pt(), multiplicity);
              }
            }
          }
        } else if (pdgCode == 321) {
          if (isImpactParam) {
            histos.fill(HIST("MC/withPID/ka/pos/prm/pt/num"), track.pt(), impParam);
            if (!mcParticle.isPhysicalPrimary()) {
              if (mcParticle.getProcess() == 4) {
                histos.fill(HIST("MC/withPID/ka/pos/prm/pt/num_str"), track.pt(), impParam);
              } else {
                histos.fill(HIST("MC/withPID/ka/pos/prm/pt/num_mat"), track.pt(), impParam);
              }
            }
          } else {
            histos.fill(HIST("MC/withPID/ka/pos/prm/pt/num"), track.pt(), multiplicity);
            if (!mcParticle.isPhysicalPrimary()) {
              if (mcParticle.getProcess() == 4) {
                histos.fill(HIST("MC/withPID/ka/pos/prm/pt/num_str"), track.pt(), multiplicity);
              } else {
                histos.fill(HIST("MC/withPID/ka/pos/prm/pt/num_mat"), track.pt(), multiplicity);
              }
            }
          }
        } else if (pdgCode == -321) {
          if (isImpactParam) {
            histos.fill(HIST("MC/withPID/ka/neg/prm/pt/num"), track.pt(), impParam);
            if (!mcParticle.isPhysicalPrimary()) {
              if (mcParticle.getProcess() == 4) {
                histos.fill(HIST("MC/withPID/ka/neg/prm/pt/num_str"), track.pt(), impParam);
              } else {
                histos.fill(HIST("MC/withPID/ka/neg/prm/pt/num_mat"), track.pt(), impParam);
              }
            }
          } else {
            histos.fill(HIST("MC/withPID/ka/neg/prm/pt/num"), track.pt(), multiplicity);
            if (!mcParticle.isPhysicalPrimary()) {
              if (mcParticle.getProcess() == 4) {
                histos.fill(HIST("MC/withPID/ka/neg/prm/pt/num_str"), track.pt(), multiplicity);
              } else {
                histos.fill(HIST("MC/withPID/ka/neg/prm/pt/num_mat"), track.pt(), multiplicity);
              }
            }
          }
        }
      }
      if (track.hasTRD() && trkselOptions.lastRequiredTrdCluster > 0) {
        int lastLayer = 0;
        for (int l = 7; l >= 0; l--) {
          if (track.trdPattern() & (1 << l)) {
            lastLayer = l;
            break;
          }
        }
        if (lastLayer < trkselOptions.lastRequiredTrdCluster) {
          return;
        }
      }
      if (track.hasTOF()) {
        if (isPionTOF || isKaonTOF || isProtonTOF) {
          // Proton (positive)
          if (pdgCode == 2212) {
            if (isImpactParam) {
              histos.fill(HIST("MC/withPID/pr/pos/prm/pt/numtof"), track.pt(), impParam);
            } else {
              histos.fill(HIST("MC/withPID/pr/pos/prm/pt/numtof"), track.pt(), multiplicity);
            }
            // Matched proton condition
            if (!(track.mcMask() & (1 << 11))) {
              if (isImpactParam) {
                histos.fill(HIST("MC/withPID/pr/pos/prm/pt/numtof_matched"), track.pt(), impParam);
              } else {
                histos.fill(HIST("MC/withPID/pr/pos/prm/pt/numtof_matched"), track.pt(), multiplicity);
              }
            }
          } else if (pdgCode == -2212) {
            if (isImpactParam) {
              histos.fill(HIST("MC/withPID/pr/neg/prm/pt/numtof"), track.pt(), impParam);
            } else {
              histos.fill(HIST("MC/withPID/pr/neg/prm/pt/numtof"), track.pt(), multiplicity);
            }
            if (!(track.mcMask() & (1 << 11))) {
              if (isImpactParam) {
                histos.fill(HIST("MC/withPID/pr/neg/prm/pt/numtof_matched"), track.pt(), impParam);
              } else {
                histos.fill(HIST("MC/withPID/pr/neg/prm/pt/numtof_matched"), track.pt(), multiplicity);
              }
            }
          } else if (pdgCode == 211) {
            if (isImpactParam) {
              histos.fill(HIST("MC/withPID/pi/pos/prm/pt/numtof"), track.pt(), impParam);
            } else {
              histos.fill(HIST("MC/withPID/pi/pos/prm/pt/numtof"), track.pt(), multiplicity);
            }
            // Matched pion condition
            if (!(track.mcMask() & (1 << 11))) {
              if (isImpactParam) {
                histos.fill(HIST("MC/withPID/pi/pos/prm/pt/numtof_matched"), track.pt(), impParam);
              } else {
                histos.fill(HIST("MC/withPID/pi/pos/prm/pt/numtof_matched"), track.pt(), multiplicity);
              }
            }
          } else if (pdgCode == -211) {
            if (isImpactParam) {
              histos.fill(HIST("MC/withPID/pi/neg/prm/pt/numtof"), track.pt(), impParam);
            } else {
              histos.fill(HIST("MC/withPID/pi/neg/prm/pt/numtof"), track.pt(), multiplicity);
            }
            // Matched pion condition
            if (!(track.mcMask() & (1 << 11))) {
              if (isImpactParam) {
                histos.fill(HIST("MC/withPID/pi/neg/prm/pt/numtof_matched"), track.pt(), impParam);
              } else {
                histos.fill(HIST("MC/withPID/pi/neg/prm/pt/numtof_matched"), track.pt(), multiplicity);
              }
            }
          } else if (pdgCode == 321) {
            if (isImpactParam) {
              histos.fill(HIST("MC/withPID/ka/pos/prm/pt/numtof"), track.pt(), impParam);
            } else {
              histos.fill(HIST("MC/withPID/ka/pos/prm/pt/numtof"), track.pt(), multiplicity);
            }
            // Matched kaon condition
            if (!(track.mcMask() & (1 << 11))) {
              if (isImpactParam) {
                histos.fill(HIST("MC/withPID/ka/pos/prm/pt/numtof_matched"), track.pt(), impParam);
              } else {
                histos.fill(HIST("MC/withPID/ka/pos/prm/pt/numtof_matched"), track.pt(), multiplicity);
              }
            }
          } else if (pdgCode == -321) {
            if (isImpactParam) {
              histos.fill(HIST("MC/withPID/ka/neg/prm/pt/numtof"), track.pt(), impParam);
            } else {
              histos.fill(HIST("MC/withPID/ka/neg/prm/pt/numtof"), track.pt(), multiplicity);
            }
            // Matched kaon condition
            if (!(track.mcMask() & (1 << 11))) {
              if (isImpactParam) {
                histos.fill(HIST("MC/withPID/ka/neg/prm/pt/numtof_matched"), track.pt(), impParam);
              } else {
                histos.fill(HIST("MC/withPID/ka/neg/prm/pt/numtof_matched"), track.pt(), multiplicity);
              }
            }
          }
        }
        if (includeCentralityMC) {
          if (isImpactParam) {
            histos.fill(HIST(hpt_numtof_prm[i]), track.pt(), impParam, track.dcaXY(), occupancy);
          } else {
            histos.fill(HIST(hpt_numtof_prm[i]), track.pt(), multiplicity, track.dcaXY(), occupancy);
          }
        } else {
          histos.fill(HIST(hpt_numtof_prm[i]), track.pt(), multiplicity);
        }
        if (!(track.mcMask() & (1 << 11))) {
          if (includeCentralityMC) {
            histos.fill(HIST(hpt_numtofgoodmatch_prm[i]), track.pt(), multiplicity, track.eta()); // RD
          } else {
            histos.fill(HIST(hpt_numtofgoodmatch_prm[i]), track.pt(), multiplicity);
          }
        }
        // Check if the signal is compatible with the PID hypothesis
        float nsigma = 0.f;
        switch (i) {
          case 2:
          case Np + 2:
            nsigma = track.tofNSigmaPi();
            break;
          case 3:
          case Np + 3:
            nsigma = track.tofNSigmaKa();
            break;
          case 4:
          case Np + 4:
            nsigma = track.tofNSigmaPr();
            break;
          default:
            break;
        }
        if (std::abs(nsigma) <= trkselOptions.cfgCutNsigma) {
          if (hPtNumTOFMatchWithPIDSignalPrm[i])
            hPtNumTOFMatchWithPIDSignalPrm[i]->Fill(track.pt(), multiplicity);
        }
      }

      // Filling mismatched info for primary tracks
      if (isMismatchedTrack(track, 0)) {
        histos.fill(HIST(hpt_mism_its_prm[i]), track.pt());
      }
      if (isMismatchedTrack(track, 1)) {
        histos.fill(HIST(hpt_mism_tpc_prm[i]), track.pt());
      }
      if (isMismatchedTrack(track, 2)) {
        histos.fill(HIST(hpt_mism_trd_prm[i]), track.pt());
      }
      if (isMismatchedTrack(track, 3)) {
        histos.fill(HIST(hpt_mism_tof_prm[i]), track.pt());
      }
    }
  }

  template <std::size_t i, typename ParticleType>
  void fillParticleHistograms_MC(const float multiplicity, ParticleType const& mcParticle)
  {
    if (!isParticleEnabled<i>()) { // Check if the particle is enabled
      return;
    }

    if (mcParticle.pdgCode() != PDGs[i]) {
      return;
    }

    if (!mcParticle.isPhysicalPrimary()) {
      if (mcParticle.getProcess() == 4) {
        histos.fill(HIST(hpt_den_str[i]), mcParticle.pt(), multiplicity);
      } else {
        histos.fill(HIST(hpt_den_mat[i]), mcParticle.pt(), multiplicity);
      }
    } else {
      if (includeCentralityMC) {
        histos.fill(HIST(hpt_den_prm[i]), mcParticle.pt(), multiplicity);
      } else {
        histos.fill(HIST(hpt_den_prm[i]), mcParticle.pt(), multiplicity);
      }
    }
  }

  template <std::size_t i, typename ParticleType>
  void fillParticleHistograms_MCRecoEvs(ParticleType const& mcParticle, RecoMCCollisions::iterator const& collision)
  {
    if (!isParticleEnabled<i>()) { // Check if the particle is enabled
      return;
    }

    if (mcParticle.pdgCode() != PDGs[i]) {
      return;
    }

    const auto& mcCollision = collision.mcCollision_as<GenMCCollisions>();
    const float multiplicity = getMultiplicityMC(mcCollision);

    if (mcParticle.isPhysicalPrimary()) {
      if (isEventSelected<false, false>(collision)) {
        if (includeCentralityMC) {
          histos.fill(HIST(hpt_den_prm_goodev[i]), mcParticle.pt(), multiplicity, mcParticle.eta());
        } else {
          histos.fill(HIST(hpt_den_prm_goodev[i]), mcParticle.pt(), multiplicity);
        }
      } else {
        bool isSelected = collision.sel8();
        if (evselOptions.askForCustomTVX) {
          isSelected = collision.selection_bit(aod::evsel::kIsTriggerTVX);
          if (evselOptions.removeITSROFrameBorder) {
            isSelected = isSelected && collision.selection_bit(aod::evsel::kNoITSROFrameBorder);
          }
          if (evselOptions.removeNoSameBunchPileup) {
            isSelected = isSelected && collision.selection_bit(aod::evsel::kNoSameBunchPileup);
          }
          if (evselOptions.requireIsGoodZvtxFT0vsPV) {
            isSelected = isSelected && collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV);
          }
          if (evselOptions.requireIsVertexITSTPC) {
            isSelected = isSelected && collision.selection_bit(aod::evsel::kIsVertexITSTPC);
          }
          if (evselOptions.removeNoTimeFrameBorder) {
            isSelected = isSelected && collision.selection_bit(aod::evsel::kNoTimeFrameBorder);
          }
          if (isSelected) {
            if (includeCentralityMC) {
              histos.fill(HIST(hpt_den_prm_evsel[i]), mcParticle.pt(), multiplicity, mcParticle.eta());
            } else {
              histos.fill(HIST(hpt_den_prm_evsel[i]), mcParticle.pt(), multiplicity);
            }
          }
        } else {
          if (includeCentralityMC) {
            histos.fill(HIST(hpt_den_prm_recoev[i]), mcParticle.pt(), multiplicity, mcParticle.eta());
          } else {
            histos.fill(HIST(hpt_den_prm_recoev[i]), mcParticle.pt(), multiplicity);
          }
        }
      }
    }
  }

  template <std::size_t i, typename ParticleType>
  void fillParticleHistograms_MCGenEvs(ParticleType const& mcParticle, GenMCCollisions::iterator const& mcCollision)
  {

    if (!isParticleEnabled<i>()) { // Check if the particle is enabled
      return;
    }

    if (mcParticle.pdgCode() != PDGs[i]) {
      return;
    }

    const float multiplicity = getMultiplicityMC(mcCollision);
    if (mcParticle.isPhysicalPrimary()) {
      if (std::abs(mcCollision.posZ()) < evselOptions.cfgCutVertex) {
        histos.fill(HIST(hpt_den_prm_mcgoodev[i]), mcParticle.pt(), multiplicity);
      } else {
        histos.fill(HIST(hpt_den_prm_mcbadev[i]), mcParticle.pt(), multiplicity);
      }
    }
  }

  Service<o2::framework::O2DatabasePDG> pdgDB;

  Preslice<aod::McParticles> perMCCol = aod::mcparticle::mcCollisionId;
  SliceCache cache;
  void processMC(soa::Join<aod::Tracks, aod::TracksExtra,
                           aod::TracksDCA, aod::McTrackLabels,
                           aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                           aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr,
                           aod::TrackSelection> const& tracks,
                 aod::McParticles const& mcParticles,
                 GenMCCollisions const& mcCollisions,
                 RecoMCCollisions const& collisions)
  {
    // Fill number of generated and reconstructed collisions for normalization
    histos.fill(HIST("MC/GenRecoCollisions"), 1.f, mcCollisions.size());
    histos.fill(HIST("MC/GenRecoCollisions"), 2.f, collisions.size());

    for (const auto& track : tracks) {
      if (!track.has_collision()) {
        if (track.sign() > 0) {
          histos.fill(HIST("MC/no_collision/pos"), track.pt());
        } else {
          histos.fill(HIST("MC/no_collision/neg"), track.pt());
        }
        continue;
      }
      if (!isEventSelected<false, false>(track.collision_as<RecoMCCollisions>())) {
        continue;
      }
      if (!passesCutWoDCA(track)) {
        continue;
      }
      if (!track.has_mcParticle()) {
        if (track.sign() > 0) {
          histos.fill(HIST("MC/fake/pos"), track.pt());
        } else {
          histos.fill(HIST("MC/fake/neg"), track.pt());
        }
        continue;
      }
      const auto& mcParticle = track.mcParticle();

      static_for<0, 17>([&](auto i) {
        fillTrackHistograms_MC<i>(track, mcParticle, track.collision_as<RecoMCCollisions>(), mcParticles);
      });
    }
    if (includeCentralityMC) {
      for (const auto& collision : collisions) {
        if (!collision.has_mcCollision()) {
          continue;
        }
        const auto& mcCollision = collision.mcCollision_as<GenMCCollisions>();
        const auto& particlesInCollision = mcParticles.sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);
        const float multiplicity = getMultiplicity(collision);
        for (const auto& mcParticle : particlesInCollision) {

          if (std::abs(mcParticle.y()) > trkselOptions.cfgCutY) {
            continue;
          }
          static_for<0, 17>([&](auto i) {
            fillParticleHistograms_MC<i>(multiplicity, mcParticle);
          });
        }
      }
    } else {
      for (const auto& mcParticle : mcParticles) {
        if (std::abs(mcParticle.y()) > trkselOptions.cfgCutY) {
          continue;
        }

        const auto& mcCollision = mcParticle.mcCollision_as<GenMCCollisions>();
        const float multiplicity = getMultiplicityMC(mcCollision);

        static_for<0, 17>([&](auto i) {
          fillParticleHistograms_MC<i>(multiplicity, mcParticle);
        });
      }
    }
    // Loop on reconstructed collisions
    for (const auto& collision : collisions) {
      if (!collision.has_mcCollision()) {
        continue;
      }
      const auto& mcCollision = collision.mcCollision_as<GenMCCollisions>();
      const auto& particlesInCollision = mcParticles.sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);
      if (evselOptions.cfgINELCut.value == 1) {
        if (!o2::pwglf::isINELgt0mc(particlesInCollision, pdgDB)) {
          continue;
        }
      }
      if (evselOptions.cfgINELCut.value == 2) {
        if (!o2::pwglf::isINELgt1mc(particlesInCollision, pdgDB)) {
          continue;
        }
      }

      if (isEventSelected<false, false>(collision)) {
        histos.fill(HIST("MC/MultiplicityRecoEv"), getMultiplicityMC(mcCollision));
      }
      for (const auto& mcParticle : particlesInCollision) {
        if (std::abs(mcParticle.y()) > trkselOptions.cfgCutY) {
          continue;
        }
        static_for<0, 17>([&](auto i) {
          fillParticleHistograms_MCRecoEvs<i>(mcParticle, collision);
        });
      }
    }

    // Loop on generated collisions
    for (const auto& mcCollision : mcCollisions) {
      if (std::abs(mcCollision.posZ()) > evselOptions.cfgCutVertex) {
        continue;
      }
      histos.fill(HIST("MC/Multiplicity"), getMultiplicityMC(mcCollision));
      const auto& particlesInCollision = mcParticles.sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);
      bool hasParticleInFT0C = false;
      bool hasParticleInFT0A = false;
      if (evselOptions.cfgINELCut.value == 1) {
        if (!o2::pwglf::isINELgt0mc(particlesInCollision, pdgDB)) {
          continue;
        }
      }
      histos.fill(HIST("MC/MultiplicityMCINELgt0"), getMultiplicityMC(mcCollision));
      if (evselOptions.cfgINELCut.value == 2) {
        if (!o2::pwglf::isINELgt1mc(particlesInCollision, pdgDB)) {
          continue;
        }
      }
      histos.fill(HIST("MC/MultiplicityMCINELgt1"), getMultiplicityMC(mcCollision));
      for (const auto& mcParticle : particlesInCollision) {
        if (std::abs(mcParticle.y()) > trkselOptions.cfgCutY) {
          continue;
        }
        static_for<0, 17>([&](auto i) {
          fillParticleHistograms_MCGenEvs<i>(mcParticle, mcCollision);
        });
      }
      if (mcCollision.isInelGt0()) {
        histos.fill(HIST("MC/GenRecoCollisions"), 3.f);
      }
      if (mcCollision.isInelGt1()) {
        histos.fill(HIST("MC/GenRecoCollisions"), 4.f);
      }
      if (hasParticleInFT0C && hasParticleInFT0A) {
        histos.fill(HIST("MC/GenRecoCollisions"), 5.f);
      }
    }
  }
  PROCESS_SWITCH(tofSpectra, processMC, "Process MC", false);

  void processMCgen(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles)
  {
    histos.fill(HIST("Vertex/histGenVtxMC"), mcCollision.posZ());
    histos.fill(HIST("Centrality/ImpParm"), mcCollision.impactParameter());
    const float multiplicity = mcCollision.impactParameter();
    for (const auto& mcParticleGen : mcParticles) {
      if (!mcParticleGen.isPhysicalPrimary())
        continue;
      int pdgCode = mcParticleGen.pdgCode();
      float pt = mcParticleGen.pt();
      float absY = std::abs(mcParticleGen.y());
      if (absY > trkselOptions.cfgCutY) {
        continue;
      }

      if (pdgCode == 2212) {
        histos.fill(HIST("MC/test/pr/pos/prm/pt/den"), pt, multiplicity);
      } else if (pdgCode == -2212) {
        histos.fill(HIST("MC/test/pr/neg/prm/pt/den"), pt, multiplicity);
      } else if (pdgCode == 211) {
        histos.fill(HIST("MC/test/pi/pos/prm/pt/den"), pt, multiplicity);
      } else if (pdgCode == -211) {
        histos.fill(HIST("MC/test/pi/neg/prm/pt/den"), pt, multiplicity);
      } else if (pdgCode == 321) {
        histos.fill(HIST("MC/test/ka/pos/prm/pt/den"), pt, multiplicity);
      } else if (pdgCode == -321) {
        histos.fill(HIST("MC/test/ka/neg/prm/pt/den"), pt, multiplicity);
      }
    }
  }
  PROCESS_SWITCH(tofSpectra, processMCgen, "process generated MC", false);
  void processMCgen_RecoEvs(soa::Join<TrackCandidates,
                                      aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                                      aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr> const& tracks,
                            aod::McTrackLabels const& mcTrackLabels,
                            GenMCCollisions const&,
                            RecoMCCollisions const& collisions,
                            aod::McParticles const& mcParticles)
  {
    for (const auto& collision : collisions) {
      if (!collision.has_mcCollision()) {
        continue;
      }

      const auto& mcCollision = collision.mcCollision_as<GenMCCollisions>();
      histos.fill(HIST("Vertex/RecoEvs/histGenVtxMC"), mcCollision.posZ());
      histos.fill(HIST("Centrality/RecoEvs/ImpParm"), mcCollision.impactParameter());
      const float multiplicity = mcCollision.impactParameter();

      for (const auto& track : tracks) {
        if (track.tpcNClsCrossedRows() < 70 ||
            track.tpcChi2NCl() > 4 ||
            track.tpcChi2NCl() < 0.5 ||
            track.itsChi2NCl() > 36 ||
            std::abs(track.dcaXY()) > 0.05 ||
            std::abs(track.dcaZ()) > 2.0 ||
            std::abs(track.eta()) > 0.8 ||
            track.tpcCrossedRowsOverFindableCls() < 0.8 ||
            track.tpcNClsFound() < 100 ||
            !(o2::aod::track::TPCrefit) ||
            !(o2::aod::track::ITSrefit)) {
          continue;
        }
        const auto& mcLabel = mcTrackLabels.iteratorAt(track.globalIndex());
        if (mcLabel.mcParticleId() < 0 || mcLabel.mcParticleId() >= mcParticles.size()) {
          continue;
        }
        const auto& mcParticle = mcParticles.iteratorAt(mcLabel.mcParticleId());
        if (!mcParticle.isPhysicalPrimary()) {
          continue;
        }

        int pdgCode = mcParticle.pdgCode();
        float pt = mcParticle.pt();
        float absY = std::abs(mcParticle.y());
        if (absY > trkselOptions.cfgCutY) {
          continue;
        }

        const auto& nsigmaTPCPi = o2::aod::pidutils::tpcNSigma<2>(track);
        const auto& nsigmaTPCKa = o2::aod::pidutils::tpcNSigma<3>(track);
        const auto& nsigmaTPCPr = o2::aod::pidutils::tpcNSigma<4>(track);

        const bool isPionTPC = std::abs(nsigmaTPCPi) < trkselOptions.cfgCutNsigma;
        const bool isKaonTPC = std::abs(nsigmaTPCKa) < trkselOptions.cfgCutNsigma;
        const bool isProtonTPC = std::abs(nsigmaTPCPr) < trkselOptions.cfgCutNsigma;

        const auto& nsigmaTOFPi = o2::aod::pidutils::tofNSigma<2>(track);
        const auto& nsigmaTOFKa = o2::aod::pidutils::tofNSigma<3>(track);
        const auto& nsigmaTOFPr = o2::aod::pidutils::tofNSigma<4>(track);

        const bool isPionTOF = track.hasTOF() && std::abs(nsigmaTOFPi) < trkselOptions.cfgCutNsigma;
        const bool isKaonTOF = track.hasTOF() && std::abs(nsigmaTOFKa) < trkselOptions.cfgCutNsigma;
        const bool isProtonTOF = track.hasTOF() && std::abs(nsigmaTOFPr) < trkselOptions.cfgCutNsigma;

        if (isPionTPC || isKaonTPC || isProtonTPC) {
          if (pdgCode == 2212) {
            histos.fill(HIST("MC/test/RecoEvs/pr/pos/prm/pt/num"), pt, multiplicity);
          } else if (pdgCode == -2212) {
            histos.fill(HIST("MC/test/RecoEvs/pr/neg/prm/pt/num"), pt, multiplicity);
          } else if (pdgCode == 211) {
            histos.fill(HIST("MC/test/RecoEvs/pi/pos/prm/pt/num"), pt, multiplicity);
          } else if (pdgCode == -211) {
            histos.fill(HIST("MC/test/RecoEvs/pi/neg/prm/pt/num"), pt, multiplicity);
          } else if (pdgCode == 321) {
            histos.fill(HIST("MC/test/RecoEvs/ka/pos/prm/pt/num"), pt, multiplicity);
          } else if (pdgCode == -321) {
            histos.fill(HIST("MC/test/RecoEvs/ka/neg/prm/pt/num"), pt, multiplicity);
          }
        }

        if (isPionTOF || isKaonTOF || isProtonTOF) {
          if (pdgCode == 2212) {
            histos.fill(HIST("MC/test/RecoEvs/pr/pos/prm/pt/numtof"), pt, multiplicity);
          } else if (pdgCode == -2212) {
            histos.fill(HIST("MC/test/RecoEvs/pr/neg/prm/pt/numtof"), pt, multiplicity);
          } else if (pdgCode == 211) {
            histos.fill(HIST("MC/test/RecoEvs/pi/pos/prm/pt/numtof"), pt, multiplicity);
          } else if (pdgCode == -211) {
            histos.fill(HIST("MC/test/RecoEvs/pi/neg/prm/pt/numtof"), pt, multiplicity);
          } else if (pdgCode == 321) {
            histos.fill(HIST("MC/test/RecoEvs/ka/pos/prm/pt/numtof"), pt, multiplicity);
          } else if (pdgCode == -321) {
            histos.fill(HIST("MC/test/RecoEvs/ka/neg/prm/pt/numtof"), pt, multiplicity);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(tofSpectra, processMCgen_RecoEvs, "process generated MC (reconstructed events)", false);
  void processTrackMCLabels(CollisionCandidates::iterator const& collisions,
                            soa::Join<TrackCandidates,
                                      aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                                      aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr> const& tracks,
                            aod::McTrackLabels const& mcTrackLabels, aod::McParticles const& mcParticles)
  {
    const float multiplicity = getMultiplicity(collisions);
    for (const auto& track : tracks) {
      if (!track.has_collision()) {
        continue;
      }
      const auto& collision = track.collision_as<CollisionCandidates>();
      if (!isEventSelected<false, false>(collision)) {
        continue;
      }
      if (!isTrackSelected<true>(track, collision)) {
        continue;
      }
      const auto& mcLabel = mcTrackLabels.iteratorAt(track.globalIndex());
      const auto& mcParticle = mcParticles.iteratorAt(mcLabel.mcParticleId());
      int pdgCode = mcParticle.pdgCode();
      static_for<2, 4>([&](auto par) {
        const auto& nsigmaTPCpar = o2::aod::pidutils::tpcNSigma<par>(track);
        const auto& nsigmaTOFpar = o2::aod::pidutils::tofNSigma<par>(track);
        bool isTPCpar = std::abs(nsigmaTPCpar) < trkselOptions.cfgCutNsigma;
        bool isTOFpar = std::abs(nsigmaTOFpar) < trkselOptions.cfgCutNsigma;
        // Precompute rapidity values to avoid redundant calculations
        double rapiditypar = std::abs(track.rapidity(PID::getMass(par)));
        // TPC Selection and histogram filling
        if (isTPCpar && rapiditypar <= trkselOptions.cfgCutY) {
          static_for<0, 17>([&](auto i) {
            if (pdgCode == PDGs[i]) {
              hMCpdg_nsigmaTPC[par - 2][i]->Fill(track.pt(), nsigmaTPCpar, multiplicity);
              hMCpdg_nsigmaTOF[par - 2][i]->Fill(track.pt(), nsigmaTOFpar, multiplicity);
            }
          });
        }
      });
    }
  }
  PROCESS_SWITCH(tofSpectra, processTrackMCLabels, "Fill track histograms using MC matched PDG labels", false);

}; // end of spectra task

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<tofSpectra>(cfgc)}; }
