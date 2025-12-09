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

#include "PWGLF/DataModel/LFParticleIdentification.h"

#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/McCollisionExtra.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/StaticFor.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "TPDGCode.h"

using namespace o2;
using namespace o2::framework;
using namespace std;
using namespace o2::framework::expressions;
using namespace o2::soa;
using std::array;
using namespace o2::track;
double eta2y(double pt, double m, double eta)
{
  double mt = sqrt(m * m + pt * pt);
  return asinh(pt / mt * sinh(eta));
}

struct QCspectraTPC {
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<int> cfgAverageClusterSize{"cfgAverageClusterSize", 1, "Min Average Cluster Size"};
  Configurable<float> cfgCutnSigmaTPC{"cfgCutnSigmaTPC", 10.0f, "lower bound for nsigmaTPC"};
  Configurable<float> cfgCutnSigmaTPCLow{"cfgCutnSigmaTPCLow", -3.0f, "lower bound for nsigmaTPC"};
  Configurable<float> cfgCutnSigmaTPCHigh{"cfgCutnSigmaTPCHigh", 6.0f, "upper bound for nsigmaTPC"};
  Configurable<float> cfgCutDeltaPLow{"cfgCutDeltaPLow", -2.0f, "lower bound for DeltaP"};
  Configurable<float> cfgCutDeltaPHigh{"cfgCutDeltaPHigh", 2.0f, "upper bound for DeltaP"};
  Configurable<float> cfgCutY{"cfgCutY", 0.5f, "Y range for tracks"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  ConfigurableAxis binsEta{"binsEta", {100, -1, 1}, "Binning of the eta axis"};
  Configurable<float> cfgCutDCAXY{"cfgCutDCAXY", 0.05f, "DCAXY range for tracks"};
  Configurable<float> cfgCutDCAZ{"cfgCutDCAZ", 2.0f, "DCAZ range for tracks"};
  Configurable<int> maxChi2PerClusterTPC{"maxChi2PerClusterTPC", 4, "Additional cut on the maximum value of the chi2 per cluster in the TPC"};
  Configurable<int> maxChi2PerClusterITS{"maxChi2PerClusterITS", 36, "Additional cut on the maximum value of the chi2 per cluster in the ITS"};
  Configurable<int> minTPCNClsFound{"minTPCNClsFound", 70, "Additional cut on the minimum value of the number of found clusters in the TPC"};
  Configurable<float> minNCrossedRowsOverFindableClustersTPC{"minNCrossedRowsOverFindableClustersTPC", 0.8f, "Additional cut on the minimum value of the ratio between crossed rows and findable clusters in the TPC"};
  Configurable<int> multiplicityEstimator{"multiplicityEstimator", 8, "Flag to use a multiplicity estimator: 0 no multiplicity, 1 MultFV0M, 2 MultFT0M, 3 MultFDDM, 4 MultTracklets, 5 MultTPC, 6 MultNTracksPV, 7 MultNTracksPVeta1, 8 CentralityFT0C, 9 CentralityFT0M, 10 CentralityFV0A"};
  Configurable<int> minNCrossedRowsTPC{"minNCrossedRowsTPC", 100, "Additional cut on the minimum number of crossed rows in the TPC"};
  Configurable<int> ITSNCls{"ITSNCls", 7, "Additional cut on the minimum number of ITS clusters"};
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0}, "Binning of the pT axis"};
  ConfigurableAxis binsDca{"binsDca", {VARIABLE_WIDTH, -3.0, -2.95, -2.9, -2.85, -2.8, -2.75, -2.7, -2.65, -2.6, -2.55, -2.5, -2.45, -2.4, -2.35, -2.3, -2.25, -2.2, -2.15, -2.1, -2.05, -2.0, -1.975, -1.95, -1.925, -1.9, -1.875, -1.85, -1.825, -1.8, -1.775, -1.75, -1.725, -1.7, -1.675, -1.65, -1.625, -1.6, -1.575, -1.55, -1.525, -1.5, -1.475, -1.45, -1.425, -1.4, -1.375, -1.35, -1.325, -1.3, -1.275, -1.25, -1.225, -1.2, -1.175, -1.15, -1.125, -1.1, -1.075, -1.05, -1.025, -1.0, -0.99, -0.98, -0.97, -0.96, -0.95, -0.94, -0.93, -0.92, -0.91, -0.9, -0.89, -0.88, -0.87, -0.86, -0.85, -0.84, -0.83, -0.82, -0.81, -0.8, -0.79, -0.78, -0.77, -0.76, -0.75, -0.74, -0.73, -0.72, -0.71, -0.7, -0.69, -0.68, -0.67, -0.66, -0.65, -0.64, -0.63, -0.62, -0.61, -0.6, -0.59, -0.58, -0.57, -0.56, -0.55, -0.54, -0.53, -0.52, -0.51, -0.5, -0.49, -0.48, -0.47, -0.46, -0.45, -0.44, -0.43, -0.42, -0.41, -0.4, -0.396, -0.392, -0.388, -0.384, -0.38, -0.376, -0.372, -0.368, -0.364, -0.36, -0.356, -0.352, -0.348, -0.344, -0.34, -0.336, -0.332, -0.328, -0.324, -0.32, -0.316, -0.312, -0.308, -0.304, -0.3, -0.296, -0.292, -0.288, -0.284, -0.28, -0.276, -0.272, -0.268, -0.264, -0.26, -0.256, -0.252, -0.248, -0.244, -0.24, -0.236, -0.232, -0.228, -0.224, -0.22, -0.216, -0.212, -0.208, -0.204, -0.2, -0.198, -0.196, -0.194, -0.192, -0.19, -0.188, -0.186, -0.184, -0.182, -0.18, -0.178, -0.176, -0.174, -0.172, -0.17, -0.168, -0.166, -0.164, -0.162, -0.16, -0.158, -0.156, -0.154, -0.152, -0.15, -0.148, -0.146, -0.144, -0.142, -0.14, -0.138, -0.136, -0.134, -0.132, -0.13, -0.128, -0.126, -0.124, -0.122, -0.12, -0.118, -0.116, -0.114, -0.112, -0.11, -0.108, -0.106, -0.104, -0.102, -0.1, -0.099, -0.098, -0.097, -0.096, -0.095, -0.094, -0.093, -0.092, -0.091, -0.09, -0.089, -0.088, -0.087, -0.086, -0.085, -0.084, -0.083, -0.082, -0.081, -0.08, -0.079, -0.078, -0.077, -0.076, -0.075, -0.074, -0.073, -0.072, -0.071, -0.07, -0.069, -0.068, -0.067, -0.066, -0.065, -0.064, -0.063, -0.062, -0.061, -0.06, -0.059, -0.058, -0.057, -0.056, -0.055, -0.054, -0.053, -0.052, -0.051, -0.05, -0.049, -0.048, -0.047, -0.046, -0.045, -0.044, -0.043, -0.042, -0.041, -0.04, -0.039, -0.038, -0.037, -0.036, -0.035, -0.034, -0.033, -0.032, -0.031, -0.03, -0.029, -0.028, -0.027, -0.026, -0.025, -0.024, -0.023, -0.022, -0.021, -0.02, -0.019, -0.018, -0.017, -0.016, -0.015, -0.014, -0.013, -0.012, -0.011, -0.01, -0.009, -0.008, -0.007, -0.006, -0.005, -0.004, -0.003, -0.002, -0.001, -0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.017, 0.018, 0.019, 0.02, 0.021, 0.022, 0.023, 0.024, 0.025, 0.026, 0.027, 0.028, 0.029, 0.03, 0.031, 0.032, 0.033, 0.034, 0.035, 0.036, 0.037, 0.038, 0.039, 0.04, 0.041, 0.042, 0.043, 0.044, 0.045, 0.046, 0.047, 0.048, 0.049, 0.05, 0.051, 0.052, 0.053, 0.054, 0.055, 0.056, 0.057, 0.058, 0.059, 0.06, 0.061, 0.062, 0.063, 0.064, 0.065, 0.066, 0.067, 0.068, 0.069, 0.07, 0.071, 0.072, 0.073, 0.074, 0.075, 0.076, 0.077, 0.078, 0.079, 0.08, 0.081, 0.082, 0.083, 0.084, 0.085, 0.086, 0.087, 0.088, 0.089, 0.09, 0.091, 0.092, 0.093, 0.094, 0.095, 0.096, 0.097, 0.098, 0.099, 0.1, 0.102, 0.104, 0.106, 0.108, 0.11, 0.112, 0.114, 0.116, 0.118, 0.12, 0.122, 0.124, 0.126, 0.128, 0.13, 0.132, 0.134, 0.136, 0.138, 0.14, 0.142, 0.144, 0.146, 0.148, 0.15, 0.152, 0.154, 0.156, 0.158, 0.16, 0.162, 0.164, 0.166, 0.168, 0.17, 0.172, 0.174, 0.176, 0.178, 0.18, 0.182, 0.184, 0.186, 0.188, 0.19, 0.192, 0.194, 0.196, 0.198, 0.2, 0.204, 0.208, 0.212, 0.216, 0.22, 0.224, 0.228, 0.232, 0.236, 0.24, 0.244, 0.248, 0.252, 0.256, 0.26, 0.264, 0.268, 0.272, 0.276, 0.28, 0.284, 0.288, 0.292, 0.296, 0.3, 0.304, 0.308, 0.312, 0.316, 0.32, 0.324, 0.328, 0.332, 0.336, 0.34, 0.344, 0.348, 0.352, 0.356, 0.36, 0.364, 0.368, 0.372, 0.376, 0.38, 0.384, 0.388, 0.392, 0.396, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.0, 1.025, 1.05, 1.075, 1.1, 1.125, 1.15, 1.175, 1.2, 1.225, 1.25, 1.275, 1.3, 1.325, 1.35, 1.375, 1.4, 1.425, 1.45, 1.475, 1.5, 1.525, 1.55, 1.575, 1.6, 1.625, 1.65, 1.675, 1.7, 1.725, 1.75, 1.775, 1.8, 1.825, 1.85, 1.875, 1.9, 1.925, 1.95, 1.975, 2.0, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.85, 2.9, 2.95, 3.0}, "Binning of DCA xy and z axis"};
  Configurable<bool> doProcessMC{"doProcessMC", true, "choose if produce MC"};
  ConfigurableAxis binsITSCusterSize{"binsITSCusterSize", {50, 0., 10.}, "Binning of the ITS Cluster axis"};
  ConfigurableAxis binsPGlobal{"binsPGlobal", {300, -3, 3}, "Binning of the global p axis"};
  ConfigurableAxis binsPTPC{"binsPTPC", {300, -3, 3}, "Binning of the TPC p axis"};
  ConfigurableAxis binsDeltaP{"binsDeltaP", {600, -3, 3}, "Binning of the #delta P axis"};
  ConfigurableAxis binsnsigmaTPC{"binsnsigmaTPC", {200, -10, 10}, "Binning of the n_{#sigma, TPC} axis"};
  ConfigurableAxis binsMultiplicity{"binsMultiplicity", {100, 0, 100}, "Binning for multiplicity"};
  // Histograms
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  void init(o2::framework::InitContext&)
  {
    const AxisSpec axisEta{30, -1.5, +1.5, "#eta"};
    const AxisSpec axisClusterSize{binsITSCusterSize, "ITS cluster size * cos(#lambda)"};
    const AxisSpec axisY{30, -1.5, +1.5, "y"};
    const AxisSpec axisPGlobal{binsPGlobal, "p_{Global}/z"};
    const AxisSpec axisPTPC{binsPTPC, "p_{TPC}/z"};
    const AxisSpec axisDeltaP{binsDeltaP, "p_{TPC}-p_{Global}/z"};
    const AxisSpec axisnSigmaTPC{binsnsigmaTPC, "nsigmaTPC"};
    const AxisSpec pAxis{binsPt, "#it{p_{TPC}} (GeV/#it{c})"};
    const AxisSpec ptAxis{binsPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec dcaXyAxis{binsDca, "DCA_{xy} (cm)"};
    const AxisSpec phiAxis{200, 0, 7, "#it{#varphi} (rad)"};
    const AxisSpec dcaZAxis{binsDca, "DCA_{z} (cm)"};
    const AxisSpec multAxis{binsMultiplicity, "multiplicity estimator"};
    histos.add("event/vertexz", "collision z position", HistType::kTH1F, {{100, -10., +10., "z position (cm)"}});
    histos.add("etaHistogram", "etaHistogram", kTH1F, {axisEta});
    histos.add("rapidityHistogram", "rapidityHistogram", kTH1F, {axisY});
    histos.add("ptHistogram", "ptHistogram", kTH1F, {ptAxis});
    histos.add("pHistogram", "pHistogram", kTH1F, {pAxis});
    histos.add("histTpcSignal", "Specific energy loss", HistType::kTH2F, {{600, -6., 6., "#it{p/z} (GeV/#it{c})"}, {1400, 0, 1400, "d#it{E} / d#it{X} (a. u.)"}});
    histos.add("histTpcSignal_Pr", "Specific energy loss of Proton", HistType::kTH2F, {{600, -6., 6., "#it{p/z} (GeV/#it{c})"}, {1400, 0, 1400, "d#it{E} / d#it{X} (a. u.)"}});
    histos.add("histTpcNsigmaPr", "n-sigmaPr TPC", HistType::kTH2F, {ptAxis, {200, -100., +100., "n#sigma_Pr} (a. u.)"}});

    histos.add("nsigmatpc/clusterSize/histITSClusterSize_Pr_layer0", "ITS cluster, p_{TPC}", HistType::kTH2F, {pAxis, {axisClusterSize}});
    histos.add("nsigmatpc/clusterSize/histITSClusterSize_Pr_layer1", "ITS cluster, p_{TPC}", HistType::kTH2F, {pAxis, {axisClusterSize}});
    histos.add("nsigmatpc/clusterSize/histITSClusterSize_Pr_layer2", "ITS cluster, p_{TPC}", HistType::kTH2F, {pAxis, {axisClusterSize}});
    histos.add("nsigmatpc/clusterSize/histITSClusterSize_Pr_layer3", "ITS cluster, p_{TPC}", HistType::kTH2F, {pAxis, {axisClusterSize}});
    histos.add("nsigmatpc/clusterSize/histITSClusterSize_Pr_layer4", "ITS cluster, p_{TPC}", HistType::kTH2F, {pAxis, {axisClusterSize}});
    histos.add("nsigmatpc/clusterSize/histITSClusterSize_Pr_layer5", "ITS cluster, p_{TPC}", HistType::kTH2F, {pAxis, {axisClusterSize}});
    histos.add("nsigmatpc/clusterSize/histITSClusterSize_Pr_layer6", "ITS cluster, p_{TPC}", HistType::kTH2F, {pAxis, {axisClusterSize}});
    histos.add("nsigmatpc/clusterSize/histITSClusterSize_Pr", "ITS cluster, p_{TPC},nsigma", HistType::kTH3F, {ptAxis, {axisClusterSize}, {axisnSigmaTPC}});

    histos.add("nsigmatpc/clusterSize/p/histITSClusterSize_Pr_Alllayer", "ITS cluster, p_{TPC}", HistType::kTH2F, {{axisPTPC}, {axisClusterSize}});
    histos.add("nsigmatpc/clusterSize/pt/histITSClusterSize_Pr_Alllayer", "ITS cluster, p_{TPC}", HistType::kTH2F, {{axisPTPC}, {axisClusterSize}});

    histos.add("nsigmatpcCut/clusterSize/histITSClusterSize_Pr_layer0", "ITS cluster, p_{TPC}", HistType::kTH2F, {pAxis, {axisClusterSize}});
    histos.add("nsigmatpcCut/clusterSize/histITSClusterSize_Pr_layer1", "ITS cluster, p_{TPC}", HistType::kTH2F, {pAxis, {axisClusterSize}});
    histos.add("nsigmatpcCut/clusterSize/histITSClusterSize_Pr_layer2", "ITS cluster, p_{TPC}", HistType::kTH2F, {pAxis, {axisClusterSize}});
    histos.add("nsigmatpcCut/clusterSize/histITSClusterSize_Pr_layer3", "ITS cluster, p_{TPC}", HistType::kTH2F, {pAxis, {axisClusterSize}});
    histos.add("nsigmatpcCut/clusterSize/histITSClusterSize_Pr_layer4", "ITS cluster, p_{TPC}", HistType::kTH2F, {pAxis, {axisClusterSize}});
    histos.add("nsigmatpcCut/clusterSize/histITSClusterSize_Pr_layer5", "ITS cluster, p_{TPC}", HistType::kTH2F, {pAxis, {axisClusterSize}});
    histos.add("nsigmatpcCut/clusterSize/histITSClusterSize_Pr_layer6", "ITS cluster, p_{TPC}", HistType::kTH2F, {pAxis, {axisClusterSize}});

    histos.add("nsigmatpcCut/clusterSize/p/histITSClusterSize_Pr_Alllayer", "ITS cluster, p_{TPC}", HistType::kTH2F, {{axisPTPC}, {axisClusterSize}});
    histos.add("nsigmatpcCut/clusterSize/pt/histITSClusterSize_Pr_Alllayer", "ITS cluster, p_{TPC}", HistType::kTH2F, {{axisPTPC}, {axisClusterSize}});
    histos.add("nsigmatpcCut/clusterSize/histITSClusterSize_Pr", "ITS cluster, p_{TPC},nsigma", HistType::kTH3F, {ptAxis, {axisClusterSize}, {axisnSigmaTPC}});
    histos.add("nsigmatpcCut/deltapCut/clusterSize/histITSClusterSize_Pr", "ITS cluster, p_{TPC},nsigma", HistType::kTH3F, {ptAxis, {axisClusterSize}, {axisnSigmaTPC}});

    histos.add("nsigmatpc/pr/pos", "n_#sigma TPC proton", HistType::kTH2F, {ptAxis, {axisnSigmaTPC}});
    histos.add("nsigmatpc/pr/neg", "n_#sigma TPC anti-proton", HistType::kTH2F, {ptAxis, {axisnSigmaTPC}});
    histos.add("nsigmatpcCut/pr/pos", "n_#sigma TPC proton", HistType::kTH2F, {ptAxis, {axisnSigmaTPC}});
    histos.add("nsigmatpcCut/pr/neg", "n_#sigma TPC anti-proton", HistType::kTH2F, {ptAxis, {axisnSigmaTPC}});
    histos.add("delta_p", "delta_p", HistType::kTH2F, {{axisPGlobal}, {axisDeltaP}});

    //**************************************************************************************************
    // MC hists
    if (doProcessMC) {
      histos.add("spectraGen/histGenVetxZ", "PosZ generated events", HistType::kTH1F, {{1500, -15.f, 15.f, "Vertex Z (cm)"}});
      histos.add("spectraGen/proton/pos/histGenPtProton", "generated particles", HistType::kTH1F, {ptAxis});
      histos.add("spectraGen/proton/pos/histGenPtProtonPrim", "generated particles", HistType::kTH1F, {ptAxis});
      histos.add("spectraGen/proton/pos/histGenPtProtonSec", "generated particles", HistType::kTH1F, {ptAxis});
      histos.add("spectraRec/proton/pos/histRecPtProton", "recosntructed particles", HistType::kTH1F, {ptAxis});
      histos.add("spectraRec/proton/pos/histRecPtProtonPrim", "recosntructed particles", HistType::kTH1F, {ptAxis});
      histos.add("MC/pi/pos/prm/pt/num", "recosntructed MC #pi^{+}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/ka/pos/prm/pt/num", "recosntructed MC k^{+}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/pr/pos/prm/pt/num", "recosntructed MC p", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/pi/neg/prm/pt/num", "recosntructed MC #pi^{-}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/ka/neg/prm/pt/num", "recosntructed MC k^{-}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/pr/neg/prm/pt/num", "recosntructed MC #bar{p}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/pi/pos/str/pt/num", "recosntructed MC #pi^{+}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/pi/neg/str/pt/num", "recosntructed MC #pi^{-}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/ka/pos/str/pt/num", "recosntructed MC k^{+}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/ka/neg/str/pt/num", "recosntructed MC k^{-}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/pr/pos/str/pt/num", "recosntructed MC p", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/pr/neg/str/pt/num", "recosntructed MC #bar{p}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/pi/pos/mat/pt/num", "recosntructed MC #pi^{+}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/pi/neg/mat/pt/num", "recosntructed MC #pi^{-}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/ka/pos/mat/pt/num", "recosntructed MC k^{+}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/ka/neg/mat/pt/num", "recosntructed MC k^{-}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/pr/pos/mat/pt/num", "recosntructed MC p", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/pr/neg/mat/pt/num", "recosntructed MC #bar{p}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/pi/pos/prm/pt/numtof", "recosntructed MC TOF #pi^{+}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/ka/pos/prm/pt/numtof", "recosntructed MC TOF k^{+}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/pr/pos/prm/pt/numtof", "recosntructed MC TOF p", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/pi/neg/prm/pt/numtof", "recosntructed MC TOF #pi^{-}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/ka/neg/prm/pt/numtof", "recosntructed MC TOF k^{-}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/pr/neg/prm/pt/numtof", "recosntructed MC TOF #bar{p}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/pi/pos/str/pt/numtof", "recosntructed MC TOF #pi^{+}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/pi/neg/str/pt/numtof", "recosntructed MC TOF #pi^{-}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/ka/pos/str/pt/numtof", "recosntructed MC TOF k^{+}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/ka/neg/str/pt/numtof", "recosntructed MC TOF k^{-}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/pr/pos/str/pt/numtof", "recosntructed MC TOF p", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/pr/neg/str/pt/numtof", "recosntructed MC TOF #bar{p}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/pi/pos/mat/pt/numtof", "recosntructed MC TOF #pi^{+}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/pi/neg/mat/pt/numtof", "recosntructed MC TOF #pi^{-}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/ka/pos/mat/pt/numtof", "recosntructed MC TOF k^{+}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/ka/neg/mat/pt/numtof", "recosntructed MC TOF k^{-}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/pr/pos/mat/pt/numtof", "recosntructed MC TOF p", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/pr/neg/mat/pt/numtof", "recosntructed MC TOF #bar{p}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/pi/pos/prm/pt/den", "generated MC #pi^{+}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/pi/neg/prm/pt/den", "generated MC #pi^{-}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/ka/pos/prm/pt/den", "generated MC K^{+}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/ka/neg/prm/pt/den", "generated MC K^{-}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/pr/pos/prm/pt/den", "generated MC p", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/pr/neg/prm/pt/den", "generated MC #bar{p}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/pi/pos/str/pt/den", "generated MC #pi^{+}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/pi/neg/str/pt/den", "generated MC #pi^{-}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/ka/pos/str/pt/den", "generated MC K^{+}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/ka/neg/str/pt/den", "generated MC K^{-}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/pr/pos/str/pt/den", "generated MC p", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/pr/neg/str/pt/den", "generated MC #bar{p}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/pi/pos/mat/pt/den", "generated MC #pi^{+}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/pi/neg/mat/pt/den", "generated MC #pi^{-}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/ka/pos/mat/pt/den", "generated MC K^{+}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/ka/neg/mat/pt/den", "generated MC K^{-}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/pr/pos/mat/pt/den", "generated MC p", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/pr/neg/mat/pt/den", "generated MC #bar{p}", HistType::kTHnSparseD, {ptAxis, multAxis, dcaXyAxis});
      histos.add("MC/test/pi/pos/prm/pt/den", "generated MC #pi^{+}", HistType::kTHnSparseD, {ptAxis});
      histos.add("MC/test/pi/neg/prm/pt/den", "generated MC #pi^{-}", HistType::kTHnSparseD, {ptAxis});
      histos.add("MC/test/ka/pos/prm/pt/den", "generated MC K^{+}", HistType::kTHnSparseD, {ptAxis});
      histos.add("MC/test/ka/neg/prm/pt/den", "generated MC K^{-}", HistType::kTHnSparseD, {ptAxis});
      histos.add("MC/test/pr/pos/prm/pt/den", "generated MC p", HistType::kTHnSparseD, {ptAxis});
      histos.add("MC/test/pr/neg/prm/pt/den", "generated MC #bar{p}", HistType::kTHnSparseD, {ptAxis});
    }
  } // init loop end

  using CollisionCandidate = soa::Join<aod::Collisions, aod::EvSels>;
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPr>;
  void process(CollisionCandidate::iterator const& collision, TrackCandidates const& tracks)
  {
    if (!collision.sel8()) {
      return;
    }

    if (abs(collision.posZ()) > cfgCutVertex) {
      return;
    }

    histos.fill(HIST("event/vertexz"), collision.posZ());

    for (auto& track : tracks) {

      if (abs(track.eta()) > cfgCutEta) {
        return;
      }

      auto y_proton = eta2y(track.pt(), 0.938, track.eta());
      //  auto y_deutron = eta2y(track.pt(),1.875, track.eta());

      if (abs(y_proton) > cfgCutY) {
        return;
      }

      if (track.itsNCls() != ITSNCls)
        continue;
      if (track.tpcNClsFound() < 100)
        continue;
      if (track.tpcNClsCrossedRows() < 100)
        continue;
      if (track.tpcChi2NCl() > maxChi2PerClusterTPC)
        continue;
      if (track.itsChi2NCl() > 36)
        continue;
      if (abs(track.dcaXY()) > cfgCutDCAXY)
        continue;
      if (abs(track.dcaZ()) > cfgCutDCAZ)
        continue;

      histos.fill(HIST("etaHistogram"), track.eta());
      histos.fill(HIST("rapidityHistogram"), y_proton);
      histos.fill(HIST("ptHistogram"), track.pt());
      histos.fill(HIST("pHistogram"), track.tpcInnerParam());

      auto nSigmaTPCPr = o2::aod::pidutils::tpcNSigma<PID::Proton>(track);

      uint32_t clsizeflag = track.itsClusterSizes();

      auto clSizeLayer0 = (clsizeflag >> (0 * 4)) & 0xf;
      auto clSizeLayer1 = (clsizeflag >> (1 * 4)) & 0xf;
      auto clSizeLayer2 = (clsizeflag >> (2 * 4)) & 0xf;
      auto clSizeLayer3 = (clsizeflag >> (3 * 4)) & 0xf;
      auto clSizeLayer4 = (clsizeflag >> (4 * 4)) & 0xf;
      auto clSizeLayer5 = (clsizeflag >> (5 * 4)) & 0xf;
      auto clSizeLayer6 = (clsizeflag >> (6 * 4)) & 0xf;

      int numLayers = 7; // Number of layers
      int sumClusterSizes = clSizeLayer1 + clSizeLayer2 + clSizeLayer3 + clSizeLayer4 + clSizeLayer5 + clSizeLayer6 + clSizeLayer0;
      double cos_lambda = std::cos(std::atan(track.tgl()));
      double averageClusterSize = (static_cast<double>(sumClusterSizes) / numLayers) * cos_lambda;

      auto delta_p = track.tpcInnerParam() - track.p();

      histos.fill(HIST("histTpcSignal"), track.tpcInnerParam() * track.sign(), track.tpcSignal());
      histos.fill(HIST("histTpcNsigmaPr"), track.tpcInnerParam(), nSigmaTPCPr);
      histos.fill(HIST("delta_p"), track.p() * track.sign(), delta_p * track.sign());
      //  return;

      if (averageClusterSize < cfgAverageClusterSize)
        continue;

      if (std::abs(nSigmaTPCPr) < cfgCutnSigmaTPC) { // put nsigma cut
        histos.fill(HIST("histTpcSignal_Pr"), track.tpcInnerParam() * track.sign(), track.tpcSignal());

        histos.fill(HIST("nsigmatpc/clusterSize/histITSClusterSize_Pr"), track.pt(), averageClusterSize, nSigmaTPCPr);
        histos.fill(HIST("nsigmatpc/clusterSize/histITSClusterSize_Pr_layer0"), track.tpcInnerParam(), clSizeLayer0 * cos_lambda);
        histos.fill(HIST("nsigmatpc/clusterSize/histITSClusterSize_Pr_layer1"), track.tpcInnerParam(), clSizeLayer1 * cos_lambda);
        histos.fill(HIST("nsigmatpc/clusterSize/histITSClusterSize_Pr_layer2"), track.tpcInnerParam(), clSizeLayer2 * cos_lambda);
        histos.fill(HIST("nsigmatpc/clusterSize/histITSClusterSize_Pr_layer3"), track.tpcInnerParam(), clSizeLayer3 * cos_lambda);
        histos.fill(HIST("nsigmatpc/clusterSize/histITSClusterSize_Pr_layer4"), track.tpcInnerParam(), clSizeLayer4 * cos_lambda);
        histos.fill(HIST("nsigmatpc/clusterSize/histITSClusterSize_Pr_layer5"), track.tpcInnerParam(), clSizeLayer5 * cos_lambda);
        histos.fill(HIST("nsigmatpc/clusterSize/histITSClusterSize_Pr_layer6"), track.tpcInnerParam(), clSizeLayer6 * cos_lambda);

        histos.fill(HIST("nsigmatpc/clusterSize/p/histITSClusterSize_Pr_Alllayer"), track.tpcInnerParam() * track.sign(), averageClusterSize);

        histos.fill(HIST("nsigmatpc/clusterSize/pt/histITSClusterSize_Pr_Alllayer"), track.pt() * track.sign(), averageClusterSize);

        if (track.sign() > 0) {
          histos.fill(HIST("nsigmatpc/pr/pos"), track.pt(), nSigmaTPCPr);
        } else {
          histos.fill(HIST("nsigmatpc/pr/neg"), track.pt(), nSigmaTPCPr);
        }
      }

      if (cfgCutnSigmaTPCLow < nSigmaTPCPr && nSigmaTPCPr < cfgCutnSigmaTPCHigh) {
        if (track.sign() > 0) {
          histos.fill(HIST("nsigmatpcCut/pr/pos"), track.pt(), nSigmaTPCPr);
        } else {
          histos.fill(HIST("nsigmatpcCut/pr/neg"), track.pt(), nSigmaTPCPr);
        }

        histos.fill(HIST("nsigmatpcCut/clusterSize/histITSClusterSize_Pr"), track.pt(), averageClusterSize, nSigmaTPCPr);
        histos.fill(HIST("nsigmatpcCut/clusterSize/histITSClusterSize_Pr_layer0"), track.tpcInnerParam(), clSizeLayer0 * cos_lambda);
        histos.fill(HIST("nsigmatpcCut/clusterSize/histITSClusterSize_Pr_layer1"), track.tpcInnerParam(), clSizeLayer1 * cos_lambda);
        histos.fill(HIST("nsigmatpcCut/clusterSize/histITSClusterSize_Pr_layer2"), track.tpcInnerParam(), clSizeLayer2 * cos_lambda);
        histos.fill(HIST("nsigmatpcCut/clusterSize/histITSClusterSize_Pr_layer3"), track.tpcInnerParam(), clSizeLayer3 * cos_lambda);
        histos.fill(HIST("nsigmatpcCut/clusterSize/histITSClusterSize_Pr_layer4"), track.tpcInnerParam(), clSizeLayer4 * cos_lambda);
        histos.fill(HIST("nsigmatpcCut/clusterSize/histITSClusterSize_Pr_layer5"), track.tpcInnerParam(), clSizeLayer5 * cos_lambda);
        histos.fill(HIST("nsigmatpcCut/clusterSize/histITSClusterSize_Pr_layer6"), track.tpcInnerParam(), clSizeLayer6 * cos_lambda);

        histos.fill(HIST("nsigmatpcCut/clusterSize/p/histITSClusterSize_Pr_Alllayer"), track.tpcInnerParam() * track.sign(), averageClusterSize);

        histos.fill(HIST("nsigmatpcCut/clusterSize/pt/histITSClusterSize_Pr_Alllayer"), track.pt() * track.sign(), averageClusterSize);

        if (cfgCutDeltaPLow < delta_p && delta_p < cfgCutDeltaPHigh) {

          histos.fill(HIST("nsigmatpcCut/deltapCut/clusterSize/histITSClusterSize_Pr"), track.pt(), averageClusterSize, nSigmaTPCPr);
        }
      }

    } // tracks loop

  } // collision loop

  PROCESS_SWITCH(QCspectraTPC, process, "process data", true);

  void processMCGen(aod::McCollision const& mcCollision, aod::McParticles& mcParticles)
  {
    if (doProcessMC) {
      histos.fill(HIST("spectraGen/histGenVetxZ"), mcCollision.posZ());
    }

    for (auto& mcParticleGen : mcParticles) {
      if (abs(mcParticleGen.y()) > std::abs(cfgCutY)) {
        continue;
      }
      bool isPhysPrim = mcParticleGen.isPhysicalPrimary();
      if (doProcessMC) {
        if (mcParticleGen.pdgCode() == 2212) {
          histos.fill(HIST("spectraGen/proton/pos/histGenPtProton"), mcParticleGen.pt());
          if (isPhysPrim) {
            histos.fill(HIST("spectraGen/proton/pos/histGenPtProtonPrim"), mcParticleGen.pt());
          } else {
            histos.fill(HIST("spectraGen/proton/pos/histGenPtProtonSec"), mcParticleGen.pt());
          }
        }
      }
    }
  } // process_mc loop end

  PROCESS_SWITCH(QCspectraTPC, processMCGen, "process MC Generated", true);

  using CollisionCandidateMC = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels>; // RD
  void processMCRec(soa::Join<aod::Tracks, aod::TracksExtra,
                              aod::TracksDCA, aod::McTrackLabels,
                              aod::TrackSelection> const& tracks,
                    aod::McParticles const& /*mcParticles*/,
                    aod::McCollisions const& /*mcCollisions*/,
                    CollisionCandidateMC const& /*collisions*/)
  {

    for (auto& track : tracks) {

      const auto& mcParticle = track.mcParticle();
      if (mcParticle.pdgCode() != 2212) {
        return;
      }

      if (abs(track.eta()) > cfgCutEta) {
        return;
      }

      if (std::abs(mcParticle.y()) > cfgCutY) {
        return;
      }

      uint32_t clsizeflag = track.itsClusterSizes();

      auto clSizeLayer0 = (clsizeflag >> (0 * 4)) & 0xf;
      auto clSizeLayer1 = (clsizeflag >> (1 * 4)) & 0xf;
      auto clSizeLayer2 = (clsizeflag >> (2 * 4)) & 0xf;
      auto clSizeLayer3 = (clsizeflag >> (3 * 4)) & 0xf;
      auto clSizeLayer4 = (clsizeflag >> (4 * 4)) & 0xf;
      auto clSizeLayer5 = (clsizeflag >> (5 * 4)) & 0xf;
      auto clSizeLayer6 = (clsizeflag >> (6 * 4)) & 0xf;

      int numLayers = 7; // Number of layers
      int sumClusterSizes = clSizeLayer1 + clSizeLayer2 + clSizeLayer3 + clSizeLayer4 + clSizeLayer5 + clSizeLayer6 + clSizeLayer0;
      double cos_lambda = std::cos(std::atan(track.tgl()));
      double averageClusterSize = (static_cast<double>(sumClusterSizes) / numLayers) * cos_lambda;
      if (averageClusterSize < cfgAverageClusterSize)
        continue;

      if (track.itsNCls() != ITSNCls)
        continue;
      if (track.tpcNClsFound() < 100)
        continue;
      if (track.tpcNClsCrossedRows() < 100)
        continue;
      if (track.tpcChi2NCl() > 4)
        continue;
      if (track.itsChi2NCl() > 36)
        continue;
      if (abs(track.dcaXY()) > cfgCutDCAXY)
        continue;
      if (abs(track.dcaZ()) > cfgCutDCAZ)
        continue;

      histos.fill(HIST("spectraRec/proton/pos/histRecPtProton"), track.pt());
      if (mcParticle.isPhysicalPrimary()) {
        histos.fill(HIST("spectraRec/proton/pos/histRecPtProtonPrim"), track.pt());
      }
    }
  } // process_mc loop end
  using CollisionCandidateMCRec = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0Cs>;
  void processTrackHistograms_MC(CollisionCandidateMCRec::iterator const& collision, soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::TrackSelection> const& tracks, aod::McParticles const& /*mcParticles*/, aod::McCollisions const& /*mcCollisions*/)
  {
    const float multiplicity = collision.centFT0C();
    for (auto& track : tracks) {

      const auto& mcParticle = track.mcParticle();
      if (abs(track.eta()) > cfgCutEta) {
        return;
      }
      if (track.tpcNClsCrossedRows() < minTPCNClsFound)
        continue;
      if (track.tpcCrossedRowsOverFindableCls() < minNCrossedRowsOverFindableClustersTPC)
        continue;
      if (track.tpcChi2NCl() > maxChi2PerClusterTPC)
        continue;
      if (track.itsChi2NCl() > maxChi2PerClusterITS)
        continue;
      if (abs(track.dcaXY()) > cfgCutDCAXY)
        continue;
      if (abs(track.dcaZ()) > cfgCutDCAZ)
        continue;
      if (mcParticle.isPhysicalPrimary()) {
        // pions
        if (mcParticle.pdgCode() == 211) {
          if (std::abs(mcParticle.y()) > cfgCutY) {
            return;
          }
          histos.fill(HIST("MC/pi/pos/prm/pt/num"), track.pt(), multiplicity, track.dcaXY());
          if (track.hasTOF()) {
            histos.fill(HIST("MC/pi/pos/prm/pt/numtof"), track.pt(), multiplicity, track.dcaXY());
          }
        }
        if (mcParticle.pdgCode() == -211) {
          if (std::abs(mcParticle.y()) > cfgCutY) {
            return;
          }
          histos.fill(HIST("MC/pi/neg/prm/pt/num"), track.pt(), multiplicity, track.dcaXY());
          if (track.hasTOF()) {
            histos.fill(HIST("MC/pi/neg/prm/pt/numtof"), track.pt(), multiplicity, track.dcaXY());
          }
        }
        // kaons
        if (mcParticle.pdgCode() == 321) {
          if (std::abs(mcParticle.y()) > cfgCutY) {
            return;
          }
          histos.fill(HIST("MC/ka/pos/prm/pt/num"), track.pt(), multiplicity, track.dcaXY());
          if (track.hasTOF()) {
            histos.fill(HIST("MC/ka/pos/prm/pt/numtof"), track.pt(), multiplicity, track.dcaXY());
          }
        }
        if (mcParticle.pdgCode() == -321) {
          if (std::abs(mcParticle.y()) > cfgCutY) {
            return;
          }
          histos.fill(HIST("MC/ka/neg/prm/pt/num"), track.pt(), multiplicity, track.dcaXY());
          if (track.hasTOF()) {
            histos.fill(HIST("MC/ka/neg/prm/pt/numtof"), track.pt(), multiplicity, track.dcaXY());
          }
        }
        // Protons
        if (mcParticle.pdgCode() == 2212) {
          if (std::abs(mcParticle.y()) > cfgCutY) {
            return;
          }
          histos.fill(HIST("MC/pr/pos/prm/pt/num"), track.pt(), multiplicity, track.dcaXY());
          if (track.hasTOF()) {
            histos.fill(HIST("MC/pr/pos/prm/pt/numtof"), track.pt(), multiplicity, track.dcaXY());
          }
        }
        // Anti-Protons
        if (mcParticle.pdgCode() == -2212) {
          if (std::abs(mcParticle.y()) > cfgCutY) {
            return;
          }
          histos.fill(HIST("MC/pr/neg/prm/pt/num"), track.pt(), multiplicity, track.dcaXY());
          if (track.hasTOF()) {
            histos.fill(HIST("MC/pr/neg/prm/pt/numtof"), track.pt(), multiplicity, track.dcaXY());
          }
        }
      } // primaries
      if (!mcParticle.isPhysicalPrimary()) { // secondaries loop start
        if (mcParticle.pdgCode() == 211) {
          if (std::abs(mcParticle.y()) > cfgCutY) {
            return;
          }
          if (mcParticle.getProcess() == 4) {
            histos.fill(HIST("MC/pi/pos/str/pt/num"), track.pt(), multiplicity, track.dcaXY());
            if (track.hasTOF()) {
              histos.fill(HIST("MC/pi/pos/str/pt/numtof"), track.pt(), multiplicity, track.dcaXY());
            }
          } else {
            histos.fill(HIST("MC/pi/pos/mat/pt/num"), track.pt(), multiplicity, track.dcaXY());
            if (track.hasTOF()) {
              histos.fill(HIST("MC/pi/pos/mat/pt/numtof"), track.pt(), multiplicity, track.dcaXY());
            }
          }
        }
        if (mcParticle.pdgCode() == -211) {
          if (std::abs(mcParticle.y()) > cfgCutY) {
            return;
          }
          if (mcParticle.getProcess() == 4) {
            histos.fill(HIST("MC/pi/neg/str/pt/num"), track.pt(), multiplicity, track.dcaXY());
            if (track.hasTOF()) {
              histos.fill(HIST("MC/pi/neg/str/pt/numtof"), track.pt(), multiplicity, track.dcaXY());
            }
          } else {
            histos.fill(HIST("MC/pi/neg/mat/pt/num"), track.pt(), multiplicity, track.dcaXY());
            if (track.hasTOF()) {
              histos.fill(HIST("MC/pi/neg/mat/pt/numtof"), track.pt(), multiplicity, track.dcaXY());
            }
          }
        }
        if (mcParticle.pdgCode() == 321) {
          if (std::abs(mcParticle.y()) > cfgCutY) {
            return;
          }
          if (mcParticle.getProcess() == 4) {
            histos.fill(HIST("MC/ka/pos/str/pt/num"), track.pt(), multiplicity, track.dcaXY());
            if (track.hasTOF()) {
              histos.fill(HIST("MC/ka/pos/str/pt/numtof"), track.pt(), multiplicity, track.dcaXY());
            }
          } else {
            histos.fill(HIST("MC/ka/pos/mat/pt/num"), track.pt(), multiplicity, track.dcaXY());
            if (track.hasTOF()) {
              histos.fill(HIST("MC/ka/pos/mat/pt/numtof"), track.pt(), multiplicity, track.dcaXY());
            }
          }
        }
        if (mcParticle.pdgCode() == -321) {
          if (std::abs(mcParticle.y()) > cfgCutY) {
            return;
          }
          if (mcParticle.getProcess() == 4) {
            histos.fill(HIST("MC/ka/neg/str/pt/num"), track.pt(), multiplicity, track.dcaXY());
            if (track.hasTOF()) {
              histos.fill(HIST("MC/ka/neg/str/pt/numtof"), track.pt(), multiplicity, track.dcaXY());
            }
          } else {
            histos.fill(HIST("MC/ka/neg/mat/pt/num"), track.pt(), multiplicity, track.dcaXY());
            if (track.hasTOF()) {
              histos.fill(HIST("MC/ka/neg/mat/pt/numtof"), track.pt(), multiplicity, track.dcaXY());
            }
          }
        }
        if (mcParticle.pdgCode() == 2212) {
          if (std::abs(mcParticle.y()) > cfgCutY) {
            return;
          }
          if (mcParticle.getProcess() == 4) {
            histos.fill(HIST("MC/pr/pos/str/pt/num"), track.pt(), multiplicity, track.dcaXY());
            if (track.hasTOF()) {
              histos.fill(HIST("MC/pr/pos/str/pt/numtof"), track.pt(), multiplicity, track.dcaXY());
            }
          } else {
            histos.fill(HIST("MC/pr/pos/mat/pt/num"), track.pt(), multiplicity, track.dcaXY());
            if (track.hasTOF()) {
              histos.fill(HIST("MC/pr/pos/mat/pt/numtof"), track.pt(), multiplicity, track.dcaXY());
            }
          }
        }
        if (mcParticle.pdgCode() == -2212) {
          if (std::abs(mcParticle.y()) > cfgCutY) {
            return;
          }
          if (mcParticle.getProcess() == 4) {
            histos.fill(HIST("MC/pr/neg/str/pt/num"), track.pt(), multiplicity, track.dcaXY());
            if (track.hasTOF()) {
              histos.fill(HIST("MC/pr/neg/str/pt/numtof"), track.pt(), multiplicity, track.dcaXY());
            }
          } else {
            histos.fill(HIST("MC/pr/neg/mat/pt/num"), track.pt(), multiplicity, track.dcaXY());
            if (track.hasTOF()) {
              histos.fill(HIST("MC/pr/neg/mat/pt/numtof"), track.pt(), multiplicity, track.dcaXY());
            }
          }
        }
      } // secondaries loop end
    }
  }
  PROCESS_SWITCH(QCspectraTPC, processTrackHistograms_MC, "process MC Reconstructed tracks", true);

  void processParticleHistograms_MC(CollisionCandidateMCRec::iterator const& collision, soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::TrackSelection> const& tracks,
                                    aod::McParticles const& /*mcParticles*/,
                                    aod::McCollisions const& /*mcCollision*/)
  {

    const float multiplicity = collision.centFT0C();
    for (auto& track : tracks) {
      const auto& mcParticle = track.mcParticle();
      if (mcParticle.isPhysicalPrimary()) {
        // pions
        if (mcParticle.pdgCode() == 211) {
          if (std::abs(mcParticle.y()) > cfgCutY) {
            return;
          }
          histos.fill(HIST("MC/pi/pos/prm/pt/den"), mcParticle.pt(), multiplicity, track.dcaXY());
        }
        if (mcParticle.pdgCode() == -211) {
          if (std::abs(mcParticle.y()) > cfgCutY) {
            return;
          }
          histos.fill(HIST("MC/pi/neg/prm/pt/den"), mcParticle.pt(), multiplicity, track.dcaXY());
        }
        if (mcParticle.pdgCode() == 321) {
          if (std::abs(mcParticle.y()) > cfgCutY) {
            return;
          }
          histos.fill(HIST("MC/ka/pos/prm/pt/den"), mcParticle.pt(), multiplicity, track.dcaXY());
        }
        if (mcParticle.pdgCode() == -321) {
          if (std::abs(mcParticle.y()) > cfgCutY) {
            return;
          }
          histos.fill(HIST("MC/ka/neg/prm/pt/den"), mcParticle.pt(), multiplicity, track.dcaXY());
        }

        if (mcParticle.pdgCode() == 2212) {
          if (std::abs(mcParticle.y()) > cfgCutY) {
            return;
          }
          histos.fill(HIST("MC/pr/pos/prm/pt/den"), mcParticle.pt(), multiplicity, track.dcaXY());
        }
        if (mcParticle.pdgCode() == -2212) {
          if (std::abs(mcParticle.y()) > cfgCutY) {
            return;
          }
          histos.fill(HIST("MC/pr/neg/prm/pt/den"), mcParticle.pt(), multiplicity, track.dcaXY());
        }
        if (!mcParticle.isPhysicalPrimary()) { // secondaries loop start
          if (mcParticle.pdgCode() == 211) {
            if (std::abs(mcParticle.y()) > cfgCutY) {
              return;
            }
            if (mcParticle.getProcess() == 4) {
              histos.fill(HIST("MC/pi/pos/str/pt/den"), track.pt(), multiplicity, track.dcaXY());
            } else {
              histos.fill(HIST("MC/pi/pos/mat/pt/den"), track.pt(), multiplicity, track.dcaXY());
            }
          }
          if (mcParticle.pdgCode() == -211) {
            if (std::abs(mcParticle.y()) > cfgCutY) {
              return;
            }
            if (mcParticle.getProcess() == 4) {
              histos.fill(HIST("MC/pi/neg/str/pt/den"), track.pt(), multiplicity, track.dcaXY());
            } else {
              histos.fill(HIST("MC/pi/neg/mat/pt/den"), track.pt(), multiplicity, track.dcaXY());
            }
          }
          if (mcParticle.pdgCode() == 321) {
            if (std::abs(mcParticle.y()) > cfgCutY) {
              return;
            }
            if (mcParticle.getProcess() == 4) {
              histos.fill(HIST("MC/ka/pos/str/pt/den"), track.pt(), multiplicity, track.dcaXY());
            } else {
              histos.fill(HIST("MC/ka/pos/mat/pt/den"), track.pt(), multiplicity, track.dcaXY());
            }
          }
          if (mcParticle.pdgCode() == -321) {
            if (std::abs(mcParticle.y()) > cfgCutY) {
              return;
            }
            if (mcParticle.getProcess() == 4) {
              histos.fill(HIST("MC/ka/neg/str/pt/den"), track.pt(), multiplicity, track.dcaXY());
            } else {
              histos.fill(HIST("MC/ka/neg/mat/pt/den"), track.pt(), multiplicity, track.dcaXY());
            }
          }
          if (mcParticle.pdgCode() == 2212) {
            if (std::abs(mcParticle.y()) > cfgCutY) {
              return;
            }
            if (mcParticle.getProcess() == 4) {
              histos.fill(HIST("MC/pr/pos/str/pt/den"), track.pt(), multiplicity, track.dcaXY());
            } else {
              histos.fill(HIST("MC/pr/pos/mat/pt/den"), track.pt(), multiplicity, track.dcaXY());
            }
          }
          if (mcParticle.pdgCode() == -2212) {
            if (std::abs(mcParticle.y()) > cfgCutY) {
              return;
            }
            if (mcParticle.getProcess() == 4) {
              histos.fill(HIST("MC/pr/neg/str/pt/den"), track.pt(), multiplicity, track.dcaXY());
            } else {
              histos.fill(HIST("MC/pr/neg/mat/pt/den"), track.pt(), multiplicity, track.dcaXY());
            }
          }
        } // secondaries loop end
      }
    }
  } // process_mc loop end
  PROCESS_SWITCH(QCspectraTPC, processParticleHistograms_MC, "process MC Generated tracks", true);
  void processParticleGen(aod::McCollision const& /*mcCollision*/, aod::McParticles& mcParticles)
  {

    for (auto& mcParticleGen : mcParticles) {
      bool isPhysPrim = mcParticleGen.isPhysicalPrimary();
      if (doProcessMC) {
        if (isPhysPrim) {
          if (mcParticleGen.pdgCode() == 2212) {
            if (abs(mcParticleGen.y()) > std::abs(cfgCutY)) {
              continue;
            }
            histos.fill(HIST("MC/test/pr/pos/prm/pt/den"), mcParticleGen.pt());
          }
          if (mcParticleGen.pdgCode() == -2212) {
            if (abs(mcParticleGen.y()) > std::abs(cfgCutY)) {
              continue;
            }
            histos.fill(HIST("MC/test/pr/neg/prm/pt/den"), mcParticleGen.pt());
          }
          if (mcParticleGen.pdgCode() == 211) {
            if (abs(mcParticleGen.y()) > std::abs(cfgCutY)) {
              continue;
            }
            histos.fill(HIST("MC/test/pi/pos/prm/pt/den"), mcParticleGen.pt());
          }
          if (mcParticleGen.pdgCode() == -211) {
            if (abs(mcParticleGen.y()) > std::abs(cfgCutY)) {
              continue;
            }
            histos.fill(HIST("MC/test/pi/neg/prm/pt/den"), mcParticleGen.pt());
          }
          if (mcParticleGen.pdgCode() == 321) {
            if (abs(mcParticleGen.y()) > std::abs(cfgCutY)) {
              continue;
            }
            histos.fill(HIST("MC/test/ka/pos/prm/pt/den"), mcParticleGen.pt());
          }
          if (mcParticleGen.pdgCode() == -321) {
            if (abs(mcParticleGen.y()) > std::abs(cfgCutY)) {
              continue;
            }
            histos.fill(HIST("MC/test/ka/neg/prm/pt/den"), mcParticleGen.pt());
          }
        }
      }
    }
  } // process_mc loop end
  PROCESS_SWITCH(QCspectraTPC, processParticleGen, "process MC Generated tracks", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<QCspectraTPC>(cfgc)};
}
