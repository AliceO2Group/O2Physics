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

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "TVector2.h"

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

struct deutRtTask {

  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutLeadingPtMin{"cfgCutLeadingPtMin", 0.15f, "Accepted min range for leading pt"};
  Configurable<float> cfgAverage_Nch_Transv{"cfgAverage_Nch_Transv", 6.81f, "Average charge multiplicity for pp"};
  Configurable<float> cfgCutDCAXY{"cfgCutDCAXY", 0.05f, "DCAXY range for tracks"};
  Configurable<float> cfgCutDCAZ{"cfgCutDCAZ", 2.0f, "DCAZ range for tracks"};
  Configurable<float> cfgCutY{"cfgCutY", 0.5f, "Y range for tracks"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<float> maxChi2PerClusterTPC{"maxChi2PerClusterTPC", 4.f, "Additional cut on the maximum value of the chi2 per cluster in the TPC"};
  Configurable<float> maxChi2PerClusterITS{"maxChi2PerClusterITS", 36.f, "Additional cut on the maximum value of the chi2 per cluster in the ITS"};
  Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 70.f, "Additional cut on the minimum number of crossed rows in the TPC"};
  Configurable<float> minNCrossedRowsOverFindableClustersTPC{"minNCrossedRowsOverFindableClustersTPC", 0.8f, "Additional cut on the min value of ratio between crossed rows & findable clusters in TPC"};
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2, 9.4, 9.6, 9.8, 10.0}, "Binning of the pT axis"};
  ConfigurableAxis binsDca{"binsDca", {VARIABLE_WIDTH, -3.0, -2.95, -2.9, -2.85, -2.8, -2.75, -2.7, -2.65, -2.6, -2.55, -2.5, -2.45, -2.4, -2.35, -2.3, -2.25, -2.2, -2.15, -2.1, -2.05, -2.0, -1.975, -1.95, -1.925, -1.9, -1.875, -1.85, -1.825, -1.8, -1.775, -1.75, -1.725, -1.7, -1.675, -1.65, -1.625, -1.6, -1.575, -1.55, -1.525, -1.5, -1.475, -1.45, -1.425, -1.4, -1.375, -1.35, -1.325, -1.3, -1.275, -1.25, -1.225, -1.2, -1.175, -1.15, -1.125, -1.1, -1.075, -1.05, -1.025, -1.0, -0.99, -0.98, -0.97, -0.96, -0.95, -0.94, -0.93, -0.92, -0.91, -0.9, -0.89, -0.88, -0.87, -0.86, -0.85, -0.84, -0.83, -0.82, -0.81, -0.8, -0.79, -0.78, -0.77, -0.76, -0.75, -0.74, -0.73, -0.72, -0.71, -0.7, -0.69, -0.68, -0.67, -0.66, -0.65, -0.64, -0.63, -0.62, -0.61, -0.6, -0.59, -0.58, -0.57, -0.56, -0.55, -0.54, -0.53, -0.52, -0.51, -0.5, -0.49, -0.48, -0.47, -0.46, -0.45, -0.44, -0.43, -0.42, -0.41, -0.4, -0.396, -0.392, -0.388, -0.384, -0.38, -0.376, -0.372, -0.368, -0.364, -0.36, -0.356, -0.352, -0.348, -0.344, -0.34, -0.336, -0.332, -0.328, -0.324, -0.32, -0.316, -0.312, -0.308, -0.304, -0.3, -0.296, -0.292, -0.288, -0.284, -0.28, -0.276, -0.272, -0.268, -0.264, -0.26, -0.256, -0.252, -0.248, -0.244, -0.24, -0.236, -0.232, -0.228, -0.224, -0.22, -0.216, -0.212, -0.208, -0.204, -0.2, -0.198, -0.196, -0.194, -0.192, -0.19, -0.188, -0.186, -0.184, -0.182, -0.18, -0.178, -0.176, -0.174, -0.172, -0.17, -0.168, -0.166, -0.164, -0.162, -0.16, -0.158, -0.156, -0.154, -0.152, -0.15, -0.148, -0.146, -0.144, -0.142, -0.14, -0.138, -0.136, -0.134, -0.132, -0.13, -0.128, -0.126, -0.124, -0.122, -0.12, -0.118, -0.116, -0.114, -0.112, -0.11, -0.108, -0.106, -0.104, -0.102, -0.1, -0.099, -0.098, -0.097, -0.096, -0.095, -0.094, -0.093, -0.092, -0.091, -0.09, -0.089, -0.088, -0.087, -0.086, -0.085, -0.084, -0.083, -0.082, -0.081, -0.08, -0.079, -0.078, -0.077, -0.076, -0.075, -0.074, -0.073, -0.072, -0.071, -0.07, -0.069, -0.068, -0.067, -0.066, -0.065, -0.064, -0.063, -0.062, -0.061, -0.06, -0.059, -0.058, -0.057, -0.056, -0.055, -0.054, -0.053, -0.052, -0.051, -0.05, -0.049, -0.048, -0.047, -0.046, -0.045, -0.044, -0.043, -0.042, -0.041, -0.04, -0.039, -0.038, -0.037, -0.036, -0.035, -0.034, -0.033, -0.032, -0.031, -0.03, -0.029, -0.028, -0.027, -0.026, -0.025, -0.024, -0.023, -0.022, -0.021, -0.02, -0.019, -0.018, -0.017, -0.016, -0.015, -0.014, -0.013, -0.012, -0.011, -0.01, -0.009, -0.008, -0.007, -0.006, -0.005, -0.004, -0.003, -0.002, -0.001, -0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.017, 0.018, 0.019, 0.02, 0.021, 0.022, 0.023, 0.024, 0.025, 0.026, 0.027, 0.028, 0.029, 0.03, 0.031, 0.032, 0.033, 0.034, 0.035, 0.036, 0.037, 0.038, 0.039, 0.04, 0.041, 0.042, 0.043, 0.044, 0.045, 0.046, 0.047, 0.048, 0.049, 0.05, 0.051, 0.052, 0.053, 0.054, 0.055, 0.056, 0.057, 0.058, 0.059, 0.06, 0.061, 0.062, 0.063, 0.064, 0.065, 0.066, 0.067, 0.068, 0.069, 0.07, 0.071, 0.072, 0.073, 0.074, 0.075, 0.076, 0.077, 0.078, 0.079, 0.08, 0.081, 0.082, 0.083, 0.084, 0.085, 0.086, 0.087, 0.088, 0.089, 0.09, 0.091, 0.092, 0.093, 0.094, 0.095, 0.096, 0.097, 0.098, 0.099, 0.1, 0.102, 0.104, 0.106, 0.108, 0.11, 0.112, 0.114, 0.116, 0.118, 0.12, 0.122, 0.124, 0.126, 0.128, 0.13, 0.132, 0.134, 0.136, 0.138, 0.14, 0.142, 0.144, 0.146, 0.148, 0.15, 0.152, 0.154, 0.156, 0.158, 0.16, 0.162, 0.164, 0.166, 0.168, 0.17, 0.172, 0.174, 0.176, 0.178, 0.18, 0.182, 0.184, 0.186, 0.188, 0.19, 0.192, 0.194, 0.196, 0.198, 0.2, 0.204, 0.208, 0.212, 0.216, 0.22, 0.224, 0.228, 0.232, 0.236, 0.24, 0.244, 0.248, 0.252, 0.256, 0.26, 0.264, 0.268, 0.272, 0.276, 0.28, 0.284, 0.288, 0.292, 0.296, 0.3, 0.304, 0.308, 0.312, 0.316, 0.32, 0.324, 0.328, 0.332, 0.336, 0.34, 0.344, 0.348, 0.352, 0.356, 0.36, 0.364, 0.368, 0.372, 0.376, 0.38, 0.384, 0.388, 0.392, 0.396, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.0, 1.025, 1.05, 1.075, 1.1, 1.125, 1.15, 1.175, 1.2, 1.225, 1.25, 1.275, 1.3, 1.325, 1.35, 1.375, 1.4, 1.425, 1.45, 1.475, 1.5, 1.525, 1.55, 1.575, 1.6, 1.625, 1.65, 1.675, 1.7, 1.725, 1.75, 1.775, 1.8, 1.825, 1.85, 1.875, 1.9, 1.925, 1.95, 1.975, 2.0, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.85, 2.9, 2.95, 3.0}, "Binning of DCA xy and z axis"};
  Configurable<bool> fIspPb{"fIspPb", false, "Type of collision"};

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  void init(o2::framework::InitContext&)
  {
    const AxisSpec axisEta{120, -1.5, +1.5, "#eta"};
    const AxisSpec ptAxis{binsPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec dcaXyAxis{binsDca, "DCA_{xy} (cm)"};
    const AxisSpec dcaZAxis{binsDca, "DCA_{z} (cm)"};
    const AxisSpec axisY{120, -1.5, +1.5, "y"};

    histos.add("event/vertexz", "collision z position", HistType::kTH1F, {{100, -10., +10., "z position (cm)"}});
    histos.add("etaHistogram", "etaHistogram", kTH1F, {axisEta});
    //  histos.add("rapidityHistogram", "rapidityHistogram", kTH1F, {axisY});
    histos.add("ptHistogram", "ptHistogram", kTH1F, {{200, 0, 10, "#it{p}_{T} (GeV/#it{c})"}});
    histos.add("leadingptHistogram", "leadingptHistogram", kTH1F, {{200, 0, 10, "#it{p}_{T} (GeV/#it{c})"}});
    histos.add("hMultTransverse", "hMultTransverse", kTH1I, {{200, 0, 200, ""}});
    histos.add("hMultToward", "hMultToward", kTH1I, {{200, 0, 200, ""}});
    histos.add("hMultAway", "hMultAway", kTH1I, {{200, 0, 200, ""}});
    histos.add("hRtDistribution", "hRtDistribution", kTH1F, {{200, 0, 20, ""}});

    // Arrays
    double Nch_Tr[] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 5.0, 10.0};
    double pt_TPC[] = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5};
    double pt_TOF[] = {1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 5.0};
    double pt_DCA[] = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0};
    double dca_xy[] = {-1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.18, -0.16, -0.14, -0.12, -0.10, -0.08, -0.06, -0.04, -0.02, 0.0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

    // Arrays Dimensions
    int nBins_Nch_Tr = sizeof(Nch_Tr) / sizeof(Double_t) - 1;
    int nBins_pt_TPC = sizeof(pt_TPC) / sizeof(Double_t) - 1;
    int nBins_pt_TOF = sizeof(pt_TOF) / sizeof(Double_t) - 1;
    int nBins_pt_DCA = sizeof(pt_DCA) / sizeof(Double_t) - 1;
    int nBins_dca_xy = sizeof(dca_xy) / sizeof(Double_t) - 1;

    //***************************  Nch_Transverse,      DCA_{xy},                 pt_DCA ***************************************************

    const AxisSpec Nch_TrAxis{nBins_Nch_Tr, Nch_Tr[0], Nch_Tr[nBins_Nch_Tr], "Nch_Tr"};
    const AxisSpec dcaXyAxis1{nBins_dca_xy, -1.0, +1.0, "DCA_{xy} (cm)"};
    const AxisSpec pTdcaZAxis{nBins_pt_DCA, pt_DCA[0], pt_DCA[nBins_pt_DCA], "#it{p} (GeV/#it{c})"};

    // histos.add("hDCAxy_deuterons_Toward", "", kTHnSparseD, {Nch_TrAxis, dcaXyAxis1, pTdcaZAxis});
    // histos.add("hDCAxy_antideuterons_Toward", "", kTHnSparseD, {Nch_TrAxis, dcaXyAxis1, pTdcaZAxis});

    histos.add("deuterons/hDCAxy_deuterons_Toward", "", kTHnSparseD, {Nch_TrAxis, dcaXyAxis1, pTdcaZAxis});
    histos.add("anti-deuterons/hDCAxy_antideuterons_Toward", "", kTHnSparseD, {Nch_TrAxis, dcaXyAxis1, pTdcaZAxis});
    histos.add("deuterons/hDCAxy_deuterons_Away", "", kTHnSparseD, {Nch_TrAxis, dcaXyAxis1, pTdcaZAxis});
    histos.add("anti-deuterons/hDCAxy_antideuterons_Away", "", kTHnSparseD, {Nch_TrAxis, dcaXyAxis1, pTdcaZAxis});
    histos.add("deuterons/hDCAxy_deuterons_Transverse", "", kTHnSparseD, {Nch_TrAxis, dcaXyAxis1, pTdcaZAxis});
    histos.add("anti-deuterons/hDCAxy_antideuterons_Transverse", "", kTHnSparseD, {Nch_TrAxis, dcaXyAxis1, pTdcaZAxis});

    //***************************  Nch_Transverse,  nsigma_{TPC},                 pt_TPC ***************************************************
    const AxisSpec nsigmaTPCAxis{100, -10.0, +10.0, "N_{#sigma}^{TPC}"};
    const AxisSpec pT_TPCAxis{nBins_pt_TPC, pt_TPC[0], pt_TPC[nBins_pt_TPC], "#it{p} (GeV/#it{c})"};

    // histos.add("hnsigmaTPC_deuterons_Toward", "", kTHnSparseD, {Nch_TrAxis, nsigmaTPCAxis, pT_TPCAxis});
    // histos.add("hnsigmaTPC_antideuterons_Toward", "", kTHnSparseD, {Nch_TrAxis, nsigmaTPCAxis, pT_TPCAxis});

    histos.add("deuterons/hnsigmaTPC_deuterons_Toward", "", kTHnSparseD, {Nch_TrAxis, nsigmaTPCAxis, pT_TPCAxis});
    histos.add("anti-deuterons/hnsigmaTPC_antideuterons_Toward", "", kTHnSparseD, {Nch_TrAxis, nsigmaTPCAxis, pT_TPCAxis});
    histos.add("deuterons/hnsigmaTPC_deuterons_Away", "", kTHnSparseD, {Nch_TrAxis, nsigmaTPCAxis, pT_TPCAxis});
    histos.add("anti-deuterons/hnsigmaTPC_antideuterons_Away", "", kTHnSparseD, {Nch_TrAxis, nsigmaTPCAxis, pT_TPCAxis});
    histos.add("deuterons/hnsigmaTPC_deuterons_Transverse", "", kTHnSparseD, {Nch_TrAxis, nsigmaTPCAxis, pT_TPCAxis});
    histos.add("anti-deuterons/hnsigmaTPC_antideuterons_Transverse", "", kTHnSparseD, {Nch_TrAxis, nsigmaTPCAxis, pT_TPCAxis});
    //***************************  Nch_Transverse,  nsigma_{TOF},                 pt_TOF ***************************************************

    const AxisSpec nsigmaTOFAxis{200, -10.0, +10.0, "N_{#sigma}^{TOF}"};
    const AxisSpec pT_TOFAxis{nBins_pt_TOF, pt_TOF[0], pt_TOF[nBins_pt_TOF], "#it{p} (GeV/#it{c})"};

    // histos.add("hnsigmaTOF_deuterons_Toward", "", kTHnSparseD, {Nch_TrAxis, nsigmaTOFAxis, pT_TOFAxis});
    // histos.add("hnsigmaTOF_antideuterons_Toward", "", kTHnSparseD, {Nch_TrAxis, nsigmaTOFAxis, pT_TOFAxis});
    histos.add("deuterons/hnsigmaTOF_deuterons_Toward", "", kTHnSparseD, {Nch_TrAxis, nsigmaTOFAxis, pT_TOFAxis});
    histos.add("anti-deuterons/hnsigmaTOF_antideuterons_Toward", "", kTHnSparseD, {Nch_TrAxis, nsigmaTOFAxis, pT_TOFAxis});
    histos.add("deuterons/hnsigmaTOF_deuterons_Away", "", kTHnSparseD, {Nch_TrAxis, nsigmaTOFAxis, pT_TOFAxis});
    histos.add("anti-deuterons/hnsigmaTOF_antideuterons_Away", "", kTHnSparseD, {Nch_TrAxis, nsigmaTOFAxis, pT_TOFAxis});
    histos.add("deuterons/hnsigmaTOF_deuterons_Transverse", "", kTHnSparseD, {Nch_TrAxis, nsigmaTOFAxis, pT_TOFAxis});
    histos.add("anti-deuterons/hnsigmaTOF_antideuterons_Transverse", "", kTHnSparseD, {Nch_TrAxis, nsigmaTOFAxis, pT_TOFAxis});
    histos.add("histTpcSignal_Deut", "Specific energy loss of Deuteron", HistType::kTH2F, {{600, -6., 6., "#it{p/z} (GeV/#it{c})"}, {1400, 0, 1400, "d#it{E} / d#it{X} (a. u.)"}});
    histos.add("histTpcSignal_CleanDeut", "Specific energy loss of Clean Deuteron", HistType::kTH2F, {{600, -6., 6., "#it{p/z} (GeV/#it{c})"}, {1400, 0, 1400, "d#it{E} / d#it{X} (a. u.)"}});
  }

  //***********************************************************************************************************************************
  template <typename T>
  bool PassedTrackQualityCuts_LeadingTrack(const T& track)
  {
    // Initialization
    bool isLeadingTrack = kFALSE;

    // Check additional condition on TPC chi2/N clusters
    if (track.tpcChi2NCl() >= maxChi2PerClusterTPC) {
      return isLeadingTrack; // If this condition is not met, it's not a leading track
    }

    if (track.itsChi2NCl() >= maxChi2PerClusterITS) {
      return isLeadingTrack; // If this condition is not met, it's not a leading track
    }

    if (track.tpcNClsCrossedRows() < minNCrossedRowsTPC) {
      return isLeadingTrack; // If this condition is not met, it's not a leading track
    }

    if (track.tpcCrossedRowsOverFindableCls() < minNCrossedRowsOverFindableClustersTPC) {
      return isLeadingTrack; // If this condition is not met, it's not a leading track
    }
    isLeadingTrack = kTRUE; // If all conditions are met, it's a leading track
    return isLeadingTrack;
  }

  //***********************************************************************************************************************************
  template <typename T2>
  bool PassedBasicTrackQualityCuts(const T2& track)
  {

    // Initialization
    bool passedTrkSelection = kFALSE;
    if (TMath::Abs(track.eta()) > 0.8) {
      return passedTrkSelection;
    }
    // Rapidity Cut
    Double_t mass = 1.875;
    Double_t p = track.p();
    Double_t pz = track.pz();
    Double_t E = TMath::Sqrt(mass * mass + p * p);
    if (E == TMath::Abs(pz)) {
      return passedTrkSelection;
    }
    Double_t y_lab = 0.5 * TMath::Log((E + pz) / (E - pz));
    Double_t y_cms(0);

    // Rapidity Selection in p-Pb Collisions
    if (fIspPb) {
      y_cms = y_lab - 0.465;
      if (y_cms < -1.0) {
        return passedTrkSelection;
      }
      if (y_cms > 0.0) {
        return passedTrkSelection;
      }
    }

    // Rapidity Selection in pp Collisions
    if (!fIspPb) {
      y_cms = y_lab;
      if (TMath::Abs(y_cms) > 0.5) {
        return passedTrkSelection;
      }
    }

    passedTrkSelection = kTRUE;
    return passedTrkSelection;
  }
  //***********************************************************************************************************************************
  template <typename T>
  bool IsDeuteronCandidate(const T& track)
  {
    // Initialization
    bool isDeuteronCandidate = kFALSE;

    if (!PassedBasicTrackQualityCuts(track)) {
      return isDeuteronCandidate;
    }

    double pt = track.pt();
    bool hasTOFhit = track.hasTOF();
    // double length    = track.GetIntegratedLength();
    double length = track.length();
    const auto& nsigmaTOF = o2::aod::pidutils::tofNSigma<PID::Deuteron>(track);
    //   const auto& nsigmaTPC = o2::aod::pidutils::tpcNSigma<PID::Deuteron>(track);

    if (pt > 1.5 && (!hasTOFhit)) {
      return isDeuteronCandidate;
    }
    if (pt > 1.5 && length < 350.0) {
      return isDeuteronCandidate;
    }
    if (pt > 1.5 && TMath::Abs(nsigmaTOF) > 10.0) {
      return isDeuteronCandidate;
    }
    isDeuteronCandidate = kTRUE;
    return isDeuteronCandidate;
  }
  //***********************************************************************************************************************************
  template <typename T5>
  bool IsCleanDeuteron(const T5& track)
  {
    // Initialization
    bool isCleanDeuteron = kFALSE;
    // Variables
    double nsigmaTPC = o2::aod::pidutils::tpcNSigma<PID::Deuteron>(track);
    double nsigmaTOF = o2::aod::pidutils::tofNSigma<PID::Deuteron>(track);
    double pt = track.pt();
    double hasTOFhit = track.hasTOF();
    double length = track.length();

    // Selections
    if (pt < 1.0 && abs(nsigmaTPC) < 3.0)
      isCleanDeuteron = kTRUE;
    if (pt > 1.0 && hasTOFhit && length > 350.0 && abs(nsigmaTPC) < 3.0 && abs(nsigmaTOF) < 3.0)
      isCleanDeuteron = kTRUE;

    return isCleanDeuteron;
  }
  //***********************************************************************************************************************************
  template <typename T6>
  void FillHistograms_StandardCuts(int mult_Transverse, double delta_phi, T6 const& track)
  {
    // Variables
    double nsigmaTPC = o2::aod::pidutils::tpcNSigma<PID::Deuteron>(track);
    double nsigmaTOF = o2::aod::pidutils::tofNSigma<PID::Deuteron>(track);
    double pt = track.pt();
    double charge = track.sign();
    double DCAxy = track.dcaXY();
    double DCAz = track.dcaZ();
    double hasTOFhit = track.hasTOF();
    double length = track.length();
    double Rt = static_cast<double>(mult_Transverse) / cfgAverage_Nch_Transv;
    // DCA_{z} Cut
    if (abs(DCAz) > 1.0) {
      return;
    }
    // DCA_{xy} Histograms
    if (IsCleanDeuteron(track)) {

      if (charge > 0) {
        if ((delta_phi >= 0.0 && delta_phi < 60.0) || (delta_phi >= 300.0 && delta_phi <= 360.)) {
          histos.fill(HIST("deuterons/hDCAxy_deuterons_Toward"), Rt, DCAxy, pt);
        }
        if (delta_phi >= 120.0 && delta_phi < 240.0) {
          histos.fill(HIST("deuterons/hDCAxy_deuterons_Away"), Rt, DCAxy, pt);
        }
        if ((delta_phi >= 60.0 && delta_phi < 120.0) || (delta_phi >= 240.0 && delta_phi < 300.0)) {
          histos.fill(HIST("deuterons/hDCAxy_deuterons_Transverse"), Rt, DCAxy, pt);
        }
      }

      if (charge < 0) {
        if ((delta_phi >= 0.0 && delta_phi < 60.0) || (delta_phi >= 300.0 && delta_phi <= 360.)) {
          histos.fill(HIST("anti-deuterons/hDCAxy_antideuterons_Toward"), Rt, DCAxy, pt);
        }
        if (delta_phi >= 120.0 && delta_phi < 240.0) {
          histos.fill(HIST("anti-deuterons/hDCAxy_antideuterons_Away"), Rt, DCAxy, pt);
        }
        if ((delta_phi >= 60.0 && delta_phi < 120.0) || (delta_phi >= 240.0 && delta_phi < 300.0)) {
          histos.fill(HIST("anti-deuterons/hDCAxy_antideuterons_Transverse"), Rt, DCAxy, pt);
        }
      }
    } // end of clean deuteron tracks

    // DCA_{xy} Cut
    if (TMath::Abs(DCAxy) > 0.1) {
      return;
    }
    if (abs(nsigmaTPC) > 5.0) {
      return;
    }
    // TPC-Only Analysis
    if (pt < 1.5) {

      if (charge > 0) {
        if ((delta_phi >= 0.0 && delta_phi < 60.0) || (delta_phi >= 300.0 && delta_phi <= 360.)) {
          histos.fill(HIST("deuterons/hnsigmaTPC_deuterons_Toward"), Rt, nsigmaTPC, pt);
        }
        if (delta_phi >= 120.0 && delta_phi < 240.0) {
          histos.fill(HIST("deuterons/hnsigmaTPC_deuterons_Away"), Rt, nsigmaTPC, pt);
        }
        if ((delta_phi >= 60.0 && delta_phi < 120.0) || (delta_phi >= 240.0 && delta_phi < 300.0)) {
          histos.fill(HIST("deuterons/hnsigmaTPC_deuterons_Transverse"), Rt, nsigmaTPC, pt);
        }
      }

      if (charge < 0) {
        if ((delta_phi >= 0.0 && delta_phi < 60.0) || (delta_phi >= 300.0 && delta_phi <= 360.)) {
          histos.fill(HIST("anti-deuterons/hnsigmaTPC_antideuterons_Toward"), Rt, nsigmaTPC, pt);
        }
        if (delta_phi >= 120.0 && delta_phi < 240.0) {
          histos.fill(HIST("anti-deuterons/hnsigmaTPC_antideuterons_Away"), Rt, nsigmaTPC, pt);
        }
        if ((delta_phi >= 60.0 && delta_phi < 120.0) || (delta_phi >= 240.0 && delta_phi < 300.0)) {
          histos.fill(HIST("anti-deuterons/hnsigmaTPC_antideuterons_Transverse"), Rt, nsigmaTPC, pt);
        }
      }
    }
    // TOF Analysis
    if (!hasTOFhit) {
      return;
    }
    if (length < 350.0) {
      return;
    }

    if (pt < 1.0) {
      return;
    }

    if (charge > 0) {
      if ((delta_phi >= 0.0 && delta_phi < 60.0) || (delta_phi >= 300.0 && delta_phi <= 360.)) {
        histos.fill(HIST("deuterons/hnsigmaTOF_deuterons_Toward"), Rt, nsigmaTOF, pt);
      }
      if (delta_phi >= 120.0 && delta_phi < 240.0) {
        histos.fill(HIST("deuterons/hnsigmaTOF_deuterons_Away"), Rt, nsigmaTOF, pt);
      }
      if ((delta_phi >= 60.0 && delta_phi < 120.0) || (delta_phi >= 240.0 && delta_phi < 300.0)) {
        histos.fill(HIST("deuterons/hnsigmaTOF_deuterons_Transverse"), Rt, nsigmaTOF, pt);
      }
    }

    if (charge < 0) {
      if ((delta_phi >= 0.0 && delta_phi < 60.0) || (delta_phi >= 300.0 && delta_phi <= 360.)) {
        histos.fill(HIST("anti-deuterons/hnsigmaTOF_antideuterons_Toward"), Rt, nsigmaTOF, pt);
      }
      if (delta_phi >= 120.0 && delta_phi < 240.0) {
        histos.fill(HIST("anti-deuterons/hnsigmaTOF_antideuterons_Away"), Rt, nsigmaTOF, pt);
      }
      if ((delta_phi >= 60.0 && delta_phi < 120.0) || (delta_phi >= 240.0 && delta_phi < 300.0)) {
        histos.fill(HIST("anti-deuterons/hnsigmaTOF_antideuterons_Transverse"), Rt, nsigmaTOF, pt);
      }
    }
  }
  //***********************************************************************************************************************************
  using CollisionCandidate = soa::Join<aod::Collisions, aod::EvSels /*, aod::CentFT0Cs*/>;
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullDe, aod::pidTOFFullDe>;

  void process(CollisionCandidate::iterator const& collision, TrackCandidates const& tracks)
  {

    if (!collision.sel8()) {
      return;
    }
    if (abs(collision.posZ()) > cfgCutVertex) {
      return;
    }
    histos.fill(HIST("event/vertexz"), collision.posZ());
    double leadingTrackpT = 0.0;
    double leadingTrackPhi = 0.0;
    for (auto& track : tracks) {

      if (abs(track.eta()) > cfgCutEta) {
        return;
      }

      histos.fill(HIST("ptHistogram"), track.pt());
      histos.fill(HIST("etaHistogram"), track.eta());

      if (PassedTrackQualityCuts_LeadingTrack(track)) {
        if (track.pt() > cfgCutLeadingPtMin && track.pt() < 10.0) {
          if (track.pt() > leadingTrackpT) {
            leadingTrackpT = track.pt();
            leadingTrackPhi = track.phi();
          }
        }
      }
      if (leadingTrackpT == 0.0) {
        return;
      }

      histos.fill(HIST("leadingptHistogram"), leadingTrackpT);
    }

    int mult_Transverse(0);
    int mult_Toward(0);
    int mult_Away(0);

    for (auto& track : tracks) {

      double phi_ref = TVector2::Phi_0_2pi(leadingTrackPhi);
      double phi_trk = TVector2::Phi_0_2pi(track.phi());
      double delta_phi = (180.0 / TMath::Pi()) * TVector2::Phi_0_2pi(phi_trk - phi_ref);

      if ((delta_phi >= 60.0 && delta_phi < 120.0) || (delta_phi >= 240.0 && delta_phi < 300.0)) {
        mult_Transverse++;
      }
      if ((delta_phi >= 0.0 && delta_phi < 60.0) || (delta_phi >= 300.0 && delta_phi <= 360.)) {
        mult_Toward++;
      }
      if (delta_phi >= 120.0 && delta_phi < 240.0) {
        mult_Away++;
      }
      histos.fill(HIST("hMultTransverse"), mult_Transverse);
      histos.fill(HIST("hMultToward"), mult_Toward);
      histos.fill(HIST("hMultAway"), mult_Away);

      double Rt = static_cast<double>(mult_Transverse) / cfgAverage_Nch_Transv;
      histos.fill(HIST("hRtDistribution"), Rt);
      if (IsCleanDeuteron(track)) {
        histos.fill(HIST("histTpcSignal_CleanDeut"), track.tpcInnerParam() * track.sign(), track.tpcSignal());
      }
      if (IsDeuteronCandidate(track)) {
        histos.fill(HIST("histTpcSignal_Deut"), track.tpcInnerParam() * track.sign(), track.tpcSignal());
        FillHistograms_StandardCuts(mult_Transverse, delta_phi, track);
      }
    } // end of tracks loop

  } // end of process function
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<deutRtTask>(cfgc)};
}
