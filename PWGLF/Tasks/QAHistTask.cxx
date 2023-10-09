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
//
// Authors: Rafael Manhart,
// Date: 05.10.2023

#include <cmath>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>

#include "ReconstructionDataFormats/Track.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"

#include "Framework/HistogramRegistry.h"
#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "TPDGCode.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct QAHistTask {

  // Data
  HistogramRegistry QA_reg{"all_species", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry QA_p{"proton", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry QA_antip{"antiproton", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry QA_d{"deuteron", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry QA_antid{"antideuteron", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry QA_t{"triton", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry QA_antit{"antitriton", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry QA_He3{"Heelium-3", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry QA_antiHe3{"antiHelium-3", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry QA_He4{"Helium-4", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry QA_antiHe4{"antiHelium-4", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  void init(o2::framework::InitContext&)
  {

    std::vector<double> ptBinning = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 3.6, 4., 5., 6., 8., 10., 12., 14.};
    std::vector<double> centBinning = {0., 1., 5., 10., 20., 30., 40., 50., 70., 100.};
    std::vector<double> etaBinning = {-1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec centAxis = {centBinning, "V0M (%)"};
    AxisSpec centralityAxis = {100, 0.0, 100.0, "VT0C (%)"};
    AxisSpec centralityAxis_extended = {105, 0.0, 105.0, "VT0C (%)"};
    AxisSpec etaAxis = {etaBinning, "#eta"};

    // +++++++++++++++++++++ Data ++++++++++++++++++++++++

    // QA all species
    QA_reg.add("histRecVtxZData", "collision z position", HistType::kTH1F, {{200, -20., +20., "z position (cm)"}});
    QA_reg.add("histCentrality", "Centrality", HistType::kTH1F, {centralityAxis_extended});
    QA_reg.add("histTpcSignalData", "Specific energy loss", HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}});
    QA_reg.add("histTofSignalData", "TOF signal", HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    QA_reg.add("histDcaVsPtData_particle", "dcaXY vs Pt (particle)", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    QA_reg.add("histDcaZVsPtData_particle", "dcaZ vs Pt (particle)", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    QA_reg.add("histDcaVsPtData_antiparticle", "dcaXY vs Pt (antiparticle)", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    QA_reg.add("histDcaZVsPtData_antiparticle", "dcaZ vs Pt (antiparticle)", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    QA_reg.add("histTOFm2", "TOF m^2 vs Pt", HistType::kTH2F, {ptAxis, {400, 0.0, 10.0, "m^2"}});
    QA_reg.add("histNClusterTPC", "Number of Clusters in TPC vs Pt", HistType::kTH2F, {ptAxis, {80, 0.0, 160.0, "nCluster"}});
    QA_reg.add("histNClusterITS", "Number of Clusters in ITS vs Pt", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    QA_reg.add("histNClusterITSib", "Number of Clusters in ib of ITS vs Pt", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    QA_reg.add("histTPCnClsFindable", "Findable TPC clusters", HistType::kTH2F, {ptAxis, {200, 0.0, 200.0, "nCluster"}});
    QA_reg.add("histTPCnClsFindableMinusFound", "TPC Clusters: Findable - Found", HistType::kTH2F, {ptAxis, {60, 0.0, 60.0, "nCluster"}});
    QA_reg.add("histTPCnClsFindableMinusCrossedRows", "TPC Clusters: Findable - crossed rows", HistType::kTH2F, {ptAxis, {60, 0.0, 60.0, "nCluster"}});
    QA_reg.add("histTPCnClsShared", "Number of shared TPC clusters", HistType::kTH2F, {ptAxis, {70, 0.0, 70.0, "nCluster"}});
    QA_reg.add("histTPCCrossedRowsOverFindableCls", "Ratio crossed rows over findable clusters", HistType::kTH1F, {{100, 0., 2.0, "Crossed Rows / Findable Cls"}});
    QA_reg.add("histTPCFoundOverFindable", "Ratio of found over findable clusters", HistType::kTH1F, {{100, 0., 2.0, "Found Cls / Findable Cls"}});
    QA_reg.add("histTPCFractionSharedCls", "Fraction of shared TPC clusters", HistType::kTH1F, {{100, -2.0, 2.0, "Shared Cls"}});
    QA_reg.add("histChi2TPC", "chi^2 TPC vs Pt", HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    QA_reg.add("histChi2ITS", "chi^2 ITS vs Pt", HistType::kTH2F, {ptAxis, {500, 0.0, 50.0, "chi^2"}});
    QA_reg.add("histChi2TOF", "chi^2 TOF vs Pt", HistType::kTH2F, {ptAxis, {75, 0.0, 15.0, "chi^2"}});
    QA_reg.add("histEtaWithOverFlow", "Pseudorapidity 0 - 105%% centrality", HistType::kTH1F, {etaAxis});
    QA_reg.add("histEta", "Pseudorapidity with centrality cut", HistType::kTH1F, {etaAxis});
    QA_reg.add("histEta_cent", "Pseudorapidity vs Centrality", HistType::kTH2F, {centralityAxis_extended, etaAxis});
    QA_reg.add("histTrackLength", "Track length", HistType::kTH1F, {{350, 0., 700., "length (cm)"}});


    // QA proton
    QA_p.add("histTpcSignalData", "Specific energy loss (p)", HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}});
    QA_p.add("histTofSignalData", "TOF signal (p)", HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    QA_p.add("histDcaVsPtData", "dcaXY vs Pt (p)", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    QA_p.add("histDcaZVsPtData", "dcaZ vs Pt (p)", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    QA_p.add("histTOFm2", "TOF m^2 vs Pt", HistType::kTH2F, {ptAxis, {400, 0.0, 10.0, "m^2"}});
    QA_p.add("histNClusterTPC", "Number of Clusters in TPC vs Pt (p)", HistType::kTH2F, {ptAxis, {80, 0.0, 160.0, "nCluster"}});
    QA_p.add("histNClusterITS", "Number of Clusters in ITS vs Pt (p)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    QA_p.add("histNClusterITSib", "Number of Clusters in ib of ITS vs Pt (p)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    QA_p.add("histTPCnClsFindable", "Findable TPC clusters (p)", HistType::kTH2F, {ptAxis, {200, 0.0, 200.0, "nCluster"}});
    QA_p.add("histTPCnClsFindableMinusFound", "TPC Clusters: Findable - Found (p)", HistType::kTH2F, {ptAxis, {60, 0.0, 60.0, "nCluster"}});
    QA_p.add("histTPCnClsFindableMinusCrossedRows", "TPC Clusters: Findable - crossed rows (p)", HistType::kTH2F, {ptAxis, {60, 0.0, 60.0, "nCluster"}});
    QA_p.add("histTPCnClsShared", "Number of shared TPC clusters (p)", HistType::kTH2F, {ptAxis, {70, 0.0, 70.0, "nCluster"}});
    QA_p.add("histTPCCrossedRowsOverFindableCls", "Ratio crossed rows over findable clusters (p)", HistType::kTH1F, {{100, 0., 2.0, "Crossed Rows / Findable Cls"}});
    QA_p.add("histTPCFoundOverFindable", "Ratio of found over findable clusters (p)", HistType::kTH1F, {{100, 0., 2.0, "Found Cls / Findable Cls"}});
    QA_p.add("histTPCFractionSharedCls", "Fraction of shared TPC clusters (p)", HistType::kTH1F, {{100, -2.0, 2.0, "Shared Cls"}});
    QA_p.add("histChi2TPC", "chi^2 TPC vs Pt (p)", HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    QA_p.add("histChi2ITS", "chi^2 ITS vs Pt (p)", HistType::kTH2F, {ptAxis, {500, 0.0, 50.0, "chi^2"}});
    QA_p.add("histChi2TOF", "chi^2 TOF vs Pt (p)", HistType::kTH2F, {ptAxis, {75, 0.0, 15.0, "chi^2"}});
    QA_p.add("histEtaWithOverFlow", "Pseudorapidity 0 - 105%% centrality (p)", HistType::kTH1F, {etaAxis});
    QA_p.add("histEta", "Pseudorapidity with centrality cut (p)", HistType::kTH1F, {etaAxis});
    QA_p.add("histEta_cent", "Pseudorapidity vs Centrality (p)", HistType::kTH2F, {centralityAxis_extended, etaAxis});
    QA_p.add("histTrackLength", "Track length (p)", HistType::kTH1F, {{350, 0., 700., "length (cm)"}});

    // QA antiproton
    QA_antip.add("histTpcSignalData", "Specific energy loss (antip)", HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}});
    QA_antip.add("histTofSignalData", "TOF signal (antip)", HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    QA_antip.add("histDcaVsPtData", "dcaXY vs Pt (antip)", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    QA_antip.add("histDcaZVsPtData", "dcaZ vs Pt (antip)", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    QA_antip.add("histTOFm2", "TOF m^2 vs Pt", HistType::kTH2F, {ptAxis, {400, 0.0, 10.0, "m^2"}});
    QA_antip.add("histNClusterTPC", "Number of Clusters in TPC vs Pt (antip)", HistType::kTH2F, {ptAxis, {80, 0.0, 160.0, "nCluster"}});
    QA_antip.add("histNClusterITS", "Number of Clusters in ITS vs Pt (antip)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    QA_antip.add("histNClusterITSib", "Number of Clusters in ib of ITS vs Pt (antip)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    QA_antip.add("histTPCnClsFindable", "Findable TPC clusters (antip)", HistType::kTH2F, {ptAxis, {200, 0.0, 200.0, "nCluster"}});
    QA_antip.add("histTPCnClsFindableMinusFound", "TPC Clusters: Findable - Found (antip)", HistType::kTH2F, {ptAxis, {60, 0.0, 60.0, "nCluster"}});
    QA_antip.add("histTPCnClsFindableMinusCrossedRows", "TPC Clusters: Findable - crossed rows (antip)", HistType::kTH2F, {ptAxis, {60, 0.0, 60.0, "nCluster"}});
    QA_antip.add("histTPCnClsShared", "Number of shared TPC clusters (antip)", HistType::kTH2F, {ptAxis, {70, 0.0, 70.0, "nCluster"}});
    QA_antip.add("histTPCCrossedRowsOverFindableCls", "Ratio crossed rows over findable clusters (antip)", HistType::kTH1F, {{100, 0., 2.0, "Crossed Rows / Findable Cls"}});
    QA_antip.add("histTPCFoundOverFindable", "Ratio of found over findable clusters (antip)", HistType::kTH1F, {{100, 0., 2.0, "Found Cls / Findable Cls"}});
    QA_antip.add("histTPCFractionSharedCls", "Fraction of shared TPC clusters (antip)", HistType::kTH1F, {{100, -2.0, 2.0, "Shared Cls"}});
    QA_antip.add("histChi2TPC", "chi^2 TPC vs Pt (antip)", HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    QA_antip.add("histChi2ITS", "chi^2 ITS vs Pt (antip)", HistType::kTH2F, {ptAxis, {500, 0.0, 50.0, "chi^2"}});
    QA_antip.add("histChi2TOF", "chi^2 TOF vs Pt (antip)", HistType::kTH2F, {ptAxis, {75, 0.0, 15.0, "chi^2"}});
    QA_antip.add("histEtaWithOverFlow", "Pseudorapidity 0 - 105%% centrality (antip)", HistType::kTH1F, {etaAxis});
    QA_antip.add("histEta", "Pseudorapidity with centrality cut (antip)", HistType::kTH1F, {etaAxis});
    QA_antip.add("histEta_cent", "Pseudorapidity vs Centrality (antip)", HistType::kTH2F, {centralityAxis_extended, etaAxis});
    QA_antip.add("histTrackLength", "Track length (antip)", HistType::kTH1F, {{350, 0., 700., "length (cm)"}});


    // QA deuteron
    QA_d.add("histTpcSignalData", "Specific energy loss (d)", HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}});
    QA_d.add("histTofSignalData", "TOF signal (d)", HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    QA_d.add("histDcaVsPtData", "dcaXY vs Pt (d)", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    QA_d.add("histDcaZVsPtData", "dcaZ vs Pt (d)", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    QA_d.add("histTOFm2", "TOF m^2 vs Pt", HistType::kTH2F, {ptAxis, {400, 0.0, 10.0, "m^2"}});
    QA_d.add("histNClusterTPC", "Number of Clusters in TPC vs Pt (d)", HistType::kTH2F, {ptAxis, {80, 0.0, 160.0, "nCluster"}});
    QA_d.add("histNClusterITS", "Number of Clusters in ITS vs Pt (d)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    QA_d.add("histNClusterITSib", "Number of Clusters in ib of ITS vs Pt (d)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    QA_d.add("histTPCnClsFindable", "Findable TPC clusters (d)", HistType::kTH2F, {ptAxis, {200, 0.0, 200.0, "nCluster"}});
    QA_d.add("histTPCnClsFindableMinusFound", "TPC Clusters: Findable - Found (d)", HistType::kTH2F, {ptAxis, {60, 0.0, 60.0, "nCluster"}});
    QA_d.add("histTPCnClsFindableMinusCrossedRows", "TPC Clusters: Findable - crossed rows (d)", HistType::kTH2F, {ptAxis, {60, 0.0, 60.0, "nCluster"}});
    QA_d.add("histTPCnClsShared", "Number of shared TPC clusters (d)", HistType::kTH2F, {ptAxis, {70, 0.0, 70.0, "nCluster"}});
    QA_d.add("histTPCCrossedRowsOverFindableCls", "Ratio crossed rows over findable clusters (d)", HistType::kTH1F, {{100, 0., 2.0, "Crossed Rows / Findable Cls"}});
    QA_d.add("histTPCFoundOverFindable", "Ratio of found over findable clusters (d)", HistType::kTH1F, {{100, 0., 2.0, "Found Cls / Findable Cls"}});
    QA_d.add("histTPCFractionSharedCls", "Fraction of shared TPC clusters (d)", HistType::kTH1F, {{100, -2.0, 2.0, "Shared Cls"}});
    QA_d.add("histChi2TPC", "chi^2 TPC vs Pt (d)", HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    QA_d.add("histChi2ITS", "chi^2 ITS vs Pt (d)", HistType::kTH2F, {ptAxis, {500, 0.0, 50.0, "chi^2"}});
    QA_d.add("histChi2TOF", "chi^2 TOF vs Pt (d)", HistType::kTH2F, {ptAxis, {75, 0.0, 15.0, "chi^2"}});
    QA_d.add("histEtaWithOverFlow", "Pseudorapidity 0 - 105%% centrality (d)", HistType::kTH1F, {etaAxis});
    QA_d.add("histEta", "Pseudorapidity with centrality cut (d)", HistType::kTH1F, {etaAxis});
    QA_d.add("histEta_cent", "Pseudorapidity vs Centrality (d)", HistType::kTH2F, {centralityAxis_extended, etaAxis});
    QA_d.add("histTrackLength", "Track length (d)", HistType::kTH1F, {{350, 0., 700., "length (cm)"}});

    // QA antideuteron
    QA_antid.add("histTpcSignalData", "Specific energy loss (antid)", HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}});
    QA_antid.add("histTofSignalData", "TOF signal (antid)", HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    QA_antid.add("histDcaVsPtData", "dcaXY vs Pt (antid)", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    QA_antid.add("histDcaZVsPtData", "dcaZ vs Pt (antid)", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    QA_antid.add("histTOFm2", "TOF m^2 vs Pt", HistType::kTH2F, {ptAxis, {400, 0.0, 10.0, "m^2"}});
    QA_antid.add("histNClusterTPC", "Number of Clusters in TPC vs Pt (antid)", HistType::kTH2F, {ptAxis, {80, 0.0, 160.0, "nCluster"}});
    QA_antid.add("histNClusterITS", "Number of Clusters in ITS vs Pt (antid)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    QA_antid.add("histNClusterITSib", "Number of Clusters in ib of ITS vs Pt (antid)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    QA_antid.add("histTPCnClsFindable", "Findable TPC clusters (antid)", HistType::kTH2F, {ptAxis, {200, 0.0, 200.0, "nCluster"}});
    QA_antid.add("histTPCnClsFindableMinusFound", "TPC Clusters: Findable - Found (antid)", HistType::kTH2F, {ptAxis, {60, 0.0, 60.0, "nCluster"}});
    QA_antid.add("histTPCnClsFindableMinusCrossedRows", "TPC Clusters: Findable - crossed rows (antid)", HistType::kTH2F, {ptAxis, {60, 0.0, 60.0, "nCluster"}});
    QA_antid.add("histTPCnClsShared", "Number of shared TPC clusters (antid)", HistType::kTH2F, {ptAxis, {70, 0.0, 70.0, "nCluster"}});
    QA_antid.add("histTPCCrossedRowsOverFindableCls", "Ratio crossed rows over findable clusters (antid)", HistType::kTH1F, {{100, 0., 2.0, "Crossed Rows / Findable Cls"}});
    QA_antid.add("histTPCFoundOverFindable", "Ratio of found over findable clusters (antid)", HistType::kTH1F, {{100, 0., 2.0, "Found Cls / Findable Cls"}});
    QA_antid.add("histTPCFractionSharedCls", "Fraction of shared TPC clusters (antid)", HistType::kTH1F, {{100, -2.0, 2.0, "Shared Cls"}});
    QA_antid.add("histChi2TPC", "chi^2 TPC vs Pt (antid)", HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    QA_antid.add("histChi2ITS", "chi^2 ITS vs Pt (antid)", HistType::kTH2F, {ptAxis, {500, 0.0, 50.0, "chi^2"}});
    QA_antid.add("histChi2TOF", "chi^2 TOF vs Pt (antid)", HistType::kTH2F, {ptAxis, {75, 0.0, 15.0, "chi^2"}});
    QA_antid.add("histEtaWithOverFlow", "Pseudorapidity 0 - 105%% centrality (antid)", HistType::kTH1F, {etaAxis});
    QA_antid.add("histEta", "Pseudorapidity with centrality cut (antid)", HistType::kTH1F, {etaAxis});
    QA_antid.add("histEta_cent", "Pseudorapidity vs Centrality (antid)", HistType::kTH2F, {centralityAxis_extended, etaAxis});
    QA_antid.add("histTrackLength", "Track length (antid)", HistType::kTH1F, {{350, 0., 700., "length (cm)"}});


    // QA triton
    QA_t.add("histTpcSignalData", "Specific energy loss (t)", HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}});
    QA_t.add("histTofSignalData", "TOF signal (t)", HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    QA_t.add("histDcaVsPtData", "dcaXY vs Pt (t)", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    QA_t.add("histDcaZVsPtData", "dcaZ vs Pt (t)", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    QA_t.add("histTOFm2", "TOF m^2 vs Pt", HistType::kTH2F, {ptAxis, {400, 0.0, 10.0, "m^2"}});
    QA_t.add("histNClusterTPC", "Number of Clusters in TPC vs Pt (t)", HistType::kTH2F, {ptAxis, {80, 0.0, 160.0, "nCluster"}});
    QA_t.add("histNClusterITS", "Number of Clusters in ITS vs Pt (t)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    QA_t.add("histNClusterITSib", "Number of Clusters in ib of ITS vs Pt (t)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    QA_t.add("histTPCnClsFindable", "Findable TPC clusters (t)", HistType::kTH2F, {ptAxis, {200, 0.0, 200.0, "nCluster"}});
    QA_t.add("histTPCnClsFindableMinusFound", "TPC Clusters: Findable - Found (t)", HistType::kTH2F, {ptAxis, {60, 0.0, 60.0, "nCluster"}});
    QA_t.add("histTPCnClsFindableMinusCrossedRows", "TPC Clusters: Findable - crossed rows (t)", HistType::kTH2F, {ptAxis, {60, 0.0, 60.0, "nCluster"}});
    QA_t.add("histTPCnClsShared", "Number of shared TPC clusters (t)", HistType::kTH2F, {ptAxis, {70, 0.0, 70.0, "nCluster"}});
    QA_t.add("histTPCCrossedRowsOverFindableCls", "Ratio crossed rows over findable clusters (t)", HistType::kTH1F, {{100, 0., 2.0, "Crossed Rows / Findable Cls"}});
    QA_t.add("histTPCFoundOverFindable", "Ratio of found over findable clusters (t)", HistType::kTH1F, {{100, 0., 2.0, "Found Cls / Findable Cls"}});
    QA_t.add("histTPCFractionSharedCls", "Fraction of shared TPC clusters (t)", HistType::kTH1F, {{100, -2.0, 2.0, "Shared Cls"}});
    QA_t.add("histChi2TPC", "chi^2 TPC vs Pt (t)", HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    QA_t.add("histChi2ITS", "chi^2 ITS vs Pt (t)", HistType::kTH2F, {ptAxis, {500, 0.0, 50.0, "chi^2"}});
    QA_t.add("histChi2TOF", "chi^2 TOF vs Pt (t)", HistType::kTH2F, {ptAxis, {75, 0.0, 15.0, "chi^2"}});
    QA_t.add("histEtaWithOverFlow", "Pseudorapidity 0 - 105%% centrality (t)", HistType::kTH1F, {etaAxis});
    QA_t.add("histEta", "Pseudorapidity with centrality cut (t)", HistType::kTH1F, {etaAxis});
    QA_t.add("histEta_cent", "Pseudorapidity vs Centrality (t)", HistType::kTH2F, {centralityAxis_extended, etaAxis});
    QA_t.add("histTrackLength", "Track length (t)", HistType::kTH1F, {{350, 0., 700., "length (cm)"}});

    // QA antitriton
    QA_antit.add("histTpcSignalData", "Specific energy loss (antit)", HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}});
    QA_antit.add("histTofSignalData", "TOF signal (antit)", HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    QA_antit.add("histDcaVsPtData", "dcaXY vs Pt (antit)", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    QA_antit.add("histDcaZVsPtData", "dcaZ vs Pt (antit)", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    QA_antit.add("histTOFm2", "TOF m^2 vs Pt", HistType::kTH2F, {ptAxis, {400, 0.0, 10.0, "m^2"}});
    QA_antit.add("histNClusterTPC", "Number of Clusters in TPC vs Pt (antit)", HistType::kTH2F, {ptAxis, {80, 0.0, 160.0, "nCluster"}});
    QA_antit.add("histNClusterITS", "Number of Clusters in ITS vs Pt (antit)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    QA_antit.add("histNClusterITSib", "Number of Clusters in ib of ITS vs Pt (antit)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    QA_antit.add("histTPCnClsFindable", "Findable TPC clusters (antit)", HistType::kTH2F, {ptAxis, {200, 0.0, 200.0, "nCluster"}});
    QA_antit.add("histTPCnClsFindableMinusFound", "TPC Clusters: Findable - Found (antit)", HistType::kTH2F, {ptAxis, {60, 0.0, 60.0, "nCluster"}});
    QA_antit.add("histTPCnClsFindableMinusCrossedRows", "TPC Clusters: Findable - crossed rows (antit)", HistType::kTH2F, {ptAxis, {60, 0.0, 60.0, "nCluster"}});
    QA_antit.add("histTPCnClsShared", "Number of shared TPC clusters (antit)", HistType::kTH2F, {ptAxis, {70, 0.0, 70.0, "nCluster"}});
    QA_antit.add("histTPCCrossedRowsOverFindableCls", "Ratio crossed rows over findable clusters (antit)", HistType::kTH1F, {{100, 0., 2.0, "Crossed Rows / Findable Cls"}});
    QA_antit.add("histTPCFoundOverFindable", "Ratio of found over findable clusters (antit)", HistType::kTH1F, {{100, 0., 2.0, "Found Cls / Findable Cls"}});
    QA_antit.add("histTPCFractionSharedCls", "Fraction of shared TPC clusters (antit)", HistType::kTH1F, {{100, -2.0, 2.0, "Shared Cls"}});
    QA_antit.add("histChi2TPC", "chi^2 TPC vs Pt (antit)", HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    QA_antit.add("histChi2ITS", "chi^2 ITS vs Pt (antit)", HistType::kTH2F, {ptAxis, {500, 0.0, 50.0, "chi^2"}});
    QA_antit.add("histChi2TOF", "chi^2 TOF vs Pt (antip)", HistType::kTH2F, {ptAxis, {75, 0.0, 15.0, "chi^2"}});
    QA_antit.add("histEtaWithOverFlow", "Pseudorapidity 0 - 105%% centrality (antit)", HistType::kTH1F, {etaAxis});
    QA_antit.add("histEta", "Pseudorapidity with centrality cut (antit)", HistType::kTH1F, {etaAxis});
    QA_antit.add("histEta_cent", "Pseudorapidity vs Centrality (antit)", HistType::kTH2F, {centralityAxis_extended, etaAxis});
    QA_antit.add("histTrackLength", "Track length (antit)", HistType::kTH1F, {{350, 0., 700., "length (cm)"}});


    // QA Helium-3
    QA_He3.add("histTpcSignalData", "Specific energy loss (He3)", HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}});
    QA_He3.add("histTofSignalData", "TOF signal (He3)", HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    QA_He3.add("histDcaVsPtData", "dcaXY vs Pt (He3)", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    QA_He3.add("histDcaZVsPtData", "dcaZ vs Pt (He3)", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    QA_He3.add("histTOFm2", "TOF m^2 vs Pt", HistType::kTH2F, {ptAxis, {400, 0.0, 10.0, "m^2"}});
    QA_He3.add("histNClusterTPC", "Number of Clusters in TPC vs Pt (He3)", HistType::kTH2F, {ptAxis, {80, 0.0, 160.0, "nCluster"}});
    QA_He3.add("histNClusterITS", "Number of Clusters in ITS vs Pt (He3)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    QA_He3.add("histNClusterITSib", "Number of Clusters in ib of ITS vs Pt (He3)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    QA_He3.add("histTPCnClsFindable", "Findable TPC clusters (He3)", HistType::kTH2F, {ptAxis, {200, 0.0, 200.0, "nCluster"}});
    QA_He3.add("histTPCnClsFindableMinusFound", "TPC Clusters: Findable - Found (He3)", HistType::kTH2F, {ptAxis, {60, 0.0, 60.0, "nCluster"}});
    QA_He3.add("histTPCnClsFindableMinusCrossedRows", "TPC Clusters: Findable - crossed rows (He3)", HistType::kTH2F, {ptAxis, {60, 0.0, 60.0, "nCluster"}});
    QA_He3.add("histTPCnClsShared", "Number of shared TPC clusters (He3)", HistType::kTH2F, {ptAxis, {70, 0.0, 70.0, "nCluster"}});
    QA_He3.add("histTPCCrossedRowsOverFindableCls", "Ratio crossed rows over findable clusters (He3)", HistType::kTH1F, {{100, 0., 2.0, "Crossed Rows / Findable Cls"}});
    QA_He3.add("histTPCFoundOverFindable", "Ratio of found over findable clusters (He3)", HistType::kTH1F, {{100, 0., 2.0, "Found Cls / Findable Cls"}});
    QA_He3.add("histTPCFractionSharedCls", "Fraction of shared TPC clusters (He3)", HistType::kTH1F, {{100, -2.0, 2.0, "Shared Cls"}});
    QA_He3.add("histChi2TPC", "chi^2 TPC vs Pt (He3)", HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    QA_He3.add("histChi2ITS", "chi^2 ITS vs Pt (He3)", HistType::kTH2F, {ptAxis, {500, 0.0, 50.0, "chi^2"}});
    QA_He3.add("histChi2TOF", "chi^2 TOF vs Pt (He3)", HistType::kTH2F, {ptAxis, {75, 0.0, 15.0, "chi^2"}});
    QA_He3.add("histEtaWithOverFlow", "Pseudorapidity 0 - 105%% centrality (He3)", HistType::kTH1F, {etaAxis});
    QA_He3.add("histEta", "Pseudorapidity with centrality cut (He3)", HistType::kTH1F, {etaAxis});
    QA_He3.add("histEta_cent", "Pseudorapidity vs Centrality (He3)", HistType::kTH2F, {centralityAxis_extended, etaAxis});
    QA_He3.add("histTrackLength", "Track length (He3)", HistType::kTH1F, {{350, 0., 700., "length (cm)"}});

    // QA antiHelium-3
    QA_antiHe3.add("histTpcSignalData", "Specific energy loss (antiHe3)", HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}});
    QA_antiHe3.add("histTofSignalData", "TOF signal (antiHe3)", HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    QA_antiHe3.add("histDcaVsPtData", "dcaXY vs Pt (antiHe3)", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    QA_antiHe3.add("histDcaZVsPtData", "dcaZ vs Pt (antiHe3)", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    QA_antiHe3.add("histTOFm2", "TOF m^2 vs Pt", HistType::kTH2F, {ptAxis, {400, 0.0, 10.0, "m^2"}});
    QA_antiHe3.add("histNClusterTPC", "Number of Clusters in TPC vs Pt (antiHe3)", HistType::kTH2F, {ptAxis, {80, 0.0, 160.0, "nCluster"}});
    QA_antiHe3.add("histNClusterITS", "Number of Clusters in ITS vs Pt (antiHe3)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    QA_antiHe3.add("histNClusterITSib", "Number of Clusters in ib of ITS vs Pt (antiHe3)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    QA_antiHe3.add("histTPCnClsFindable", "Findable TPC clusters (antiHe3)", HistType::kTH2F, {ptAxis, {200, 0.0, 200.0, "nCluster"}});
    QA_antiHe3.add("histTPCnClsFindableMinusFound", "TPC Clusters: Findable - Found (antiHe3)", HistType::kTH2F, {ptAxis, {60, 0.0, 60.0, "nCluster"}});
    QA_antiHe3.add("histTPCnClsFindableMinusCrossedRows", "TPC Clusters: Findable - crossed rows (antiHe3)", HistType::kTH2F, {ptAxis, {60, 0.0, 60.0, "nCluster"}});
    QA_antiHe3.add("histTPCnClsShared", "Number of shared TPC clusters (antiHe3)", HistType::kTH2F, {ptAxis, {70, 0.0, 70.0, "nCluster"}});
    QA_antiHe3.add("histTPCCrossedRowsOverFindableCls", "Ratio crossed rows over findable clusters (antiHe3)", HistType::kTH1F, {{100, 0., 2.0, "Crossed Rows / Findable Cls"}});
    QA_antiHe3.add("histTPCFoundOverFindable", "Ratio of found over findable clusters (antiHe3)", HistType::kTH1F, {{100, 0., 2.0, "Found Cls / Findable Cls"}});
    QA_antiHe3.add("histTPCFractionSharedCls", "Fraction of shared TPC clusters (antiHe3)", HistType::kTH1F, {{100, -2.0, 2.0, "Shared Cls"}});
    QA_antiHe3.add("histChi2TPC", "chi^2 TPC vs Pt (antiHe3)", HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    QA_antiHe3.add("histChi2ITS", "chi^2 ITS vs Pt (antiHe3)", HistType::kTH2F, {ptAxis, {500, 0.0, 50.0, "chi^2"}});
    QA_antiHe3.add("histChi2TOF", "chi^2 TOF vs Pt (antiHe3)", HistType::kTH2F, {ptAxis, {75, 0.0, 15.0, "chi^2"}});
    QA_antiHe3.add("histEtaWithOverFlow", "Pseudorapidity 0 - 105%% centrality (antiHe3)", HistType::kTH1F, {etaAxis});
    QA_antiHe3.add("histEta", "Pseudorapidity with centrality cut (antiHe3)", HistType::kTH1F, {etaAxis});
    QA_antiHe3.add("histEta_cent", "Pseudorapidity vs Centrality (antiHe3)", HistType::kTH2F, {centralityAxis_extended, etaAxis});
    QA_antiHe3.add("histTrackLength", "Track length (antiHe3)", HistType::kTH1F, {{350, 0., 700., "length (cm)"}});


    // QA Helium-4
    QA_He4.add("histTpcSignalData", "Specific energy loss (He4)", HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}});
    QA_He4.add("histTofSignalData", "TOF signal (He4)", HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    QA_He4.add("histDcaVsPtData", "dcaXY vs Pt (He4)", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    QA_He4.add("histDcaZVsPtData", "dcaZ vs Pt (He4)", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    QA_He4.add("histTOFm2", "TOF m^2 vs Pt", HistType::kTH2F, {ptAxis, {400, 0.0, 10.0, "m^2"}});
    QA_He4.add("histNClusterTPC", "Number of Clusters in TPC vs Pt (He4)", HistType::kTH2F, {ptAxis, {80, 0.0, 160.0, "nCluster"}});
    QA_He4.add("histNClusterITS", "Number of Clusters in ITS vs Pt (He4)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    QA_He4.add("histNClusterITSib", "Number of Clusters in ib of ITS vs Pt (He4)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    QA_He4.add("histTPCnClsFindable", "Findable TPC clusters (He4)", HistType::kTH2F, {ptAxis, {200, 0.0, 200.0, "nCluster"}});
    QA_He4.add("histTPCnClsFindableMinusFound", "TPC Clusters: Findable - Found (He4)", HistType::kTH2F, {ptAxis, {60, 0.0, 60.0, "nCluster"}});
    QA_He4.add("histTPCnClsFindableMinusCrossedRows", "TPC Clusters: Findable - crossed rows (He4)", HistType::kTH2F, {ptAxis, {60, 0.0, 60.0, "nCluster"}});
    QA_He4.add("histTPCnClsShared", "Number of shared TPC clusters (He4)", HistType::kTH2F, {ptAxis, {70, 0.0, 70.0, "nCluster"}});
    QA_He4.add("histTPCCrossedRowsOverFindableCls", "Ratio crossed rows over findable clusters (He4)", HistType::kTH1F, {{100, 0., 2.0, "Crossed Rows / Findable Cls"}});
    QA_He4.add("histTPCFoundOverFindable", "Ratio of found over findable clusters (He4)", HistType::kTH1F, {{100, 0., 2.0, "Found Cls / Findable Cls"}});
    QA_He4.add("histTPCFractionSharedCls", "Fraction of shared TPC clusters (He4)", HistType::kTH1F, {{100, -2.0, 2.0, "Shared Cls"}});
    QA_He4.add("histChi2TPC", "chi^2 TPC vs Pt (He4)", HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    QA_He4.add("histChi2ITS", "chi^2 ITS vs Pt (He4)", HistType::kTH2F, {ptAxis, {500, 0.0, 50.0, "chi^2"}});
    QA_He4.add("histChi2TOF", "chi^2 TOF vs Pt (He4)", HistType::kTH2F, {ptAxis, {75, 0.0, 15.0, "chi^2"}});
    QA_He4.add("histEtaWithOverFlow", "Pseudorapidity 0 - 105%% centrality (He4)", HistType::kTH1F, {etaAxis});
    QA_He4.add("histEta", "Pseudorapidity with centrality cut (He4)", HistType::kTH1F, {etaAxis});
    QA_He4.add("histEta_cent", "Pseudorapidity vs Centrality (He4)", HistType::kTH2F, {centralityAxis_extended, etaAxis});
    QA_He4.add("histTrackLength", "Track length (He4)", HistType::kTH1F, {{350, 0., 700., "length (cm)"}});

    // QA antiHelium-4
    QA_antiHe4.add("histTpcSignalData", "Specific energy loss (antiHe4)", HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}});
    QA_antiHe4.add("histTofSignalData", "TOF signal (antiHe4)", HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    QA_antiHe4.add("histDcaVsPtData", "dcaXY vs Pt (antiHe4)", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    QA_antiHe4.add("histDcaZVsPtData", "dcaZ vs Pt (antiHe4)", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    QA_antiHe4.add("histTOFm2", "TOF m^2 vs Pt", HistType::kTH2F, {ptAxis, {400, 0.0, 10.0, "m^2"}});
    QA_antiHe4.add("histNClusterTPC", "Number of Clusters in TPC vs Pt (antiHe4)", HistType::kTH2F, {ptAxis, {80, 0.0, 160.0, "nCluster"}});
    QA_antiHe4.add("histNClusterITS", "Number of Clusters in ITS vs Pt (antiHe4)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    QA_antiHe4.add("histNClusterITSib", "Number of Clusters in ib of ITS vs Pt (antiHe4)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    QA_antiHe4.add("histTPCnClsFindable", "Findable TPC clusters (antiHe4)", HistType::kTH2F, {ptAxis, {200, 0.0, 200.0, "nCluster"}});
    QA_antiHe4.add("histTPCnClsFindableMinusFound", "TPC Clusters: Findable - Found (antiHe4)", HistType::kTH2F, {ptAxis, {60, 0.0, 60.0, "nCluster"}});
    QA_antiHe4.add("histTPCnClsFindableMinusCrossedRows", "TPC Clusters: Findable - crossed rows (antiHe4)", HistType::kTH2F, {ptAxis, {60, 0.0, 60.0, "nCluster"}});
    QA_antiHe4.add("histTPCnClsShared", "Number of shared TPC clusters (antiHe4)", HistType::kTH2F, {ptAxis, {70, 0.0, 70.0, "nCluster"}});
    QA_antiHe4.add("histTPCCrossedRowsOverFindableCls", "Ratio crossed rows over findable clusters (antiHe4)", HistType::kTH1F, {{100, 0., 2.0, "Crossed Rows / Findable Cls"}});
    QA_antiHe4.add("histTPCFoundOverFindable", "Ratio of found over findable clusters (antiHe4)", HistType::kTH1F, {{100, 0., 2.0, "Found Cls / Findable Cls"}});
    QA_antiHe4.add("histTPCFractionSharedCls", "Fraction of shared TPC clusters (antiHe4)", HistType::kTH1F, {{100, -2.0, 2.0, "Shared Cls"}});
    QA_antiHe4.add("histChi2TPC", "chi^2 TPC vs Pt (antiHe4)", HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    QA_antiHe4.add("histChi2ITS", "chi^2 ITS vs Pt (antiHe4)", HistType::kTH2F, {ptAxis, {500, 0.0, 50.0, "chi^2"}});
    QA_antiHe4.add("histChi2TOF", "chi^2 TOF vs Pt (antiHe4)", HistType::kTH2F, {ptAxis, {75, 0.0, 15.0, "chi^2"}});
    QA_antiHe4.add("histEtaWithOverFlow", "Pseudorapidity 0 - 105%% centrality (antiHe4)", HistType::kTH1F, {etaAxis});
    QA_antiHe4.add("histEta", "Pseudorapidity with centrality cut (antiHe4)", HistType::kTH1F, {etaAxis});
    QA_antiHe4.add("histEta_cent", "Pseudorapidity vs Centrality (antiHe4)", HistType::kTH2F, {centralityAxis_extended, etaAxis});
    QA_antiHe4.add("histTrackLength", "Track length (antiHe4)", HistType::kTH1F, {{350, 0., 700., "length (cm)"}});

  }

  // Configurables
  Configurable<bool> enable_proton_processing{"enable_proton_processing", false, "0: disabled, 1: enabled"};
  Configurable<bool> enable_deuteron_processing{"enable_deuteron_processing", false, "0: disabled, 1: enabled"};
  Configurable<bool> enable_triton_processing{"enable_triton_processing", false, "0: disabled, 1: enabled"};
  Configurable<bool> enable_Helium3_processing{"enable_Helium3_processing", false, "0: disabled, 1: enabled"};
  Configurable<bool> enable_Helium4_processing{"enable_Helium4_processing", false, "0: disabled, 1: enabled"};
  
  Configurable<bool> event_selection_sel8{"event_selection_sel8", true, "0: disabled, 1: enabled"};
  Configurable<float> yMin{"yMin", -0.5, "Maximum rapidity"};
  Configurable<float> yMax{"yMax", 0.5, "Minimum rapidity"};
  Configurable<float> pTmin{"pTmin", 0.1f, "min pT"};
  Configurable<float> pTmax{"pTmax", 1e+10f, "max pT"};
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<bool> enable_Centrality_cut{"enable_Centrality_cut", false, "Centrality for histEta; 0: disabled, 1: enabled"};
  Configurable<float> minCentrality{"minCentrality", 0.0, "min Centrality used (histEta)"};
  Configurable<float> maxCentrality{"maxCentrality", 80.0, "max Centrality used (histEta)"};
  Configurable<bool> enable_PVcontributor{"enable_PVcontributor", false, "0: disabled, 1: enabled"};

  Configurable<float> nsigmacut{"nsigmacut", 3.0, "absolute value of the Nsigma cut"};

  Configurable<bool> custom_Track_selection{"custom_Track_selection", false, "0: disabled, 1: enabled"};
  // custom Track selection cuts
  Configurable<float> minReqClusterITS{"minReqClusterITS", 0.0, "min number of clusters required in ITS"};                            // ITS_nCls
  Configurable<float> minReqClusterITSib{"minReqClusterITSib", 0.0, "min number of clusters required in ITS inner barrel"};           // ITS_nCls
  Configurable<float> minTPCnClsFound{"minTPCnClsFound", 0.0f, "min number of crossed rows TPC"};                                     // TPC_nCls_found
  Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 0.0f, "min number of crossed rows TPC"};                               // TPC_nCls_crossed_Rows
  Configurable<float> minRatioCrossedRowsTPC{"minRatioCrossedRowsTPC", 0.0f, "min ratio of crossed rows over findable clusters TPC"}; // TPC_crossed_Rows_over_findable_Cls_min
  Configurable<float> maxRatioCrossedRowsTPC{"maxRatioCrossedRowsTPC", 2.0f, "max ratio of crossed rows over findable clusters TPC"}; // TPC_crossed_Rows_over_findable_Cls_max
  Configurable<float> maxChi2ITS{"maxChi2ITS", 100.0f, "max chi2 per cluster ITS"};
  Configurable<float> maxChi2TPC{"maxChi2TPC", 100.0f, "max chi2 per cluster TPC"};
  Configurable<float> maxDCA_XY{"maxDCA_XY", 10.0f, "max DCA to vertex xy"};
  Configurable<float> maxDCA_Z{"maxDCA_Z", 10.0f, "max DCA to vertex z"};
  Configurable<int> lastRequiredTrdCluster{"lastRequiredTrdCluster", -1, "Last cluster to required in TRD for track selection. -1 does not require any TRD cluster"};
  Configurable<float> minTrackLength{"minTrackLength", 0.0f, "min track length"};
  Configurable<float> maxTrackLength{"maxTrackLength", 999999.0f, "max track length"};
  Configurable<float> minTPCNClsFindable{"minTPCNClsFindable", 0.0f, "min findable TPC clusters for this track geometry"};
  Configurable<float> maxTPCNClsShared{"maxTPCNClsShared", 999999.0f, "max number of shared TPC clusters"};
  Configurable<float> maxChi2TOF{"maxChi2TOF", 100.0f, "max chi2 for the TOF track segment"};
  Configurable<float> minTPCFoundOverFindable{"minTPCFoundOverFindable", 0.0f, "min ratio of found over findable clusters TPC"};
  Configurable<float> maxTPCFoundOverFindable{"maxTPCFoundOverFindable", 2.0f, "max ratio of found over findable clusters TPC"};


  template <typename CollisionType, typename TracksType>
  void fillHistograms(const CollisionType& event, const TracksType& tracks)
  {

    if (event_selection_sel8 && event.sel8()) {
      QA_reg.fill(HIST("histRecVtxZData"), event.posZ());
      QA_reg.fill(HIST("histCentrality"), event.centFT0C());
    }

    if (!event_selection_sel8) {
      QA_reg.fill(HIST("histRecVtxZData"), event.posZ());
      QA_reg.fill(HIST("histCentrality"), event.centFT0C());
    }

    for (auto track : tracks) { // start loop over tracks

      if (event_selection_sel8 && !event.sel8()) {
        continue;
      }

      QA_reg.fill(HIST("histEtaWithOverFlow"), track.eta());
      QA_reg.fill(HIST("histEta_cent"), event.centFT0C(), track.eta());

      if ((event.centFT0C() > minCentrality) && (event.centFT0C() < maxCentrality) && enable_Centrality_cut) {
        QA_reg.fill(HIST("histEta"), track.eta());
      }

      float TPCnumberClsFound = track.tpcNClsFound();
      float TPC_nCls_Crossed_Rows = track.tpcNClsCrossedRows();
      float RatioCrossedRowsOverFindableTPC = track.tpcCrossedRowsOverFindableCls();
      float Chi2perClusterTPC = track.tpcChi2NCl();
      float Chi2perClusterITS = track.itsChi2NCl();

      if (track.sign() > 0) {
        QA_reg.fill(HIST("histDcaVsPtData_particle"), track.pt(), track.dcaXY());
        QA_reg.fill(HIST("histDcaZVsPtData_particle"), track.pt(), track.dcaZ());
      }

      if (track.sign() < 0) {
        QA_reg.fill(HIST("histDcaVsPtData_antiparticle"), track.pt(), track.dcaXY());
        QA_reg.fill(HIST("histDcaZVsPtData_antiparticle"), track.pt(), track.dcaZ());
      }

      if (custom_Track_selection && (TMath::Abs(track.dcaXY()) > maxDCA_XY || TMath::Abs(track.dcaZ()) > maxDCA_Z || TPCnumberClsFound < minTPCnClsFound || TPC_nCls_Crossed_Rows < minNCrossedRowsTPC || RatioCrossedRowsOverFindableTPC < minRatioCrossedRowsTPC || RatioCrossedRowsOverFindableTPC > maxRatioCrossedRowsTPC || Chi2perClusterTPC > maxChi2TPC || Chi2perClusterITS > maxChi2ITS || !(track.passedTPCRefit()) || !(track.passedITSRefit()) || (track.itsNClsInnerBarrel()) < minReqClusterITSib || (track.itsNCls()) < minReqClusterITS || track.length() < minTrackLength || track.length() > maxTrackLength || track.tpcNClsFindable() < minTPCNClsFindable || track.tpcNClsShared() > maxTPCNClsShared || track.tpcFoundOverFindableCls() < minTPCFoundOverFindable || track.tpcFoundOverFindableCls() > maxTPCFoundOverFindable)) {

        if (track.hasTOF() && track.tofChi2() > maxChi2TOF) continue;
        continue;
      }

      if (enable_PVcontributor && !(track.isPVContributor())) {
        continue;
      }

      // cut on rapidity
      TLorentzVector lorentzVector_proton{};
      TLorentzVector lorentzVector_deuteron{};
      TLorentzVector lorentzVector_triton{};
      TLorentzVector lorentzVector_He3{};
      TLorentzVector lorentzVector_He4{};

      lorentzVector_proton.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassProton);
      lorentzVector_deuteron.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassDeuteron);
      lorentzVector_triton.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassTriton);
      lorentzVector_He3.SetPtEtaPhiM(track.pt() * 2.0, track.eta(), track.phi(), constants::physics::MassHelium3);
      lorentzVector_He4.SetPtEtaPhiM(track.pt() * 2.0, track.eta(), track.phi(), constants::physics::MassAlpha);

      if (lorentzVector_proton.Rapidity() < yMin || lorentzVector_proton.Rapidity() > yMax ||
          lorentzVector_deuteron.Rapidity() < yMin || lorentzVector_deuteron.Rapidity() > yMax ||
          lorentzVector_triton.Rapidity() < yMin || lorentzVector_triton.Rapidity() > yMax ||
          lorentzVector_He3.Rapidity() < yMin || lorentzVector_He3.Rapidity() > yMax ||
          lorentzVector_He4.Rapidity() < yMin || lorentzVector_He4.Rapidity() > yMax) {
        continue;
      }

      if (track.pt() < pTmin || track.pt() > pTmax) continue;

      // fill QA histograms (all species)
      float nSigmaProton = track.tpcNSigmaPr();
      float nSigmaDeut = track.tpcNSigmaDe();
      float nSigmaTriton = track.tpcNSigmaTr();
      float nSigmaHe3 = track.tpcNSigmaHe();
      float nSigmaHe4 = track.tpcNSigmaAl();

      QA_reg.fill(HIST("histTpcSignalData"), track.tpcInnerParam() * track.sign(), track.tpcSignal());
      QA_reg.fill(HIST("histNClusterTPC"), track.pt(), track.tpcNClsCrossedRows());
      QA_reg.fill(HIST("histNClusterITS"), track.pt(), track.itsNCls());
      QA_reg.fill(HIST("histNClusterITSib"), track.pt(), track.itsNClsInnerBarrel());
      QA_reg.fill(HIST("histChi2TPC"), track.pt(), track.tpcChi2NCl());
      QA_reg.fill(HIST("histChi2ITS"), track.pt(), track.itsChi2NCl());
      QA_reg.fill(HIST("histTPCnClsFindable"), track.pt(), track.tpcNClsFindable());
      QA_reg.fill(HIST("histTPCnClsFindableMinusFound"), track.pt(), track.tpcNClsFindableMinusFound());
      QA_reg.fill(HIST("histTPCnClsFindableMinusCrossedRows"), track.pt(), track.tpcNClsFindableMinusCrossedRows());
      QA_reg.fill(HIST("histTPCnClsShared"), track.pt(), track.tpcNClsShared());
      QA_reg.fill(HIST("histChi2TOF"), track.pt(), track.tofChi2());
      QA_reg.fill(HIST("histTrackLength"), track.length());
      QA_reg.fill(HIST("histTPCCrossedRowsOverFindableCls"), track.tpcCrossedRowsOverFindableCls());
      QA_reg.fill(HIST("histTPCFoundOverFindable"), track.tpcFoundOverFindableCls());
      QA_reg.fill(HIST("histTPCFractionSharedCls"), track.tpcFractionSharedCls());

      if (track.hasTOF()) {

        if (track.hasTRD() && (lastRequiredTrdCluster > 0) && custom_Track_selection) {
          int lastLayer = 0;
          for (int l = 7; l >= 0; l--) {
            if (track.trdPattern() & (1 << l)) {
              lastLayer = l;
              break;
            }
          }
          if (lastLayer < lastRequiredTrdCluster) {
            continue;
          }
        }

        Float_t TOFmass2 = ((track.mass()) * (track.mass()));

        QA_reg.fill(HIST("histTOFm2"), track.pt(), TOFmass2);
        QA_reg.fill(HIST("histTofSignalData"), track.tpcInnerParam() * track.sign(), track.beta());
      }


      // fill QA histograms (proton)
      if (TMath::Abs(nSigmaProton) < nsigmacut && enable_proton_processing) {
        if (track.sign() > 0) {
          QA_p.fill(HIST("histTpcSignalData"), track.tpcInnerParam(), track.tpcSignal());
          QA_p.fill(HIST("histNClusterTPC"), track.pt(), track.tpcNClsCrossedRows());
          QA_p.fill(HIST("histNClusterITS"), track.pt(), track.itsNCls());
          QA_p.fill(HIST("histNClusterITSib"), track.pt(), track.itsNClsInnerBarrel());
          QA_p.fill(HIST("histChi2TPC"), track.pt(), track.tpcChi2NCl());
          QA_p.fill(HIST("histChi2ITS"), track.pt(), track.itsChi2NCl());
          QA_p.fill(HIST("histTPCnClsFindable"), track.pt(), track.tpcNClsFindable());
          QA_p.fill(HIST("histTPCnClsFindableMinusFound"), track.pt(), track.tpcNClsFindableMinusFound());
          QA_p.fill(HIST("histTPCnClsFindableMinusCrossedRows"), track.pt(), track.tpcNClsFindableMinusCrossedRows());
          QA_p.fill(HIST("histTPCnClsShared"), track.pt(), track.tpcNClsShared());
          QA_p.fill(HIST("histChi2TOF"), track.pt(), track.tofChi2());
          QA_p.fill(HIST("histTrackLength"), track.length());
          QA_p.fill(HIST("histTPCCrossedRowsOverFindableCls"), track.tpcCrossedRowsOverFindableCls());
          QA_p.fill(HIST("histTPCFoundOverFindable"), track.tpcFoundOverFindableCls());
          QA_p.fill(HIST("histTPCFractionSharedCls"), track.tpcFractionSharedCls());

          if (track.hasTOF()) {

            if (track.hasTRD() && (lastRequiredTrdCluster > 0) && custom_Track_selection) {
              int lastLayer = 0;
              for (int l = 7; l >= 0; l--) {
                if (track.trdPattern() & (1 << l)) {
                  lastLayer = l;
                  break;
                }
              }
              if (lastLayer < lastRequiredTrdCluster) {
                continue;
              }
            }

            Float_t TOFmass2 = ((track.mass()) * (track.mass()));

            QA_p.fill(HIST("histTOFm2"), track.pt(), TOFmass2);
            QA_p.fill(HIST("histTofSignalData"), track.tpcInnerParam() * track.sign(), track.beta());
          }
        }
        if (track.sign() < 0) {
          QA_antip.fill(HIST("histTpcSignalData"), track.tpcInnerParam(), track.tpcSignal());
          QA_antip.fill(HIST("histNClusterTPC"), track.pt(), track.tpcNClsCrossedRows());
          QA_antip.fill(HIST("histNClusterITS"), track.pt(), track.itsNCls());
          QA_antip.fill(HIST("histNClusterITSib"), track.pt(), track.itsNClsInnerBarrel());
          QA_antip.fill(HIST("histChi2TPC"), track.pt(), track.tpcChi2NCl());
          QA_antip.fill(HIST("histChi2ITS"), track.pt(), track.itsChi2NCl());
          QA_antip.fill(HIST("histTPCnClsFindable"), track.pt(), track.tpcNClsFindable());
          QA_antip.fill(HIST("histTPCnClsFindableMinusFound"), track.pt(), track.tpcNClsFindableMinusFound());
          QA_antip.fill(HIST("histTPCnClsFindableMinusCrossedRows"), track.pt(), track.tpcNClsFindableMinusCrossedRows());
          QA_antip.fill(HIST("histTPCnClsShared"), track.pt(), track.tpcNClsShared());
          QA_antip.fill(HIST("histChi2TOF"), track.pt(), track.tofChi2());
          QA_antip.fill(HIST("histTrackLength"), track.length());
          QA_antip.fill(HIST("histTPCCrossedRowsOverFindableCls"), track.tpcCrossedRowsOverFindableCls());
          QA_antip.fill(HIST("histTPCFoundOverFindable"), track.tpcFoundOverFindableCls());
          QA_antip.fill(HIST("histTPCFractionSharedCls"), track.tpcFractionSharedCls());

          if (track.hasTOF()) {

            if (track.hasTRD() && (lastRequiredTrdCluster > 0) && custom_Track_selection) {
              int lastLayer = 0;
              for (int l = 7; l >= 0; l--) {
                if (track.trdPattern() & (1 << l)) {
                  lastLayer = l;
                  break;
                }
              }
              if (lastLayer < lastRequiredTrdCluster) {
                continue;
              }
            }

            Float_t TOFmass2 = ((track.mass()) * (track.mass()));

            QA_antip.fill(HIST("histTOFm2"), track.pt(), TOFmass2);
            QA_antip.fill(HIST("histTofSignalData"), track.tpcInnerParam() * track.sign(), track.beta());
          }
        }
      }


      // fill QA histograms (deuteron)
      if (TMath::Abs(nSigmaDeut) < nsigmacut && enable_deuteron_processing) {
        if (track.sign() > 0) {
          QA_d.fill(HIST("histTpcSignalData"), track.tpcInnerParam(), track.tpcSignal());
          QA_d.fill(HIST("histNClusterTPC"), track.pt(), track.tpcNClsCrossedRows());
          QA_d.fill(HIST("histNClusterITS"), track.pt(), track.itsNCls());
          QA_d.fill(HIST("histNClusterITSib"), track.pt(), track.itsNClsInnerBarrel());
          QA_d.fill(HIST("histChi2TPC"), track.pt(), track.tpcChi2NCl());
          QA_d.fill(HIST("histChi2ITS"), track.pt(), track.itsChi2NCl());
          QA_d.fill(HIST("histTPCnClsFindable"), track.pt(), track.tpcNClsFindable());
          QA_d.fill(HIST("histTPCnClsFindableMinusFound"), track.pt(), track.tpcNClsFindableMinusFound());
          QA_d.fill(HIST("histTPCnClsFindableMinusCrossedRows"), track.pt(), track.tpcNClsFindableMinusCrossedRows());
          QA_d.fill(HIST("histTPCnClsShared"), track.pt(), track.tpcNClsShared());
          QA_d.fill(HIST("histChi2TOF"), track.pt(), track.tofChi2());
          QA_d.fill(HIST("histTrackLength"), track.length());
          QA_d.fill(HIST("histTPCCrossedRowsOverFindableCls"), track.tpcCrossedRowsOverFindableCls());
          QA_d.fill(HIST("histTPCFoundOverFindable"), track.tpcFoundOverFindableCls());
          QA_d.fill(HIST("histTPCFractionSharedCls"), track.tpcFractionSharedCls());

          if (track.hasTOF()) {

            if (track.hasTRD() && (lastRequiredTrdCluster > 0) && custom_Track_selection) {
              int lastLayer = 0;
              for (int l = 7; l >= 0; l--) {
                if (track.trdPattern() & (1 << l)) {
                  lastLayer = l;
                  break;
                }
              }
              if (lastLayer < lastRequiredTrdCluster) {
                continue;
              }
            }

            Float_t TOFmass2 = ((track.mass()) * (track.mass()));

            QA_d.fill(HIST("histTOFm2"), track.pt(), TOFmass2);
            QA_d.fill(HIST("histTofSignalData"), track.tpcInnerParam() * track.sign(), track.beta());
          }
        }
        if (track.sign() < 0) {
          QA_antid.fill(HIST("histTpcSignalData"), track.tpcInnerParam(), track.tpcSignal());
          QA_antid.fill(HIST("histNClusterTPC"), track.pt(), track.tpcNClsCrossedRows());
          QA_antid.fill(HIST("histNClusterITS"), track.pt(), track.itsNCls());
          QA_antid.fill(HIST("histNClusterITSib"), track.pt(), track.itsNClsInnerBarrel());
          QA_antid.fill(HIST("histChi2TPC"), track.pt(), track.tpcChi2NCl());
          QA_antid.fill(HIST("histChi2ITS"), track.pt(), track.itsChi2NCl());
          QA_antid.fill(HIST("histTPCnClsFindable"), track.pt(), track.tpcNClsFindable());
          QA_antid.fill(HIST("histTPCnClsFindableMinusFound"), track.pt(), track.tpcNClsFindableMinusFound());
          QA_antid.fill(HIST("histTPCnClsFindableMinusCrossedRows"), track.pt(), track.tpcNClsFindableMinusCrossedRows());
          QA_antid.fill(HIST("histTPCnClsShared"), track.pt(), track.tpcNClsShared());
          QA_antid.fill(HIST("histChi2TOF"), track.pt(), track.tofChi2());
          QA_antid.fill(HIST("histTrackLength"), track.length());
          QA_antid.fill(HIST("histTPCCrossedRowsOverFindableCls"), track.tpcCrossedRowsOverFindableCls());
          QA_antid.fill(HIST("histTPCFoundOverFindable"), track.tpcFoundOverFindableCls());
          QA_antid.fill(HIST("histTPCFractionSharedCls"), track.tpcFractionSharedCls());

          if (track.hasTOF()) {

            if (track.hasTRD() && (lastRequiredTrdCluster > 0) && custom_Track_selection) {
              int lastLayer = 0;
              for (int l = 7; l >= 0; l--) {
                if (track.trdPattern() & (1 << l)) {
                  lastLayer = l;
                  break;
                }
              }
              if (lastLayer < lastRequiredTrdCluster) {
                continue;
              }
            }

            Float_t TOFmass2 = ((track.mass()) * (track.mass()));

            QA_antid.fill(HIST("histTOFm2"), track.pt(), TOFmass2);
            QA_antid.fill(HIST("histTofSignalData"), track.tpcInnerParam() * track.sign(), track.beta());
          }
        }
      }


      // fill QA histograms (triton)
      if (TMath::Abs(nSigmaTriton) < nsigmacut && enable_triton_processing) {
        if (track.sign() > 0) {
          QA_t.fill(HIST("histTpcSignalData"), track.tpcInnerParam(), track.tpcSignal());
          QA_t.fill(HIST("histNClusterTPC"), track.pt(), track.tpcNClsCrossedRows());
          QA_t.fill(HIST("histNClusterITS"), track.pt(), track.itsNCls());
          QA_t.fill(HIST("histNClusterITSib"), track.pt(), track.itsNClsInnerBarrel());
          QA_t.fill(HIST("histChi2TPC"), track.pt(), track.tpcChi2NCl());
          QA_t.fill(HIST("histChi2ITS"), track.pt(), track.itsChi2NCl());
          QA_t.fill(HIST("histTPCnClsFindable"), track.pt(), track.tpcNClsFindable());
          QA_t.fill(HIST("histTPCnClsFindableMinusFound"), track.pt(), track.tpcNClsFindableMinusFound());
          QA_t.fill(HIST("histTPCnClsFindableMinusCrossedRows"), track.pt(), track.tpcNClsFindableMinusCrossedRows());
          QA_t.fill(HIST("histTPCnClsShared"), track.pt(), track.tpcNClsShared());
          QA_t.fill(HIST("histChi2TOF"), track.pt(), track.tofChi2());
          QA_t.fill(HIST("histTrackLength"), track.length());
          QA_t.fill(HIST("histTPCCrossedRowsOverFindableCls"), track.tpcCrossedRowsOverFindableCls());
          QA_t.fill(HIST("histTPCFoundOverFindable"), track.tpcFoundOverFindableCls());
          QA_t.fill(HIST("histTPCFractionSharedCls"), track.tpcFractionSharedCls());

          if (track.hasTOF()) {

            if (track.hasTRD() && (lastRequiredTrdCluster > 0) && custom_Track_selection) {
              int lastLayer = 0;
              for (int l = 7; l >= 0; l--) {
                if (track.trdPattern() & (1 << l)) {
                  lastLayer = l;
                  break;
                }
              }
              if (lastLayer < lastRequiredTrdCluster) {
                continue;
              }
            }

            Float_t TOFmass2 = ((track.mass()) * (track.mass()));

            QA_t.fill(HIST("histTOFm2"), track.pt(), TOFmass2);
            QA_t.fill(HIST("histTofSignalData"), track.tpcInnerParam() * track.sign(), track.beta());
          }
        }
        if (track.sign() < 0) {
          QA_antit.fill(HIST("histTpcSignalData"), track.tpcInnerParam(), track.tpcSignal());
          QA_antit.fill(HIST("histNClusterTPC"), track.pt(), track.tpcNClsCrossedRows());
          QA_antit.fill(HIST("histNClusterITS"), track.pt(), track.itsNCls());
          QA_antit.fill(HIST("histNClusterITSib"), track.pt(), track.itsNClsInnerBarrel());
          QA_antit.fill(HIST("histChi2TPC"), track.pt(), track.tpcChi2NCl());
          QA_antit.fill(HIST("histChi2ITS"), track.pt(), track.itsChi2NCl());
          QA_antit.fill(HIST("histTPCnClsFindable"), track.pt(), track.tpcNClsFindable());
          QA_antit.fill(HIST("histTPCnClsFindableMinusFound"), track.pt(), track.tpcNClsFindableMinusFound());
          QA_antit.fill(HIST("histTPCnClsFindableMinusCrossedRows"), track.pt(), track.tpcNClsFindableMinusCrossedRows());
          QA_antit.fill(HIST("histTPCnClsShared"), track.pt(), track.tpcNClsShared());
          QA_antit.fill(HIST("histChi2TOF"), track.pt(), track.tofChi2());
          QA_antit.fill(HIST("histTrackLength"), track.length());
          QA_antit.fill(HIST("histTPCCrossedRowsOverFindableCls"), track.tpcCrossedRowsOverFindableCls());
          QA_antit.fill(HIST("histTPCFoundOverFindable"), track.tpcFoundOverFindableCls());
          QA_antit.fill(HIST("histTPCFractionSharedCls"), track.tpcFractionSharedCls());

          if (track.hasTOF()) {

            if (track.hasTRD() && (lastRequiredTrdCluster > 0) && custom_Track_selection) {
              int lastLayer = 0;
              for (int l = 7; l >= 0; l--) {
                if (track.trdPattern() & (1 << l)) {
                  lastLayer = l;
                  break;
                }
              }
              if (lastLayer < lastRequiredTrdCluster) {
                continue;
              }
            }

            Float_t TOFmass2 = ((track.mass()) * (track.mass()));

            QA_antit.fill(HIST("histTOFm2"), track.pt(), TOFmass2);
            QA_antit.fill(HIST("histTofSignalData"), track.tpcInnerParam() * track.sign(), track.beta());
          }
        }
      }


      // fill QA histograms (Helium-3)
      if (TMath::Abs(nSigmaHe3) < nsigmacut && enable_Helium3_processing) {
        if (track.sign() > 0) {
          QA_He3.fill(HIST("histTpcSignalData"), 2 * track.tpcInnerParam(), track.tpcSignal());
          QA_He3.fill(HIST("histNClusterTPC"), 2 * track.pt(), track.tpcNClsCrossedRows());
          QA_He3.fill(HIST("histNClusterITS"), 2 * track.pt(), track.itsNCls());
          QA_He3.fill(HIST("histNClusterITSib"), 2 * track.pt(), track.itsNClsInnerBarrel());
          QA_He3.fill(HIST("histChi2TPC"), 2 * track.pt(), track.tpcChi2NCl());
          QA_He3.fill(HIST("histChi2ITS"), 2 * track.pt(), track.itsChi2NCl());
          QA_He3.fill(HIST("histTPCnClsFindable"), 2 * track.pt(), track.tpcNClsFindable());
          QA_He3.fill(HIST("histTPCnClsFindableMinusFound"), 2 * track.pt(), track.tpcNClsFindableMinusFound());
          QA_He3.fill(HIST("histTPCnClsFindableMinusCrossedRows"), 2 * track.pt(), track.tpcNClsFindableMinusCrossedRows());
          QA_He3.fill(HIST("histTPCnClsShared"), 2 * track.pt(), track.tpcNClsShared());
          QA_He3.fill(HIST("histChi2TOF"), 2 * track.pt(), track.tofChi2());
          QA_He3.fill(HIST("histTrackLength"), track.length());
          QA_He3.fill(HIST("histTPCCrossedRowsOverFindableCls"), track.tpcCrossedRowsOverFindableCls());
          QA_He3.fill(HIST("histTPCFoundOverFindable"), track.tpcFoundOverFindableCls());
          QA_He3.fill(HIST("histTPCFractionSharedCls"), track.tpcFractionSharedCls());

          if (track.hasTOF()) {

            if (track.hasTRD() && (lastRequiredTrdCluster > 0) && custom_Track_selection) {
              int lastLayer = 0;
              for (int l = 7; l >= 0; l--) {
                if (track.trdPattern() & (1 << l)) {
                  lastLayer = l;
                  break;
                }
              }
              if (lastLayer < lastRequiredTrdCluster) {
                continue;
              }
            }

            Float_t TOFmass2 = ((track.mass()) * (track.mass()));

            QA_He3.fill(HIST("histTOFm2"), 2 * track.pt(), TOFmass2);
            QA_He3.fill(HIST("histTofSignalData"), 2 * track.tpcInnerParam() * track.sign(), track.beta());
          }
        }
        if (track.sign() < 0) {
          QA_antiHe3.fill(HIST("histTpcSignalData"), 2 * track.tpcInnerParam(), track.tpcSignal());
          QA_antiHe3.fill(HIST("histNClusterTPC"), 2 * track.pt(), track.tpcNClsCrossedRows());
          QA_antiHe3.fill(HIST("histNClusterITS"), 2 * track.pt(), track.itsNCls());
          QA_antiHe3.fill(HIST("histNClusterITSib"), 2 * track.pt(), track.itsNClsInnerBarrel());
          QA_antiHe3.fill(HIST("histChi2TPC"), 2 * track.pt(), track.tpcChi2NCl());
          QA_antiHe3.fill(HIST("histChi2ITS"), 2 * track.pt(), track.itsChi2NCl());
          QA_antiHe3.fill(HIST("histTPCnClsFindable"), 2 * track.pt(), track.tpcNClsFindable());
          QA_antiHe3.fill(HIST("histTPCnClsFindableMinusFound"), 2 * track.pt(), track.tpcNClsFindableMinusFound());
          QA_antiHe3.fill(HIST("histTPCnClsFindableMinusCrossedRows"), 2 * track.pt(), track.tpcNClsFindableMinusCrossedRows());
          QA_antiHe3.fill(HIST("histTPCnClsShared"), 2 * track.pt(), track.tpcNClsShared());
          QA_antiHe3.fill(HIST("histChi2TOF"), 2 * track.pt(), track.tofChi2());
          QA_antiHe3.fill(HIST("histTrackLength"), track.length());
          QA_antiHe3.fill(HIST("histTPCCrossedRowsOverFindableCls"), track.tpcCrossedRowsOverFindableCls());
          QA_antiHe3.fill(HIST("histTPCFoundOverFindable"), track.tpcFoundOverFindableCls());
          QA_antiHe3.fill(HIST("histTPCFractionSharedCls"), track.tpcFractionSharedCls());

          if (track.hasTOF()) {

            if (track.hasTRD() && (lastRequiredTrdCluster > 0) && custom_Track_selection) {
              int lastLayer = 0;
              for (int l = 7; l >= 0; l--) {
                if (track.trdPattern() & (1 << l)) {
                  lastLayer = l;
                  break;
                }
              }
              if (lastLayer < lastRequiredTrdCluster) {
                continue;
              }
            }

            Float_t TOFmass2 = ((track.mass()) * (track.mass()));

            QA_antiHe3.fill(HIST("histTOFm2"), 2 * track.pt(), TOFmass2);
            QA_antiHe3.fill(HIST("histTofSignalData"), 2 * track.tpcInnerParam() * track.sign(), track.beta());
          }
        }
      }


      // fill QA histograms (Helium-4)
      if (TMath::Abs(nSigmaHe4) < nsigmacut && enable_Helium4_processing) {
        if (track.sign() > 0) {
          QA_He4.fill(HIST("histTpcSignalData"), 2 * track.tpcInnerParam(), track.tpcSignal());
          QA_He4.fill(HIST("histNClusterTPC"), 2 * track.pt(), track.tpcNClsCrossedRows());
          QA_He4.fill(HIST("histNClusterITS"), 2 * track.pt(), track.itsNCls());
          QA_He4.fill(HIST("histNClusterITSib"), 2 * track.pt(), track.itsNClsInnerBarrel());
          QA_He4.fill(HIST("histChi2TPC"), 2 * track.pt(), track.tpcChi2NCl());
          QA_He4.fill(HIST("histChi2ITS"), 2 * track.pt(), track.itsChi2NCl());
          QA_He4.fill(HIST("histTPCnClsFindable"), 2 * track.pt(), track.tpcNClsFindable());
          QA_He4.fill(HIST("histTPCnClsFindableMinusFound"), 2 * track.pt(), track.tpcNClsFindableMinusFound());
          QA_He4.fill(HIST("histTPCnClsFindableMinusCrossedRows"), 2 * track.pt(), track.tpcNClsFindableMinusCrossedRows());
          QA_He4.fill(HIST("histTPCnClsShared"), 2 * track.pt(), track.tpcNClsShared());
          QA_He4.fill(HIST("histChi2TOF"), 2 * track.pt(), track.tofChi2());
          QA_He4.fill(HIST("histTrackLength"), track.length());
          QA_He4.fill(HIST("histTPCCrossedRowsOverFindableCls"), track.tpcCrossedRowsOverFindableCls());
          QA_He4.fill(HIST("histTPCFoundOverFindable"), track.tpcFoundOverFindableCls());
          QA_He4.fill(HIST("histTPCFractionSharedCls"), track.tpcFractionSharedCls());

          if (track.hasTOF()) {

            if (track.hasTRD() && (lastRequiredTrdCluster > 0) && custom_Track_selection) {
              int lastLayer = 0;
              for (int l = 7; l >= 0; l--) {
                if (track.trdPattern() & (1 << l)) {
                  lastLayer = l;
                  break;
                }
              }
              if (lastLayer < lastRequiredTrdCluster) {
                continue;
              }
            }

            Float_t TOFmass2 = ((track.mass()) * (track.mass()));

            QA_He4.fill(HIST("histTOFm2"), 2 * track.pt(), TOFmass2);
            QA_He4.fill(HIST("histTofSignalData"), 2 * track.tpcInnerParam() * track.sign(), track.beta());
          }
        }
        if (track.sign() < 0) {
          QA_antiHe4.fill(HIST("histTpcSignalData"), 2 * track.tpcInnerParam(), track.tpcSignal());
          QA_antiHe4.fill(HIST("histNClusterTPC"), 2 * track.pt(), track.tpcNClsCrossedRows());
          QA_antiHe4.fill(HIST("histNClusterITS"), 2 * track.pt(), track.itsNCls());
          QA_antiHe4.fill(HIST("histNClusterITSib"), 2 * track.pt(), track.itsNClsInnerBarrel());
          QA_antiHe4.fill(HIST("histChi2TPC"), 2 * track.pt(), track.tpcChi2NCl());
          QA_antiHe4.fill(HIST("histChi2ITS"), 2 * track.pt(), track.itsChi2NCl());
          QA_antiHe4.fill(HIST("histTPCnClsFindable"), 2 * track.pt(), track.tpcNClsFindable());
          QA_antiHe4.fill(HIST("histTPCnClsFindableMinusFound"), 2 * track.pt(), track.tpcNClsFindableMinusFound());
          QA_antiHe4.fill(HIST("histTPCnClsFindableMinusCrossedRows"), 2 * track.pt(), track.tpcNClsFindableMinusCrossedRows());
          QA_antiHe4.fill(HIST("histTPCnClsShared"), 2 * track.pt(), track.tpcNClsShared());
          QA_antiHe4.fill(HIST("histChi2TOF"), 2 * track.pt(), track.tofChi2());
          QA_antiHe4.fill(HIST("histTrackLength"), track.length());
          QA_antiHe4.fill(HIST("histTPCCrossedRowsOverFindableCls"), track.tpcCrossedRowsOverFindableCls());
          QA_antiHe4.fill(HIST("histTPCFoundOverFindable"), track.tpcFoundOverFindableCls());
          QA_antiHe4.fill(HIST("histTPCFractionSharedCls"), track.tpcFractionSharedCls());

          if (track.hasTOF()) {

            if (track.hasTRD() && (lastRequiredTrdCluster > 0) && custom_Track_selection) {
              int lastLayer = 0;
              for (int l = 7; l >= 0; l--) {
                if (track.trdPattern() & (1 << l)) {
                  lastLayer = l;
                  break;
                }
              }
              if (lastLayer < lastRequiredTrdCluster) {
                continue;
              }
            }

            Float_t TOFmass2 = ((track.mass()) * (track.mass()));

            QA_antiHe4.fill(HIST("histTOFm2"), 2 * track.pt(), TOFmass2);
            QA_antiHe4.fill(HIST("histTofSignalData"), 2 * track.tpcInnerParam() * track.sign(), track.beta());
          }
        }
      }



    } // end loop over tracks
  }




  //****************************************************************************************************

  Filter collisionFilter = (nabs(aod::collision::posZ) < cfgCutVertex);
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta && requireGlobalTrackWoDCAInFilter());

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>;

  using EventCandidatesCent = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>>;

  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCLfFullPr, aod::pidTOFFullPr, aod::pidTPCLfFullDe, aod::pidTOFFullDe, aod::pidTPCLfFullTr, aod::pidTOFFullTr, aod::pidTPCLfFullHe, aod::pidTOFFullHe, aod::pidTPCLfFullAl, aod::pidTOFFullAl, aod::TrackSelection, aod::TrackSelectionExtension, aod::TOFSignal, aod::pidTOFmass, aod::pidTOFbeta>>;

  void processData(EventCandidatesCent::iterator const& event, TrackCandidates const& tracks)
  {
    fillHistograms(event, tracks);
  }
  PROCESS_SWITCH(QAHistTask, processData, "process data", true);

};

//****************************************************************************************************

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<QAHistTask>(cfgc, TaskName{"qa-hist"})};
}
