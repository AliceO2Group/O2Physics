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
#include <TH1F.h>
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct QAHistTask {

  // Data
  HistogramRegistry QA_reg{"data_all_species", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry QA_species{"data_species", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry QA_species_pos{"data_positive", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry QA_species_neg{"data_negative", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry MC_truth_reg{"MC_particles", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry MC_recon_reg{"MC_particles_reconstructed", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry MC_recon_diff_reg{"MC_reconstructed_diff", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  OutputObj<TH1F> histPDG{TH1F("PDG", "PDG;PDG code", 100, 0.0, 100.0)};

  void init(o2::framework::InitContext&)
  {

    if ((process_pion == true && (process_kaon == true || process_proton == true || process_deuteron == true || process_triton == true || process_He3 == true || process_He4 == true)) || (process_kaon == true && (process_proton == true || process_deuteron == true || process_triton == true || process_He3 == true || process_He4 == true)) || (process_proton == true && (process_deuteron == true || process_triton == true || process_He3 == true || process_He4 == true)) || (process_deuteron == true && (process_triton == true || process_He3 == true || process_He4 == true)) || (process_triton == true && (process_He3 == true || process_He4 == true)) || (process_He3 == true && process_He4 == true)) {
      LOG(fatal) << "++++++++ Can't enable more than one species at a time, use subwagons for that purpose. ++++++++";
    }

    std::string species;

    if (process_pion)
      species = "pi";
    if (process_kaon)
      species = "ka";
    if (process_proton)
      species = "p";
    if (process_deuteron)
      species = "d";
    if (process_triton)
      species = "t";
    if (process_He3)
      species = "He3";
    if (process_He4)
      species = "He4";

    std::vector<double> ptBinning = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 3.6, 4., 5., 6., 8., 10., 12., 14.};
    std::vector<double> ptBinning_short = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2};
    std::vector<double> ptBinning_diff = {-14.0, -12.0, -10.0, -8.0, -6.0, -5.0, -4.0, -3.6, -3.2, -2.8, -2.4, -2.2, -2.0, -1.8, -1.6, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 3.6, 4., 5., 6., 8., 10., 12., 14.};
    std::vector<double> centBinning = {0., 1., 5., 10., 20., 30., 40., 50., 70., 100.};
    std::vector<double> etaBinning = {-1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec ptAxis_short = {ptBinning_short, "Global #it{p}_{T} (GeV/#it{c})"};
    AxisSpec ptAxisITS_TPC = {ptBinning_short, "ITS-TPC #it{p}_{T} (GeV/#it{c})"};
    AxisSpec centAxis = {centBinning, "V0M (%)"};
    AxisSpec centralityAxis = {100, 0.0, 100.0, "VT0C (%)"};
    AxisSpec centralityAxis_extended = {105, 0.0, 105.0, "VT0C (%)"};
    AxisSpec etaAxis = {etaBinning, "#eta"};
    AxisSpec PDGBINNING = {20, 0.0, 20.0, "PDG code"};
    AxisSpec ptAxis_diff = {ptBinning_diff, "#it{p}_{T,diff} (GeV/#it{c})"};

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
    QA_reg.add("histChi2ITSvsITSnCls", "chi^2 ITS vs ITS nCls", HistType::kTH2F, {{125, 0.0, 50.0, "chi^2"}, {10, 0.0, 10.0, "ITS nCls"}});
    QA_reg.add("histChi2ITSvsITSibnCls", "chi^2 ITS vs ITS ib nCls", HistType::kTH2F, {{125, 0.0, 50.0, "chi^2"}, {10, 0.0, 10.0, "ITS ib nCls"}});
    QA_reg.add("histChi2ITSvsITSnClsAll", "chi^2 ITS vs ITS nCls", HistType::kTH2F, {{125, 0.0, 50.0, "chi^2"}, {10, 0.0, 10.0, "ITS nCls"}});
    QA_reg.add("histChi2TOF", "chi^2 TOF vs Pt", HistType::kTH2F, {ptAxis, {75, 0.0, 15.0, "chi^2"}});
    QA_reg.add("histEtaWithOverFlow", "Pseudorapidity 0 - 105%% centrality", HistType::kTH1F, {etaAxis});
    QA_reg.add("histEta", "Pseudorapidity with centrality cut", HistType::kTH1F, {etaAxis});
    QA_reg.add("histEta_cent", "Pseudorapidity vs Centrality", HistType::kTH2F, {centralityAxis_extended, etaAxis});
    QA_reg.add("histTrackLength", "Track length", HistType::kTH1F, {{350, 0., 700., "length (cm)"}});
    QA_reg.add("histDcaVsPID", "DCA XY vs PID hypothesis", HistType::kTH2F, {{1000, -2.0, 2.0, "dca XY"}, {10, 0.0, 10.0, "PID ID"}});
    QA_reg.add("histDcaZVsPID", "DCA Z vs PID hypothesis", HistType::kTH2F, {{1000, -2.0, 2.0, "dca Z"}, {10, 0.0, 10.0, "PID ID"}});
    QA_reg.add("histpTCorralation", "TPC-glo pT vs glo pT", HistType::kTH2F, {{100, -5.0, 5.0, "#it{p}^{global} (GeV/#it{c})"}, {80, -4.0, 4.0, "#it{p}^{TPC} - #it{p}^{global} (GeV/#it{c})"}});

    // QA choosen species (all charges)
    QA_species.add("histpTCorralation", "TPC-glo pT vs glo pT", HistType::kTH2F, {{100, -5.0, 5.0, "#it{p}^{global} (GeV/#it{c})"}, {80, -4.0, 4.0, "#it{p}^{TPC} - #it{p}^{global} (GeV/#it{c})"}});

    // QA choosen species (positive)
    QA_species_pos.add("histTpcSignalData", Form("Specific energy loss (%s)", species.c_str()), HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}});
    QA_species_pos.add("histTofSignalData", Form("TOF signal (%s)", species.c_str()), HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    QA_species_pos.add("histDcaVsPtData", Form("dcaXY vs Pt (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    QA_species_pos.add("histDcaZVsPtData", Form("dcaZ vs Pt (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    QA_species_pos.add("histTOFm2", Form("TOF m^2 vs Pt (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {400, 0.0, 10.0, "m^2"}});
    QA_species_pos.add("histNClusterTPC", Form("Number of Clusters in TPC vs Pt (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {80, 0.0, 160.0, "nCluster"}});
    QA_species_pos.add("histNClusterITS", Form("Number of Clusters in ITS vs Pt (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "ITS nCls"}});
    QA_species_pos.add("histNClusterITSib", Form("Number of Clusters in ib of ITS vs Pt (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "ITS ib nCls"}});
    QA_species_pos.add("histTPCnClsFindable", Form("Findable TPC clusters (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {200, 0.0, 200.0, "nCluster"}});
    QA_species_pos.add("histTPCnClsFindableMinusFound", Form("TPC Clusters: Findable - Found (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {60, 0.0, 60.0, "nCluster"}});
    QA_species_pos.add("histTPCnClsFindableMinusCrossedRows", Form("TPC Clusters: Findable - crossed rows (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {60, 0.0, 60.0, "nCluster"}});
    QA_species_pos.add("histTPCnClsShared", Form("Number of shared TPC clusters (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {70, 0.0, 70.0, "nCluster"}});
    QA_species_pos.add("histTPCCrossedRowsOverFindableCls", Form("Ratio crossed rows over findable clusters (%s)", species.c_str()), HistType::kTH1F, {{100, 0., 2.0, "Crossed Rows / Findable Cls"}});
    QA_species_pos.add("histTPCFoundOverFindable", Form("Ratio of found over findable clusters (%s)", species.c_str()), HistType::kTH1F, {{100, 0., 2.0, "Found Cls / Findable Cls"}});
    QA_species_pos.add("histTPCFractionSharedCls", Form("Fraction of shared TPC clusters (%s)", species.c_str()), HistType::kTH1F, {{100, -2.0, 2.0, "Shared Cls"}});
    QA_species_pos.add("histChi2TPC", Form("chi^2 TPC vs Pt (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    QA_species_pos.add("histChi2ITS", Form("chi^2 ITS vs Pt (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {125, 0.0, 50.0, "chi^2"}});
    QA_species_pos.add("histChi2ITSvsITSnCls", "chi^2 ITS vs ITS nCls", HistType::kTH2F, {{125, 0.0, 50.0, "chi^2"}, {10, 0.0, 10.0, "ITS nCls"}});
    QA_species_pos.add("histChi2ITSvsITSibnCls", "chi^2 ITS vs ITS ib nCls", HistType::kTH2F, {{125, 0.0, 50.0, "chi^2"}, {10, 0.0, 10.0, "ITS ib nCls"}});
    QA_species_pos.add("histChi2ITSvsITSnClsAll", "chi^2 ITS vs ITS nCls", HistType::kTH2F, {{125, 0.0, 50.0, "chi^2"}, {10, 0.0, 10.0, "ITS nCls"}});
    QA_species_pos.add("histChi2TOF", Form("chi^2 TOF vs Pt (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {75, 0.0, 15.0, "chi^2"}});
    QA_species_pos.add("histEtaWithOverFlow", Form("Pseudorapidity 0 - 105%% centrality (%s)", species.c_str()), HistType::kTH1F, {etaAxis});
    QA_species_pos.add("histEta", Form("Pseudorapidity with centrality cut (%s)", species.c_str()), HistType::kTH1F, {etaAxis});
    QA_species_pos.add("histEta_cent", Form("Pseudorapidity vs Centrality (%s)", species.c_str()), HistType::kTH2F, {centralityAxis_extended, etaAxis});
    QA_species_pos.add("histTrackLength", Form("Track length (%s)", species.c_str()), HistType::kTH1F, {{350, 0., 700., "length (cm)"}});
    QA_species_pos.add("histpTCorralation", "TPC-glo pT vs glo pT", HistType::kTH2F, {{100, -5.0, 5.0, "#it{p}^{global} (GeV/#it{c})"}, {80, -4.0, 4.0, "#it{p}^{TPC} - #it{p}^{global} (GeV/#it{c})"}});

    // QA choosen species (negative)
    QA_species_neg.add("histTpcSignalData", Form("Specific energy loss (anti-%s)", species.c_str()), HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}});
    QA_species_neg.add("histTofSignalData", Form("TOF signal (%s)", species.c_str()), HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    QA_species_neg.add("histDcaVsPtData", Form("dcaXY vs Pt (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    QA_species_neg.add("histDcaZVsPtData", Form("dcaZ vs Pt (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    QA_species_neg.add("histTOFm2", Form("TOF m^2 vs Pt (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {400, 0.0, 10.0, "m^2"}});
    QA_species_neg.add("histNClusterTPC", Form("Number of Clusters in TPC vs Pt (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {80, 0.0, 160.0, "nCluster"}});
    QA_species_neg.add("histNClusterITS", Form("Number of Clusters in ITS vs Pt (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "ITS nCls"}});
    QA_species_neg.add("histNClusterITSib", Form("Number of Clusters in ib of ITS vs Pt (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "ITS ib nCls"}});
    QA_species_neg.add("histTPCnClsFindable", Form("Findable TPC clusters (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {200, 0.0, 200.0, "nCluster"}});
    QA_species_neg.add("histTPCnClsFindableMinusFound", Form("TPC Clusters: Findable - Found (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {60, 0.0, 60.0, "nCluster"}});
    QA_species_neg.add("histTPCnClsFindableMinusCrossedRows", Form("TPC Clusters: Findable - crossed rows (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {60, 0.0, 60.0, "nCluster"}});
    QA_species_neg.add("histTPCnClsShared", Form("Number of shared TPC clusters (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {70, 0.0, 70.0, "nCluster"}});
    QA_species_neg.add("histTPCCrossedRowsOverFindableCls", Form("Ratio crossed rows over findable clusters (%s)", species.c_str()), HistType::kTH1F, {{100, 0., 2.0, "Crossed Rows / Findable Cls"}});
    QA_species_neg.add("histTPCFoundOverFindable", Form("Ratio of found over findable clusters (%s)", species.c_str()), HistType::kTH1F, {{100, 0., 2.0, "Found Cls / Findable Cls"}});
    QA_species_neg.add("histTPCFractionSharedCls", Form("Fraction of shared TPC clusters (%s)", species.c_str()), HistType::kTH1F, {{100, -2.0, 2.0, "Shared Cls"}});
    QA_species_neg.add("histChi2TPC", Form("chi^2 TPC vs Pt (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    QA_species_neg.add("histChi2ITS", Form("chi^2 ITS vs Pt (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {125, 0.0, 50.0, "chi^2"}});
    QA_species_neg.add("histChi2ITSvsITSnCls", "chi^2 ITS vs ITS nCls", HistType::kTH2F, {{125, 0.0, 50.0, "chi^2"}, {10, 0.0, 10.0, "ITS nCls"}});
    QA_species_neg.add("histChi2ITSvsITSibnCls", "chi^2 ITS vs ITS ib nCls", HistType::kTH2F, {{125, 0.0, 50.0, "chi^2"}, {10, 0.0, 10.0, "ITS ib nCls"}});
    QA_species_neg.add("histChi2ITSvsITSnClsAll", "chi^2 ITS vs ITS nCls", HistType::kTH2F, {{125, 0.0, 50.0, "chi^2"}, {10, 0.0, 10.0, "ITS nCls"}});
    QA_species_neg.add("histChi2TOF", Form("chi^2 TOF vs Pt (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {75, 0.0, 15.0, "chi^2"}});
    QA_species_neg.add("histEtaWithOverFlow", Form("Pseudorapidity 0 - 105%% centrality (%s)", species.c_str()), HistType::kTH1F, {etaAxis});
    QA_species_neg.add("histEta", Form("Pseudorapidity with centrality cut (%s)", species.c_str()), HistType::kTH1F, {etaAxis});
    QA_species_neg.add("histEta_cent", Form("Pseudorapidity vs Centrality (%s)", species.c_str()), HistType::kTH2F, {centralityAxis_extended, etaAxis});
    QA_species_neg.add("histTrackLength", Form("Track length (%s)", species.c_str()), HistType::kTH1F, {{350, 0., 700., "length (cm)"}});
    QA_species_neg.add("histpTCorralation", "TPC-glo pT vs glo pT", HistType::kTH2F, {{100, -5.0, 5.0, "#it{p}^{global} (GeV/#it{c})"}, {80, -4.0, 4.0, "#it{p}^{TPC} - #it{p}^{global} (GeV/#it{c})"}});

    // MC truth
    MC_truth_reg.add("histPhi", "#phi", HistType::kTH2F, {{100, 0., 2. * TMath::Pi()}, PDGBINNING});
    MC_truth_reg.add("histEta", "#eta", HistType::kTH2F, {{102, -2.01, 2.01}, PDGBINNING});
    MC_truth_reg.add("histPt", "p_{t}", HistType::kTH2F, {ptAxis, PDGBINNING});
    MC_truth_reg.add("histRecVtxMC", "MC collision z position", HistType::kTH1F, {{400, -40., +40., "z position (cm)"}});
    MC_truth_reg.add("histCentrality", "Centrality", HistType::kTH1F, {centralityAxis_extended});

    // MC reconstructed
    MC_recon_reg.add("histPhi", "#phi", HistType::kTH2F, {{100, 0., 2. * TMath::Pi()}, PDGBINNING});
    MC_recon_reg.add("histEta", "#eta", HistType::kTH2F, {{102, -2.01, 2.01}, PDGBINNING});
    MC_recon_reg.add("histPt", "p_{t}", HistType::kTH2F, {ptAxis, PDGBINNING});
    MC_recon_reg.add("histDCA", "DCA xy", HistType::kTH3F, {ptAxis, {250, -0.5, 0.5, "dca"}, PDGBINNING});
    MC_recon_reg.add("histDCAz", "DCA z", HistType::kTH3F, {ptAxis, {1000, -2.0, 2.0, "dca"}, PDGBINNING});
    MC_recon_reg.add("histTpcSignalData", "Specific energy loss", HistType::kTH3F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}, PDGBINNING});
    MC_recon_reg.add("histTofSignalData", "TOF signal", HistType::kTH3F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}, PDGBINNING});
    MC_recon_reg.add("histTOFm2", "TOF m^2 vs Pt", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2"}, PDGBINNING});
    MC_recon_reg.add("histNClusterTPC", "Number of Clusters in TPC vs Pt", HistType::kTH3F, {ptAxis, {80, 0.0, 160.0, "nCluster"}, PDGBINNING});
    MC_recon_reg.add("histNClusterITS", "Number of Clusters in ITS vs Pt", HistType::kTH3F, {ptAxis, {10, 0.0, 10.0, "ITS nCls"}, PDGBINNING});
    MC_recon_reg.add("histNClusterITSib", "Number of Clusters in ib of ITS vs Pt", HistType::kTH3F, {ptAxis, {10, 0.0, 10.0, "ITS ib nCls"}, PDGBINNING});
    MC_recon_reg.add("histTPCnClsFindable", "Findable TPC clusters", HistType::kTH3F, {ptAxis, {200, 0.0, 200.0, "nCluster"}, PDGBINNING});
    MC_recon_reg.add("histTPCnClsFindableMinusFound", "TPC Clusters: Findable - Found", HistType::kTH3F, {ptAxis, {60, 0.0, 60.0, "nCluster"}, PDGBINNING});
    MC_recon_reg.add("histTPCnClsFindableMinusCrossedRows", "TPC Clusters: Findable - crossed rows", HistType::kTH3F, {ptAxis, {60, 0.0, 60.0, "nCluster"}, PDGBINNING});
    MC_recon_reg.add("histTPCnClsShared", "Number of shared TPC clusters", HistType::kTH3F, {ptAxis, {70, 0.0, 70.0, "nCluster"}, PDGBINNING});
    MC_recon_reg.add("histTPCCrossedRowsOverFindableCls", "Ratio crossed rows over findable clusters", HistType::kTH2F, {{100, 0., 2.0, "Crossed Rows / Findable Cls"}, PDGBINNING});
    MC_recon_reg.add("histTPCFoundOverFindable", "Ratio of found over findable clusters", HistType::kTH2F, {{100, 0., 2.0, "Found Cls / Findable Cls"}, PDGBINNING});
    MC_recon_reg.add("histChi2TPC", "chi^2 TPC vs Pt", HistType::kTH3F, {ptAxis, {100, 0.0, 5.0, "chi^2"}, PDGBINNING});
    MC_recon_reg.add("histChi2ITS", "chi^2 ITS vs Pt", HistType::kTH3F, {ptAxis, {500, 0.0, 50.0, "chi^2"}, PDGBINNING});
    MC_recon_reg.add("histChi2ITSvsITSnCls", "chi^2 ITS vs ITS nCls", HistType::kTH2F, {{125, 0.0, 50.0, "chi^2"}, {10, 0.0, 10.0, "ITS nCls"}});
    MC_recon_reg.add("histChi2ITSvsITSibnCls", "chi^2 ITS vs ITS ib nCls", HistType::kTH2F, {{125, 0.0, 50.0, "chi^2"}, {10, 0.0, 10.0, "ITS ib nCls"}});
    MC_recon_reg.add("histChi2ITSvsITSnClsAll", "chi^2 ITS vs ITS nCls", HistType::kTH2F, {{125, 0.0, 50.0, "chi^2"}, {10, 0.0, 10.0, "ITS nCls"}});
    MC_recon_reg.add("histChi2TOF", "chi^2 TOF vs Pt", HistType::kTH3F, {ptAxis, {75, 0.0, 15.0, "chi^2"}, PDGBINNING});
    MC_recon_reg.add("histTrackLength", "Track length", HistType::kTH2F, {{350, 0., 700., "length (cm)"}, PDGBINNING});
    MC_recon_reg.add("histGlobalpDist", "Global p distribution", HistType::kTH2F, {{150, -10, +10, "#it{p} (GeV/c)"}, PDGBINNING});
    MC_recon_reg.add("histTPCpDist", "TPC p distribution", HistType::kTH2F, {{150, -10, +10, "#it{p} (GeV/c)"}, PDGBINNING});
    MC_recon_reg.add("histTPCFractionSharedCls", "Fraction of shared TPC clusters", HistType::kTH2F, {{100, -2.0, 2.0, "Shared Cls"}, PDGBINNING});
    MC_recon_reg.add("histpTCorralation", "TPC-glo p vs glo p", HistType::kTH2F, {{100, -5.0, 5.0, "#it{p}^{global} (GeV/#it{c})"}, {80, -4.0, 4.0, "#it{p}^{TPC} - #it{p}^{global} (GeV/#it{c})"}});
    MC_recon_reg.add("histpTCorralation_PDG", "TPC-glo p vs glo p", HistType::kTH3F, {{100, -5.0, 5.0, "#it{p}^{global} (GeV/#it{c})"}, {80, -4.0, 4.0, "#it{p}^{TPC} - #it{p}^{global} (GeV/#it{c})"}, PDGBINNING});

    MC_recon_reg.add("histpTCorralation_pion", "TPC-glo p vs glo p", HistType::kTH2F, {{100, -5.0, 5.0, "#it{p}^{global} (GeV/#it{c})"}, {80, -4.0, 4.0, "#it{p}^{TPC} - #it{p}^{global} (GeV/#it{c})"}});
    MC_recon_reg.add("histpTCorralation_apion", "TPC-glo p vs glo p", HistType::kTH2F, {{100, -5.0, 5.0, "#it{p}^{global} (GeV/#it{c})"}, {80, -4.0, 4.0, "#it{p}^{TPC} - #it{p}^{global} (GeV/#it{c})"}});
    MC_recon_reg.add("histpTCorralation_2", "TPC p vs global p", HistType::kTH3F, {{110, -0.5, 5.0, "#it{p}^{global} (GeV/#it{c})"}, {110, -0.5, 5.0, "#it{p}^{TPC} (GeV/#it{c})"}, PDGBINNING});

    // MC diff (truth - reconstructed)
    MC_recon_diff_reg.add("histPhiDiff", "MC t", HistType::kTH2F, {ptAxis_diff, PDGBINNING});
    MC_recon_diff_reg.add("histEtaDiff", "MC t", HistType::kTH2F, {ptAxis_diff, PDGBINNING});
    MC_recon_diff_reg.add("histPtDiff", "MC t", HistType::kTH2F, {ptAxis_diff, PDGBINNING});
  }

  // Configurables
  Configurable<bool> process_pion{"process_pion", false, "0: disabled, 1: enabled"};
  Configurable<bool> process_kaon{"process_kaon", false, "0: disabled, 1: enabled"};
  Configurable<bool> process_proton{"process_proton", false, "0: disabled, 1: enabled"};
  Configurable<bool> process_deuteron{"process_deuteron", false, "0: disabled, 1: enabled"};
  Configurable<bool> process_triton{"process_triton", false, "0: disabled, 1: enabled"};
  Configurable<bool> process_He3{"process_He3", false, "0: disabled, 1: enabled"};
  Configurable<bool> process_He4{"process_He4", false, "0: disabled, 1: enabled"};

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
  void fillDataHistograms(const CollisionType& event, const TracksType& tracks)
  {

    if (event_selection_sel8 && event.sel8()) {
      QA_reg.fill(HIST("histRecVtxZData"), event.posZ());
    }

    if (!event_selection_sel8) {
      QA_reg.fill(HIST("histRecVtxZData"), event.posZ());
    }

    for (auto track : tracks) { // start loop over tracks

      QA_reg.fill(HIST("histDcaVsPID"), track.dcaXY(), track.pidForTracking());
      QA_reg.fill(HIST("histDcaZVsPID"), track.dcaZ(), track.pidForTracking());
      QA_reg.fill(HIST("histpTCorralation"), track.sign() * track.pt(), track.tpcInnerParam() - track.pt());

      float nSigmaSpecies = 999.0;

      if (process_pion)
        nSigmaSpecies = track.tpcNSigmaPi();
      if (process_kaon)
        nSigmaSpecies = track.tpcNSigmaKa();
      if (process_proton)
        nSigmaSpecies = track.tpcNSigmaPr();
      if (process_deuteron)
        nSigmaSpecies = track.tpcNSigmaDe();
      if (process_triton)
        nSigmaSpecies = track.tpcNSigmaTr();
      if (process_He3)
        nSigmaSpecies = track.tpcNSigmaHe();
      if (process_He4)
        nSigmaSpecies = track.tpcNSigmaAl();

      if (event_selection_sel8 && !event.sel8()) {
        continue;
      }

      QA_reg.fill(HIST("histEtaWithOverFlow"), track.eta());

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

        if (track.hasTOF() && track.tofChi2() > maxChi2TOF)
          continue;
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

      if (track.pt() < pTmin || track.pt() > pTmax)
        continue;

      QA_reg.fill(HIST("histTpcSignalData"), track.pt() * track.sign(), track.tpcSignal());
      QA_reg.fill(HIST("histNClusterTPC"), track.pt(), track.tpcNClsCrossedRows());
      QA_reg.fill(HIST("histNClusterITS"), track.pt(), track.itsNCls());
      QA_reg.fill(HIST("histNClusterITSib"), track.pt(), track.itsNClsInnerBarrel());
      QA_reg.fill(HIST("histChi2TPC"), track.pt(), track.tpcChi2NCl());
      QA_reg.fill(HIST("histChi2ITS"), track.pt(), track.itsChi2NCl());
      QA_reg.fill(HIST("histChi2ITSvsITSnCls"), track.itsChi2NCl(), track.itsNCls());
      QA_reg.fill(HIST("histChi2ITSvsITSibnCls"), track.itsChi2NCl(), track.itsNClsInnerBarrel());
      QA_reg.fill(HIST("histChi2ITSvsITSnClsAll"), track.itsChi2NCl(), track.itsNCls());
      QA_reg.fill(HIST("histChi2ITSvsITSnClsAll"), track.itsChi2NCl(), track.itsNClsInnerBarrel());
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
        QA_reg.fill(HIST("histTofSignalData"), track.pt() * track.sign(), track.beta());
      }

      // fill QA histograms
      if (TMath::Abs(nSigmaSpecies) < nsigmacut) {

        QA_species.fill(HIST("histpTCorralation"), track.sign() * track.pt(), track.tpcInnerParam() - track.pt());

        if (track.sign() > 0) {
          QA_species_pos.fill(HIST("histDcaVsPtData"), track.pt(), track.dcaXY());
          QA_species_pos.fill(HIST("histDcaZVsPtData"), track.pt(), track.dcaZ());
          QA_species_pos.fill(HIST("histTpcSignalData"), track.pt(), track.tpcSignal());
          QA_species_pos.fill(HIST("histNClusterTPC"), track.pt(), track.tpcNClsCrossedRows());
          QA_species_pos.fill(HIST("histNClusterITS"), track.pt(), track.itsNCls());
          QA_species_pos.fill(HIST("histNClusterITSib"), track.pt(), track.itsNClsInnerBarrel());
          QA_species_pos.fill(HIST("histChi2TPC"), track.pt(), track.tpcChi2NCl());
          QA_species_pos.fill(HIST("histChi2ITS"), track.pt(), track.itsChi2NCl());
          QA_species_pos.fill(HIST("histChi2ITSvsITSnCls"), track.itsChi2NCl(), track.itsNCls());
          QA_species_pos.fill(HIST("histChi2ITSvsITSibnCls"), track.itsChi2NCl(), track.itsNClsInnerBarrel());
          QA_species_pos.fill(HIST("histChi2ITSvsITSnClsAll"), track.itsChi2NCl(), track.itsNCls());
          QA_species_pos.fill(HIST("histChi2ITSvsITSnClsAll"), track.itsChi2NCl(), track.itsNClsInnerBarrel());
          QA_species_pos.fill(HIST("histTPCnClsFindable"), track.pt(), track.tpcNClsFindable());
          QA_species_pos.fill(HIST("histTPCnClsFindableMinusFound"), track.pt(), track.tpcNClsFindableMinusFound());
          QA_species_pos.fill(HIST("histTPCnClsFindableMinusCrossedRows"), track.pt(), track.tpcNClsFindableMinusCrossedRows());
          QA_species_pos.fill(HIST("histTPCnClsShared"), track.pt(), track.tpcNClsShared());
          QA_species_pos.fill(HIST("histChi2TOF"), track.pt(), track.tofChi2());
          QA_species_pos.fill(HIST("histTrackLength"), track.length());
          QA_species_pos.fill(HIST("histTPCCrossedRowsOverFindableCls"), track.tpcCrossedRowsOverFindableCls());
          QA_species_pos.fill(HIST("histTPCFoundOverFindable"), track.tpcFoundOverFindableCls());
          QA_species_pos.fill(HIST("histTPCFractionSharedCls"), track.tpcFractionSharedCls());
          QA_species_pos.fill(HIST("histpTCorralation"), track.sign() * track.pt(), track.tpcInnerParam() - track.pt());

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

            QA_species_pos.fill(HIST("histTOFm2"), track.pt(), TOFmass2);
            QA_species_pos.fill(HIST("histTofSignalData"), track.pt() * track.sign(), track.beta());
          }
        }
        if (track.sign() < 0) {
          QA_species_neg.fill(HIST("histDcaVsPtData"), track.pt(), track.dcaXY());
          QA_species_neg.fill(HIST("histDcaZVsPtData"), track.pt(), track.dcaZ());
          QA_species_neg.fill(HIST("histTpcSignalData"), track.pt(), track.tpcSignal());
          QA_species_neg.fill(HIST("histNClusterTPC"), track.pt(), track.tpcNClsCrossedRows());
          QA_species_neg.fill(HIST("histNClusterITS"), track.pt(), track.itsNCls());
          QA_species_neg.fill(HIST("histNClusterITSib"), track.pt(), track.itsNClsInnerBarrel());
          QA_species_neg.fill(HIST("histChi2TPC"), track.pt(), track.tpcChi2NCl());
          QA_species_neg.fill(HIST("histChi2ITS"), track.pt(), track.itsChi2NCl());
          QA_species_neg.fill(HIST("histChi2ITSvsITSnCls"), track.itsChi2NCl(), track.itsNCls());
          QA_species_neg.fill(HIST("histChi2ITSvsITSibnCls"), track.itsChi2NCl(), track.itsNClsInnerBarrel());
          QA_species_neg.fill(HIST("histChi2ITSvsITSnClsAll"), track.itsChi2NCl(), track.itsNCls());
          QA_species_neg.fill(HIST("histChi2ITSvsITSnClsAll"), track.itsChi2NCl(), track.itsNClsInnerBarrel());
          QA_species_neg.fill(HIST("histTPCnClsFindable"), track.pt(), track.tpcNClsFindable());
          QA_species_neg.fill(HIST("histTPCnClsFindableMinusFound"), track.pt(), track.tpcNClsFindableMinusFound());
          QA_species_neg.fill(HIST("histTPCnClsFindableMinusCrossedRows"), track.pt(), track.tpcNClsFindableMinusCrossedRows());
          QA_species_neg.fill(HIST("histTPCnClsShared"), track.pt(), track.tpcNClsShared());
          QA_species_neg.fill(HIST("histChi2TOF"), track.pt(), track.tofChi2());
          QA_species_neg.fill(HIST("histTrackLength"), track.length());
          QA_species_neg.fill(HIST("histTPCCrossedRowsOverFindableCls"), track.tpcCrossedRowsOverFindableCls());
          QA_species_neg.fill(HIST("histTPCFoundOverFindable"), track.tpcFoundOverFindableCls());
          QA_species_neg.fill(HIST("histTPCFractionSharedCls"), track.tpcFractionSharedCls());
          QA_species_neg.fill(HIST("histpTCorralation"), track.sign() * track.pt(), track.tpcInnerParam() - track.pt());

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

            QA_species_neg.fill(HIST("histTOFm2"), track.pt(), TOFmass2);
            QA_species_neg.fill(HIST("histTofSignalData"), track.pt() * track.sign(), track.beta());
          }
        }
      }
    } // end loop over tracks
  }

  //*******************************************************************************************************

  template <typename CollisionType, typename TracksType>
  void fillDataCentHistorgrams(const CollisionType& event, const TracksType& tracks)
  {

    if (event_selection_sel8 && event.sel8()) {
      QA_reg.fill(HIST("histCentrality"), event.centFT0C());
    }

    if (!event_selection_sel8) {
      QA_reg.fill(HIST("histCentrality"), event.centFT0C());
    }

    for (auto track : tracks) { // start loop over tracks

      if (event_selection_sel8 && !event.sel8()) {
        continue;
      }

      QA_reg.fill(HIST("histEta_cent"), event.centFT0C(), track.eta());

      if ((event.centFT0C() > minCentrality) && (event.centFT0C() < maxCentrality) && enable_Centrality_cut) {
        QA_reg.fill(HIST("histEta"), track.eta());
      }
    }
  }

  //*******************************************************************************************************

  Filter collisionFilter = (nabs(aod::collision::posZ) < cfgCutVertex);
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta && requireGlobalTrackWoDCAInFilter());

  using EventCandidatesData = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>;

  using EventCandidatesDataCent = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>>;

  using TrackCandidatesData = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCLfFullPi, aod::pidTOFFullPi, aod::pidTPCLfFullKa, aod::pidTOFFullKa, aod::pidTPCLfFullPr, aod::pidTOFFullPr, aod::pidTPCLfFullDe, aod::pidTOFFullDe, aod::pidTPCLfFullTr, aod::pidTOFFullTr, aod::pidTPCLfFullHe, aod::pidTOFFullHe, aod::pidTPCLfFullAl, aod::pidTOFFullAl, aod::TrackSelection, aod::TrackSelectionExtension, aod::TOFSignal, aod::pidTOFmass, aod::pidTOFbeta>>;

  void processData(EventCandidatesData::iterator const& event, TrackCandidatesData const& tracks)
  {
    fillDataHistograms(event, tracks);
  }
  PROCESS_SWITCH(QAHistTask, processData, "process data", true);

  void processDataCent(EventCandidatesDataCent::iterator const& event, TrackCandidatesData const& tracks)
  {
    fillDataHistograms(event, tracks);
    fillDataCentHistorgrams(event, tracks);
  }
  PROCESS_SWITCH(QAHistTask, processDataCent, "process data containing centralities", false);

  void processMC(soa::Join<aod::Collisions, aod::McCollisionLabels>::iterator const& collisions, soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels, aod::TrackSelection, aod::TrackSelectionExtension, aod::TOFSignal, aod::pidTOFmass, aod::pidTOFbeta>> const& tracks,
                 aod::McParticles& mcParticles, aod::McCollisions const& mcCollisions)
  {

    MC_truth_reg.fill(HIST("histRecVtxMC"), collisions.posZ());

    for (auto& track : tracks) {
      const auto particle = track.mcParticle();
      const auto pdg = Form("%i", particle.pdgCode());

      if (!particle.isPhysicalPrimary())
        continue;

      histPDG->Fill(pdg, 1);
      const float pdgbin = histPDG->GetXaxis()->GetBinCenter(histPDG->GetXaxis()->FindBin(pdg));

      MC_truth_reg.fill(HIST("histPhi"), particle.phi(), pdgbin);
      MC_truth_reg.fill(HIST("histEta"), particle.eta(), pdgbin);
      MC_truth_reg.fill(HIST("histPt"), particle.pt(), pdgbin);

      MC_recon_reg.fill(HIST("histPhi"), track.phi(), pdgbin);
      MC_recon_reg.fill(HIST("histEta"), track.eta(), pdgbin);
      MC_recon_reg.fill(HIST("histPt"), track.pt(), pdgbin);
      MC_recon_reg.fill(HIST("histDCA"), track.pt(), track.dcaXY(), pdgbin);
      MC_recon_reg.fill(HIST("histDCAz"), track.pt(), track.dcaZ(), pdgbin);
      MC_recon_reg.fill(HIST("histTpcSignalData"), track.pt(), track.tpcSignal(), pdgbin);
      MC_recon_reg.fill(HIST("histNClusterTPC"), track.pt(), track.tpcNClsCrossedRows(), pdgbin);
      MC_recon_reg.fill(HIST("histNClusterITS"), track.pt(), track.itsNCls(), pdgbin);
      MC_recon_reg.fill(HIST("histNClusterITSib"), track.pt(), track.itsNClsInnerBarrel(), pdgbin);
      MC_recon_reg.fill(HIST("histTPCnClsFindable"), track.pt(), track.tpcNClsFindable(), pdgbin);
      MC_recon_reg.fill(HIST("histTPCnClsFindableMinusFound"), track.pt(), track.tpcNClsFindableMinusFound(), pdgbin);
      MC_recon_reg.fill(HIST("histTPCnClsFindableMinusCrossedRows"), track.pt(), track.tpcNClsFindableMinusCrossedRows(), pdgbin);
      MC_recon_reg.fill(HIST("histTPCnClsShared"), track.pt(), track.tpcNClsShared(), pdgbin);
      MC_recon_reg.fill(HIST("histTPCCrossedRowsOverFindableCls"), track.tpcCrossedRowsOverFindableCls(), pdgbin);
      MC_recon_reg.fill(HIST("histTPCFoundOverFindable"), track.tpcFoundOverFindableCls(), pdgbin);
      MC_recon_reg.fill(HIST("histChi2TPC"), track.pt(), track.tpcChi2NCl(), pdgbin);
      MC_recon_reg.fill(HIST("histChi2ITS"), track.pt(), track.itsChi2NCl(), pdgbin);
      MC_recon_reg.fill(HIST("histChi2ITSvsITSnCls"), track.itsChi2NCl(), track.itsNCls());
      MC_recon_reg.fill(HIST("histChi2ITSvsITSibnCls"), track.itsChi2NCl(), track.itsNClsInnerBarrel());
      MC_recon_reg.fill(HIST("histChi2ITSvsITSnClsAll"), track.itsChi2NCl(), track.itsNCls());
      MC_recon_reg.fill(HIST("histChi2ITSvsITSnClsAll"), track.itsChi2NCl(), track.itsNClsInnerBarrel());
      MC_recon_reg.fill(HIST("histChi2TOF"), track.pt(), track.tofChi2(), pdgbin);
      MC_recon_reg.fill(HIST("histTrackLength"), track.length(), pdgbin);
      MC_recon_reg.fill(HIST("histTPCFractionSharedCls"), track.tpcFractionSharedCls(), pdgbin);
      MC_recon_reg.fill(HIST("histpTCorralation"), track.sign() * track.p(), track.tpcInnerParam() - track.p());
      MC_recon_reg.fill(HIST("histpTCorralation_PDG"), track.sign() * track.p(), track.tpcInnerParam() - track.p(), pdgbin);

      MC_recon_reg.fill(HIST("histGlobalpDist"), track.p(), pdgbin);
      MC_recon_reg.fill(HIST("histTPCpDist"), track.tpcInnerParam(), pdgbin);
      MC_recon_reg.fill(HIST("histpTCorralation_2"), track.p(), track.tpcInnerParam(), pdgbin);

      if (track.hasTOF()) {
        Float_t TOFmass2 = ((track.mass()) * (track.mass()));

        MC_recon_reg.fill(HIST("histTOFm2"), track.pt(), TOFmass2, pdgbin);
        MC_recon_reg.fill(HIST("histTofSignalData"), track.pt(), track.beta(), pdgbin);
      }

      MC_recon_diff_reg.fill(HIST("histEtaDiff"), particle.eta() - track.eta(), pdgbin);
      auto delta = particle.phi() - track.phi();
      if (delta > TMath::Pi()) {
        delta -= 2 * TMath::Pi();
      }
      if (delta < -TMath::Pi()) {
        delta += 2 * TMath::Pi();
      }

      MC_recon_diff_reg.fill(HIST("histPhiDiff"), delta, pdgbin);
      MC_recon_diff_reg.fill(HIST("histPtDiff"), particle.pt() - track.pt(), pdgbin);
    }
  }
  PROCESS_SWITCH(QAHistTask, processMC, "process MC", false);

  void processMCCent(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::CentFT0Cs>::iterator const& collisions)
  {

    MC_truth_reg.fill(HIST("histCentrality"), collisions.centFT0C());
  }
  PROCESS_SWITCH(QAHistTask, processMCCent, "process MC with centrality", false);
};

//****************************************************************************************************

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<QAHistTask>(cfgc, TaskName{"qa-hist"})};
}
