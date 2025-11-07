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

#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGLF/DataModel/LFParticleIdentification.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "TPDGCode.h"
#include <TH1F.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>

#include <cmath>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct QAHistTask {

  // Data
  HistogramRegistry QA_reg{"data_all_species", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry QA_species{"data_species", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry particle_reg{"data_positive", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry aparticle_reg{"data_negative", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry MC_recon_reg{"MC_particles_reco", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  OutputObj<TH1I> histPDG_reco{TH1I("PDG reconstructed", "PDG;PDG code", 18, 0.0, 18)};

  void init(o2::framework::InitContext&)
  {

    if ((do_pion == true && (do_kaon == true || do_proton == true || do_deuteron == true || do_triton == true || do_He3 == true || do_He4 == true)) || (do_kaon == true && (do_proton == true || do_deuteron == true || do_triton == true || do_He3 == true || do_He4 == true)) || (do_proton == true && (do_deuteron == true || do_triton == true || do_He3 == true || do_He4 == true)) || (do_deuteron == true && (do_triton == true || do_He3 == true || do_He4 == true)) || (do_triton == true && (do_He3 == true || do_He4 == true)) || (do_He3 == true && do_He4 == true)) {
      LOG(fatal) << "++++++++ Can't enable more than one species at a time, use subwagons for that purpose. ++++++++";
    }

    std::string species;

    if (do_pion)
      species = "pi";
    if (do_kaon)
      species = "ka";
    if (do_proton)
      species = "p";
    if (do_deuteron)
      species = "d";
    if (do_triton)
      species = "t";
    if (do_He3)
      species = "He3";
    if (do_He4)
      species = "He4";

    std::vector<double> ptBinning = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 3.6, 4., 5., 6., 8., 10., 12., 14.};
    std::vector<double> ptBinning_short = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2};
    std::vector<double> ptBinning_diff = {-14.0, -12.0, -10.0, -8.0, -6.0, -5.0, -4.0, -3.6, -3.2, -2.8, -2.4, -2.2, -2.0, -1.8, -1.6, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 3.6, 4., 5., 6., 8., 10., 12., 14.};
    std::vector<double> centBinning = {0., 1., 5., 10., 20., 30., 40., 50., 70., 100.};
    std::vector<double> etaBinning = {-1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    std::vector<double> PDGBinning = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0};

    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec ptAxis_short = {ptBinning_short, "Global #it{p}_{T} (GeV/#it{c})"};
    AxisSpec ptAxisITS_TPC = {ptBinning_short, "ITS-TPC #it{p}_{T} (GeV/#it{c})"};
    AxisSpec centAxis = {centBinning, "V0M (%)"};
    AxisSpec centralityAxis = {100, 0.0, 100.0, "VT0C (%)"};
    AxisSpec centralityAxis_extended = {105, 0.0, 105.0, "VT0C (%)"};
    AxisSpec etaAxis = {etaBinning, "#eta"};
    AxisSpec PDGBINNING = {PDGBinning, "PDG code"};

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
    particle_reg.add("histTpcSignalData", Form("Specific energy loss (%s)", species.c_str()), HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}});
    particle_reg.add("histTofSignalData", Form("TOF signal (%s)", species.c_str()), HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    particle_reg.add("histDcaVsPtData", Form("dcaXY vs Pt (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    particle_reg.add("histDcaZVsPtData", Form("dcaZ vs Pt (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    particle_reg.add("histTOFm2", Form("TOF m^2 vs Pt (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {400, 0.0, 10.0, "m^2"}});
    particle_reg.add("histNClusterTPC", Form("Number of Clusters in TPC vs Pt (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {80, 0.0, 160.0, "nCluster"}});
    particle_reg.add("histNClusterITS", Form("Number of Clusters in ITS vs Pt (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "ITS nCls"}});
    particle_reg.add("histNClusterITSib", Form("Number of Clusters in ib of ITS vs Pt (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "ITS ib nCls"}});
    particle_reg.add("histTPCnClsFindable", Form("Findable TPC clusters (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {200, 0.0, 200.0, "nCluster"}});
    particle_reg.add("histTPCnClsFindableMinusFound", Form("TPC Clusters: Findable - Found (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {60, 0.0, 60.0, "nCluster"}});
    particle_reg.add("histTPCnClsFindableMinusCrossedRows", Form("TPC Clusters: Findable - crossed rows (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {60, 0.0, 60.0, "nCluster"}});
    particle_reg.add("histTPCnClsShared", Form("Number of shared TPC clusters (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {70, 0.0, 70.0, "nCluster"}});
    particle_reg.add("histTPCCrossedRowsOverFindableCls", Form("Ratio crossed rows over findable clusters (%s)", species.c_str()), HistType::kTH1F, {{100, 0., 2.0, "Crossed Rows / Findable Cls"}});
    particle_reg.add("histTPCFoundOverFindable", Form("Ratio of found over findable clusters (%s)", species.c_str()), HistType::kTH1F, {{100, 0., 2.0, "Found Cls / Findable Cls"}});
    particle_reg.add("histTPCFractionSharedCls", Form("Fraction of shared TPC clusters (%s)", species.c_str()), HistType::kTH1F, {{100, -2.0, 2.0, "Shared Cls"}});
    particle_reg.add("histChi2TPC", Form("chi^2 TPC vs Pt (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    particle_reg.add("histChi2ITS", Form("chi^2 ITS vs Pt (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {125, 0.0, 50.0, "chi^2"}});
    particle_reg.add("histChi2ITSvsITSnCls", "chi^2 ITS vs ITS nCls", HistType::kTH2F, {{125, 0.0, 50.0, "chi^2"}, {10, 0.0, 10.0, "ITS nCls"}});
    particle_reg.add("histChi2ITSvsITSibnCls", "chi^2 ITS vs ITS ib nCls", HistType::kTH2F, {{125, 0.0, 50.0, "chi^2"}, {10, 0.0, 10.0, "ITS ib nCls"}});
    particle_reg.add("histChi2ITSvsITSnClsAll", "chi^2 ITS vs ITS nCls", HistType::kTH2F, {{125, 0.0, 50.0, "chi^2"}, {10, 0.0, 10.0, "ITS nCls"}});
    particle_reg.add("histChi2TOF", Form("chi^2 TOF vs Pt (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {75, 0.0, 15.0, "chi^2"}});
    particle_reg.add("histEtaWithOverFlow", Form("Pseudorapidity 0 - 105%% centrality (%s)", species.c_str()), HistType::kTH1F, {etaAxis});
    particle_reg.add("histEta", Form("Pseudorapidity with centrality cut (%s)", species.c_str()), HistType::kTH1F, {etaAxis});
    particle_reg.add("histEta_cent", Form("Pseudorapidity vs Centrality (%s)", species.c_str()), HistType::kTH2F, {centralityAxis_extended, etaAxis});
    particle_reg.add("histTrackLength", Form("Track length (%s)", species.c_str()), HistType::kTH1F, {{350, 0., 700., "length (cm)"}});
    particle_reg.add("histpTCorralation", "TPC-glo pT vs glo pT", HistType::kTH2F, {{100, -5.0, 5.0, "#it{p}^{global} (GeV/#it{c})"}, {80, -4.0, 4.0, "#it{p}^{TPC} - #it{p}^{global} (GeV/#it{c})"}});

    // QA choosen species (negative)
    aparticle_reg.add("histTpcSignalData", Form("Specific energy loss (anti-%s)", species.c_str()), HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}});
    aparticle_reg.add("histTofSignalData", Form("TOF signal (%s)", species.c_str()), HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    aparticle_reg.add("histDcaVsPtData", Form("dcaXY vs Pt (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    aparticle_reg.add("histDcaZVsPtData", Form("dcaZ vs Pt (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    aparticle_reg.add("histTOFm2", Form("TOF m^2 vs Pt (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {400, 0.0, 10.0, "m^2"}});
    aparticle_reg.add("histNClusterTPC", Form("Number of Clusters in TPC vs Pt (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {80, 0.0, 160.0, "nCluster"}});
    aparticle_reg.add("histNClusterITS", Form("Number of Clusters in ITS vs Pt (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "ITS nCls"}});
    aparticle_reg.add("histNClusterITSib", Form("Number of Clusters in ib of ITS vs Pt (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "ITS ib nCls"}});
    aparticle_reg.add("histTPCnClsFindable", Form("Findable TPC clusters (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {200, 0.0, 200.0, "nCluster"}});
    aparticle_reg.add("histTPCnClsFindableMinusFound", Form("TPC Clusters: Findable - Found (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {60, 0.0, 60.0, "nCluster"}});
    aparticle_reg.add("histTPCnClsFindableMinusCrossedRows", Form("TPC Clusters: Findable - crossed rows (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {60, 0.0, 60.0, "nCluster"}});
    aparticle_reg.add("histTPCnClsShared", Form("Number of shared TPC clusters (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {70, 0.0, 70.0, "nCluster"}});
    aparticle_reg.add("histTPCCrossedRowsOverFindableCls", Form("Ratio crossed rows over findable clusters (%s)", species.c_str()), HistType::kTH1F, {{100, 0., 2.0, "Crossed Rows / Findable Cls"}});
    aparticle_reg.add("histTPCFoundOverFindable", Form("Ratio of found over findable clusters (%s)", species.c_str()), HistType::kTH1F, {{100, 0., 2.0, "Found Cls / Findable Cls"}});
    aparticle_reg.add("histTPCFractionSharedCls", Form("Fraction of shared TPC clusters (%s)", species.c_str()), HistType::kTH1F, {{100, -2.0, 2.0, "Shared Cls"}});
    aparticle_reg.add("histChi2TPC", Form("chi^2 TPC vs Pt (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    aparticle_reg.add("histChi2ITS", Form("chi^2 ITS vs Pt (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {125, 0.0, 50.0, "chi^2"}});
    aparticle_reg.add("histChi2ITSvsITSnCls", "chi^2 ITS vs ITS nCls", HistType::kTH2F, {{125, 0.0, 50.0, "chi^2"}, {10, 0.0, 10.0, "ITS nCls"}});
    aparticle_reg.add("histChi2ITSvsITSibnCls", "chi^2 ITS vs ITS ib nCls", HistType::kTH2F, {{125, 0.0, 50.0, "chi^2"}, {10, 0.0, 10.0, "ITS ib nCls"}});
    aparticle_reg.add("histChi2ITSvsITSnClsAll", "chi^2 ITS vs ITS nCls", HistType::kTH2F, {{125, 0.0, 50.0, "chi^2"}, {10, 0.0, 10.0, "ITS nCls"}});
    aparticle_reg.add("histChi2TOF", Form("chi^2 TOF vs Pt (%s)", species.c_str()), HistType::kTH2F, {ptAxis, {75, 0.0, 15.0, "chi^2"}});
    aparticle_reg.add("histEtaWithOverFlow", Form("Pseudorapidity 0 - 105%% centrality (%s)", species.c_str()), HistType::kTH1F, {etaAxis});
    aparticle_reg.add("histEta", Form("Pseudorapidity with centrality cut (%s)", species.c_str()), HistType::kTH1F, {etaAxis});
    aparticle_reg.add("histEta_cent", Form("Pseudorapidity vs Centrality (%s)", species.c_str()), HistType::kTH2F, {centralityAxis_extended, etaAxis});
    aparticle_reg.add("histTrackLength", Form("Track length (%s)", species.c_str()), HistType::kTH1F, {{350, 0., 700., "length (cm)"}});
    aparticle_reg.add("histpTCorralation", "TPC-glo pT vs glo pT", HistType::kTH2F, {{100, -5.0, 5.0, "#it{p}^{global} (GeV/#it{c})"}, {80, -4.0, 4.0, "#it{p}^{TPC} - #it{p}^{global} (GeV/#it{c})"}});

    // +++++++++++++++++++++ MC ++++++++++++++++++++++++++

    // MC reconstructed
    MC_recon_reg.add("histRecVtxMC", "MC reconstructed vertex z position", HistType::kTH1F, {{400, -40., +40., "z position (cm)"}});
    MC_recon_reg.add("histCentrality", "Centrality", HistType::kTH1F, {centralityAxis_extended});
    MC_recon_reg.add("histPhi", "#phi", HistType::kTH2F, {{100, 0., 2. * TMath::Pi()}, PDGBINNING});
    MC_recon_reg.add("histEta", "#eta", HistType::kTH2F, {{102, -2.01, 2.01}, PDGBINNING});
    MC_recon_reg.add("histPt", "p_{t}", HistType::kTH2F, {ptAxis, PDGBINNING});
    MC_recon_reg.add("histDCA", "DCA xy", HistType::kTH3F, {ptAxis, {250, -0.5, 0.5, "dca"}, PDGBINNING});
    MC_recon_reg.add("histDCAz", "DCA z", HistType::kTH3F, {ptAxis, {1000, -2.0, 2.0, "dca"}, PDGBINNING});
    MC_recon_reg.add("histTpcSignalData", "Specific energy loss", HistType::kTH3F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}, PDGBINNING});
    MC_recon_reg.add("histTpcSignalData_all_species", "Specific energy loss", HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}});
    MC_recon_reg.add("histTofSignalData", "TOF signal", HistType::kTH3F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}, PDGBINNING});
    MC_recon_reg.add("histTofSignalData_all_species", "TOF signal", HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
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
    MC_recon_reg.add("histTPCpDist", "TPC p distribution", HistType::kTH2F, {{500, -20, +20, "#it{p} (GeV/c)"}, PDGBINNING});
    MC_recon_reg.add("histTPCFractionSharedCls", "Fraction of shared TPC clusters", HistType::kTH2F, {{100, -2.0, 2.0, "Shared Cls"}, PDGBINNING});
    MC_recon_reg.add("histpTCorralation", "TPC-glo p vs glo p", HistType::kTH2F, {{100, -5.0, 5.0, "#it{p}^{global} (GeV/#it{c})"}, {80, -4.0, 4.0, "#it{p}^{TPC} - #it{p}^{global} (GeV/#it{c})"}});
    MC_recon_reg.add("histpTCorralation_PDG", "TPC-glo p vs glo p", HistType::kTH3F, {{100, -5.0, 5.0, "#it{p}^{global} (GeV/#it{c})"}, {80, -4.0, 4.0, "#it{p}^{TPC} - #it{p}^{global} (GeV/#it{c})"}, PDGBINNING});
  }

  // Configurables
  Configurable<bool> do_pion{"do_pion", false, "0: disabled, 1: enabled"};
  Configurable<bool> do_kaon{"do_kaon", false, "0: disabled, 1: enabled"};
  Configurable<bool> do_proton{"do_proton", false, "0: disabled, 1: enabled"};
  Configurable<bool> do_deuteron{"do_deuteron", false, "0: disabled, 1: enabled"};
  Configurable<bool> do_triton{"do_triton", false, "0: disabled, 1: enabled"};
  Configurable<bool> do_He3{"do_He3", false, "0: disabled, 1: enabled"};
  Configurable<bool> do_He4{"do_He4", false, "0: disabled, 1: enabled"};

  Configurable<bool> event_selection_sel8{"event_selection_sel8", true, "0: disabled, 1: enabled"};
  Configurable<bool> event_selection_MC_sel8{"event_selection_MC_sel8", true, "Enable sel8 event selection in MC processing"};
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
  Configurable<bool> removeITSROFrameBorder{"removeITSROFrameBorder", false, "Remove TF border"};
  Configurable<bool> removeNoSameBunchPileup{"removeNoSameBunchPileup", false, "Remove TF border"};
  Configurable<bool> requireIsGoodZvtxFT0vsPV{"requireIsGoodZvtxFT0vsPV", false, "Remove TF border"};
  Configurable<bool> requireIsVertexITSTPC{"requireIsVertexITSTPC", false, "Remove TF border"};
  Configurable<bool> removeNoTimeFrameBorder{"removeNoTimeFrameBorder", false, "Remove TF border"};

  //***********************************************************************************

  template <typename CollisionType>
  bool isEventSelected(CollisionType const& collision)
  {
    if (removeITSROFrameBorder && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))
      return false;
    if (removeNoSameBunchPileup && !collision.selection_bit(aod::evsel::kNoSameBunchPileup))
      return false;
    if (requireIsGoodZvtxFT0vsPV && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))
      return false;
    if (requireIsVertexITSTPC && !collision.selection_bit(aod::evsel::kIsVertexITSTPC))
      return false;
    if (removeNoTimeFrameBorder && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder))
      return false;
    return true;
  }

  //***********************************************************************************

  template <typename CollisionType, typename TracksType>
  void fillDataHistograms(const CollisionType& event, const TracksType& tracks, const int Partilce_type)
  {

    if (event_selection_sel8 && event.sel8()) {
      QA_reg.fill(HIST("histRecVtxZData"), event.posZ());
    }

    if (!event_selection_sel8) {
      QA_reg.fill(HIST("histRecVtxZData"), event.posZ());
    }

    if (!isEventSelected(event))
      return;

    for (auto track : tracks) { // start loop over all tracks

      if (event_selection_sel8 && !event.sel8())
        continue;

      QA_reg.fill(HIST("histDcaVsPID"), track.dcaXY(), track.pidForTracking());
      QA_reg.fill(HIST("histDcaZVsPID"), track.dcaZ(), track.pidForTracking());
      QA_reg.fill(HIST("histpTCorralation"), track.sign() * track.pt(), track.tpcInnerParam() - track.pt());

      float TPCnSigma_particle = -100;

      float momentum;
      TLorentzVector lorentzVector{};

      switch (Partilce_type) {
        case 0: // pi plus/minus
          lorentzVector.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassPiPlus);
          TPCnSigma_particle = track.tpcNSigmaPi();
          momentum = track.pt();
          break;
        case 1: // (anti)kaon
          lorentzVector.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassKPlus);
          TPCnSigma_particle = track.tpcNSigmaKa();
          momentum = track.pt();
          break;
        case 2: // (anti)proton
          lorentzVector.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassProton);
          TPCnSigma_particle = track.tpcNSigmaPr();
          momentum = track.pt();
          break;
        case 3: // (anti)deuteron
          lorentzVector.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassDeuteron);
          TPCnSigma_particle = track.tpcNSigmaDe();
          momentum = track.pt();
          break;
        case 4: // (anti)triton
          lorentzVector.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassTriton);
          TPCnSigma_particle = track.tpcNSigmaTr();
          momentum = track.pt();
          break;
        case 5: // (anti)Helium-3
          lorentzVector.SetPtEtaPhiM(track.pt() * 2.0, track.eta(), track.phi(), constants::physics::MassHelium3);
          TPCnSigma_particle = track.tpcNSigmaHe();
          momentum = track.pt() * 2.0;
          break;
        case 6: // (anti)Helium-4
          lorentzVector.SetPtEtaPhiM(track.pt() * 2.0, track.eta(), track.phi(), constants::physics::MassAlpha);
          TPCnSigma_particle = track.tpcNSigmaAl();
          momentum = track.pt() * 2.0;
          break;
        default:
          continue;
          break;
      }

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
        QA_reg.fill(HIST("histDcaVsPtData_particle"), momentum, track.dcaXY());
        QA_reg.fill(HIST("histDcaZVsPtData_particle"), momentum, track.dcaZ());
      }

      if (track.sign() < 0) {
        QA_reg.fill(HIST("histDcaVsPtData_antiparticle"), momentum, track.dcaXY());
        QA_reg.fill(HIST("histDcaZVsPtData_antiparticle"), momentum, track.dcaZ());
      }

      if (custom_Track_selection && (TMath::Abs(track.dcaXY()) > maxDCA_XY || TMath::Abs(track.dcaZ()) > maxDCA_Z || TPCnumberClsFound < minTPCnClsFound || TPC_nCls_Crossed_Rows < minNCrossedRowsTPC || RatioCrossedRowsOverFindableTPC < minRatioCrossedRowsTPC || RatioCrossedRowsOverFindableTPC > maxRatioCrossedRowsTPC || Chi2perClusterTPC > maxChi2TPC || Chi2perClusterITS > maxChi2ITS || !(track.passedTPCRefit()) || !(track.passedITSRefit()) || (track.itsNClsInnerBarrel()) < minReqClusterITSib || (track.itsNCls()) < minReqClusterITS || track.length() < minTrackLength || track.length() > maxTrackLength || track.tpcNClsFindable() < minTPCNClsFindable || track.tpcNClsShared() > maxTPCNClsShared || track.tpcFoundOverFindableCls() < minTPCFoundOverFindable || track.tpcFoundOverFindableCls() > maxTPCFoundOverFindable)) {
        if (track.hasTOF() && track.tofChi2() > maxChi2TOF)
          continue;
        continue;
      }

      if (enable_PVcontributor && !(track.isPVContributor())) {
        continue;
      }

      if (lorentzVector.Rapidity() < yMin || lorentzVector.Rapidity() > yMax) {
        continue;
      }

      if (momentum < pTmin || momentum > pTmax)
        continue;

      QA_reg.fill(HIST("histTpcSignalData"), momentum * track.sign(), track.tpcSignal());
      QA_reg.fill(HIST("histNClusterTPC"), momentum, track.tpcNClsCrossedRows());
      QA_reg.fill(HIST("histNClusterITS"), momentum, track.itsNCls());
      QA_reg.fill(HIST("histNClusterITSib"), momentum, track.itsNClsInnerBarrel());
      QA_reg.fill(HIST("histChi2TPC"), momentum, track.tpcChi2NCl());
      QA_reg.fill(HIST("histChi2ITS"), momentum, track.itsChi2NCl());
      QA_reg.fill(HIST("histChi2ITSvsITSnCls"), track.itsChi2NCl(), track.itsNCls());
      QA_reg.fill(HIST("histChi2ITSvsITSibnCls"), track.itsChi2NCl(), track.itsNClsInnerBarrel());
      QA_reg.fill(HIST("histChi2ITSvsITSnClsAll"), track.itsChi2NCl(), track.itsNCls());
      QA_reg.fill(HIST("histChi2ITSvsITSnClsAll"), track.itsChi2NCl(), track.itsNClsInnerBarrel());
      QA_reg.fill(HIST("histTPCnClsFindable"), momentum, track.tpcNClsFindable());
      QA_reg.fill(HIST("histTPCnClsFindableMinusFound"), momentum, track.tpcNClsFindableMinusFound());
      QA_reg.fill(HIST("histTPCnClsFindableMinusCrossedRows"), momentum, track.tpcNClsFindableMinusCrossedRows());
      QA_reg.fill(HIST("histTPCnClsShared"), momentum, track.tpcNClsShared());
      QA_reg.fill(HIST("histChi2TOF"), momentum, track.tofChi2());
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

        QA_reg.fill(HIST("histTOFm2"), momentum, TOFmass2);
        QA_reg.fill(HIST("histTofSignalData"), momentum * track.sign(), track.beta());
      }

      if (TMath::Abs(TPCnSigma_particle) < nsigmacut) {

        QA_species.fill(HIST("histpTCorralation"), track.sign() * momentum, track.tpcInnerParam() - momentum);

        if (track.sign() > 0) {
          particle_reg.fill(HIST("histDcaVsPtData"), momentum, track.dcaXY());
          particle_reg.fill(HIST("histDcaZVsPtData"), momentum, track.dcaZ());
          particle_reg.fill(HIST("histTpcSignalData"), momentum, track.tpcSignal());
          particle_reg.fill(HIST("histNClusterTPC"), momentum, track.tpcNClsCrossedRows());
          particle_reg.fill(HIST("histNClusterITS"), momentum, track.itsNCls());
          particle_reg.fill(HIST("histNClusterITSib"), momentum, track.itsNClsInnerBarrel());
          particle_reg.fill(HIST("histChi2TPC"), momentum, track.tpcChi2NCl());
          particle_reg.fill(HIST("histChi2ITS"), momentum, track.itsChi2NCl());
          particle_reg.fill(HIST("histChi2ITSvsITSnCls"), track.itsChi2NCl(), track.itsNCls());
          particle_reg.fill(HIST("histChi2ITSvsITSibnCls"), track.itsChi2NCl(), track.itsNClsInnerBarrel());
          particle_reg.fill(HIST("histChi2ITSvsITSnClsAll"), track.itsChi2NCl(), track.itsNCls());
          particle_reg.fill(HIST("histChi2ITSvsITSnClsAll"), track.itsChi2NCl(), track.itsNClsInnerBarrel());
          particle_reg.fill(HIST("histTPCnClsFindable"), momentum, track.tpcNClsFindable());
          particle_reg.fill(HIST("histTPCnClsFindableMinusFound"), momentum, track.tpcNClsFindableMinusFound());
          particle_reg.fill(HIST("histTPCnClsFindableMinusCrossedRows"), momentum, track.tpcNClsFindableMinusCrossedRows());
          particle_reg.fill(HIST("histTPCnClsShared"), momentum, track.tpcNClsShared());
          particle_reg.fill(HIST("histChi2TOF"), momentum, track.tofChi2());
          particle_reg.fill(HIST("histTrackLength"), track.length());
          particle_reg.fill(HIST("histTPCCrossedRowsOverFindableCls"), track.tpcCrossedRowsOverFindableCls());
          particle_reg.fill(HIST("histTPCFoundOverFindable"), track.tpcFoundOverFindableCls());
          particle_reg.fill(HIST("histTPCFractionSharedCls"), track.tpcFractionSharedCls());
          particle_reg.fill(HIST("histpTCorralation"), track.sign() * momentum, track.tpcInnerParam() - momentum);

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

            particle_reg.fill(HIST("histTOFm2"), momentum, TOFmass2);
            particle_reg.fill(HIST("histTofSignalData"), momentum * track.sign(), track.beta());
          }
        }
        if (track.sign() < 0) {
          aparticle_reg.fill(HIST("histDcaVsPtData"), momentum, track.dcaXY());
          aparticle_reg.fill(HIST("histDcaZVsPtData"), momentum, track.dcaZ());
          aparticle_reg.fill(HIST("histTpcSignalData"), momentum, track.tpcSignal());
          aparticle_reg.fill(HIST("histNClusterTPC"), momentum, track.tpcNClsCrossedRows());
          aparticle_reg.fill(HIST("histNClusterITS"), momentum, track.itsNCls());
          aparticle_reg.fill(HIST("histNClusterITSib"), momentum, track.itsNClsInnerBarrel());
          aparticle_reg.fill(HIST("histChi2TPC"), momentum, track.tpcChi2NCl());
          aparticle_reg.fill(HIST("histChi2ITS"), momentum, track.itsChi2NCl());
          aparticle_reg.fill(HIST("histChi2ITSvsITSnCls"), track.itsChi2NCl(), track.itsNCls());
          aparticle_reg.fill(HIST("histChi2ITSvsITSibnCls"), track.itsChi2NCl(), track.itsNClsInnerBarrel());
          aparticle_reg.fill(HIST("histChi2ITSvsITSnClsAll"), track.itsChi2NCl(), track.itsNCls());
          aparticle_reg.fill(HIST("histChi2ITSvsITSnClsAll"), track.itsChi2NCl(), track.itsNClsInnerBarrel());
          aparticle_reg.fill(HIST("histTPCnClsFindable"), momentum, track.tpcNClsFindable());
          aparticle_reg.fill(HIST("histTPCnClsFindableMinusFound"), momentum, track.tpcNClsFindableMinusFound());
          aparticle_reg.fill(HIST("histTPCnClsFindableMinusCrossedRows"), momentum, track.tpcNClsFindableMinusCrossedRows());
          aparticle_reg.fill(HIST("histTPCnClsShared"), momentum, track.tpcNClsShared());
          aparticle_reg.fill(HIST("histChi2TOF"), momentum, track.tofChi2());
          aparticle_reg.fill(HIST("histTrackLength"), track.length());
          aparticle_reg.fill(HIST("histTPCCrossedRowsOverFindableCls"), track.tpcCrossedRowsOverFindableCls());
          aparticle_reg.fill(HIST("histTPCFoundOverFindable"), track.tpcFoundOverFindableCls());
          aparticle_reg.fill(HIST("histTPCFractionSharedCls"), track.tpcFractionSharedCls());
          aparticle_reg.fill(HIST("histpTCorralation"), track.sign() * momentum, track.tpcInnerParam() - momentum);

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

            aparticle_reg.fill(HIST("histTOFm2"), momentum, TOFmass2);
            aparticle_reg.fill(HIST("histTofSignalData"), momentum * track.sign(), track.beta());
          }
        }
      }
    }
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

    if (!isEventSelected(event))
      return;

    for (auto track : tracks) {

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

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>;

  using EventCandidatesCent = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>>;

  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCLfFullPi, aod::pidTOFFullPi, aod::pidTPCLfFullKa, aod::pidTOFFullKa, aod::pidTPCLfFullPr, aod::pidTOFFullPr, aod::pidTPCLfFullDe, aod::pidTOFFullDe, aod::pidTPCLfFullTr, aod::pidTOFFullTr, aod::pidTPCLfFullHe, aod::pidTOFFullHe, aod::pidTPCLfFullAl, aod::pidTOFFullAl, aod::TrackSelection, aod::TrackSelectionExtension, aod::TOFSignal, aod::pidTOFmass, aod::pidTOFbeta>>;

  void processData(EventCandidates::iterator const& event, TrackCandidates const& tracks)
  {
    if (do_pion)
      fillDataHistograms(event, tracks, 0); // pion
    if (do_kaon)
      fillDataHistograms(event, tracks, 1); // kaon
    if (do_proton)
      fillDataHistograms(event, tracks, 2); // proton
    if (do_deuteron)
      fillDataHistograms(event, tracks, 3); // deuteron
    if (do_triton)
      fillDataHistograms(event, tracks, 4); // triton
    if (do_He3)
      fillDataHistograms(event, tracks, 5); // He3
    if (do_He4)
      fillDataHistograms(event, tracks, 6); // He4
  }
  PROCESS_SWITCH(QAHistTask, processData, "process data", true);

  void processDataCent(EventCandidatesCent::iterator const& event, TrackCandidates const& tracks)
  {
    if (do_pion)
      fillDataHistograms(event, tracks, 0); // pion
    if (do_kaon)
      fillDataHistograms(event, tracks, 1); // kaon
    if (do_proton)
      fillDataHistograms(event, tracks, 2); // proton
    if (do_deuteron)
      fillDataHistograms(event, tracks, 3); // deuteron
    if (do_triton)
      fillDataHistograms(event, tracks, 4); // triton
    if (do_He3)
      fillDataHistograms(event, tracks, 5); // He3
    if (do_He4)
      fillDataHistograms(event, tracks, 6); // He4
    fillDataCentHistorgrams(event, tracks);
  }
  PROCESS_SWITCH(QAHistTask, processDataCent, "process data containing centralities", false);

  void processMCreco(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0Cs>::iterator const& collisions, soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCLfFullPi, aod::pidTOFFullPi, aod::pidTPCLfFullKa, aod::pidTOFFullKa, aod::pidTPCLfFullPr, aod::pidTOFFullPr, aod::pidTPCLfFullDe, aod::pidTOFFullDe, aod::pidTPCLfFullTr, aod::pidTOFFullTr, aod::pidTPCLfFullHe, aod::pidTOFFullHe, aod::pidTPCLfFullAl, aod::pidTOFFullAl, aod::McTrackLabels, aod::TrackSelection, aod::TrackSelectionExtension, aod::TOFSignal, aod::pidTOFmass, aod::pidTOFbeta>> const& tracks,
                     aod::McParticles& /*mcParticles*/, aod::McCollisions const& /*mcCollisions*/)
  {

    if (event_selection_MC_sel8 && !collisions.sel8())
      return;
    MC_recon_reg.fill(HIST("histRecVtxMC"), collisions.posZ());
    MC_recon_reg.fill(HIST("histCentrality"), collisions.centFT0C());
    if (!isEventSelected(collisions))
      return;

    for (auto& track : tracks) {
      const auto particle = track.mcParticle();
      if (!particle.isPhysicalPrimary())
        continue;

      MC_recon_reg.fill(HIST("histChi2ITSvsITSnCls"), track.itsChi2NCl(), track.itsNCls());
      MC_recon_reg.fill(HIST("histChi2ITSvsITSibnCls"), track.itsChi2NCl(), track.itsNClsInnerBarrel());
      MC_recon_reg.fill(HIST("histChi2ITSvsITSnClsAll"), track.itsChi2NCl(), track.itsNCls());
      MC_recon_reg.fill(HIST("histChi2ITSvsITSnClsAll"), track.itsChi2NCl(), track.itsNClsInnerBarrel());
      MC_recon_reg.fill(HIST("histpTCorralation"), track.sign() * track.p(), track.tpcInnerParam() - track.p());

      if ((particle.pdgCode() == 1000020030) || (particle.pdgCode() == -1000020030) || (particle.pdgCode() == 1000020040) || (particle.pdgCode() == -1000020040)) {
        MC_recon_reg.fill(HIST("histTpcSignalData_all_species"), track.pt() * 2.0 * track.sign(), track.tpcSignal());
      } else {
        MC_recon_reg.fill(HIST("histTpcSignalData_all_species"), track.pt() * track.sign(), track.tpcSignal());
      }
      if (track.hasTOF()) {
        if ((particle.pdgCode() == 1000020030) || (particle.pdgCode() == -1000020030) || (particle.pdgCode() == 1000020040) || (particle.pdgCode() == -1000020040)) {
          MC_recon_reg.fill(HIST("histTofSignalData_all_species"), track.pt() * 2.0 * track.sign(), track.beta());
        } else {
          MC_recon_reg.fill(HIST("histTofSignalData_all_species"), track.pt() * track.sign(), track.beta());
        }
      }

      int pdgbin = -10;
      TLorentzVector lorentzVector_particle_MC{};
      switch (particle.pdgCode()) {
        case +211:
          histPDG_reco->AddBinContent(1);
          pdgbin = 0;
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassPiPlus);
          break;
        case -211:
          histPDG_reco->AddBinContent(2);
          pdgbin = 1;
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassPiPlus);
          break;
        case +321:
          histPDG_reco->AddBinContent(3);
          pdgbin = 2;
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassKPlus);
          break;
        case -321:
          histPDG_reco->AddBinContent(4);
          pdgbin = 3;
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassKPlus);
          break;
        case +2212:
          histPDG_reco->AddBinContent(5);
          pdgbin = 4;
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassProton);
          break;
        case -2212:
          histPDG_reco->AddBinContent(6);
          pdgbin = 5;
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassProton);
          break;
        case +1000010020:
          histPDG_reco->AddBinContent(7);
          pdgbin = 6;
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassDeuteron);
          break;
        case -1000010020:
          histPDG_reco->AddBinContent(8);
          pdgbin = 7;
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassDeuteron);
          break;
        case +1000010030:
          histPDG_reco->AddBinContent(9);
          pdgbin = 8;
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassTriton);
          break;
        case -1000010030:
          histPDG_reco->AddBinContent(10);
          pdgbin = 9;
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassTriton);
          break;
        case +1000020030:
          histPDG_reco->AddBinContent(11);
          pdgbin = 10;
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt() * 2.0, track.eta(), track.phi(), constants::physics::MassHelium3);
          break;
        case -1000020030:
          histPDG_reco->AddBinContent(12);
          pdgbin = 11;
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt() * 2.0, track.eta(), track.phi(), constants::physics::MassHelium3);
          break;
        case +1000020040:
          histPDG_reco->AddBinContent(13);
          pdgbin = 12;
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt() * 2.0, track.eta(), track.phi(), constants::physics::MassAlpha);
          break;
        case -1000020040:
          histPDG_reco->AddBinContent(14);
          pdgbin = 13;
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt() * 2.0, track.eta(), track.phi(), constants::physics::MassAlpha);
          break;
        default:
          pdgbin = -10;
          continue;
          break;
      }

      MC_recon_reg.fill(HIST("histPhi"), track.phi(), pdgbin);
      MC_recon_reg.fill(HIST("histEta"), track.eta(), pdgbin);
      MC_recon_reg.fill(HIST("histTPCCrossedRowsOverFindableCls"), track.tpcCrossedRowsOverFindableCls(), pdgbin);
      MC_recon_reg.fill(HIST("histTPCFoundOverFindable"), track.tpcFoundOverFindableCls(), pdgbin);
      MC_recon_reg.fill(HIST("histTrackLength"), track.length(), pdgbin);
      MC_recon_reg.fill(HIST("histTPCFractionSharedCls"), track.tpcFractionSharedCls(), pdgbin);
      MC_recon_reg.fill(HIST("histpTCorralation_PDG"), track.sign() * track.p(), track.tpcInnerParam() - track.p(), pdgbin);
      MC_recon_reg.fill(HIST("histGlobalpDist"), track.p(), pdgbin);
      MC_recon_reg.fill(HIST("histTPCpDist"), track.tpcInnerParam(), pdgbin);

      if ((particle.pdgCode() == 1000020030) || (particle.pdgCode() == -1000020030) || (particle.pdgCode() == 1000020040) || (particle.pdgCode() == -1000020040)) {
        MC_recon_reg.fill(HIST("histPt"), track.pt() * 2.0, pdgbin);
        MC_recon_reg.fill(HIST("histDCA"), track.pt() * 2.0, track.dcaXY(), pdgbin);
        MC_recon_reg.fill(HIST("histDCAz"), track.pt() * 2.0, track.dcaZ(), pdgbin);
        MC_recon_reg.fill(HIST("histTpcSignalData"), track.pt() * 2.0 * track.sign(), track.tpcSignal(), pdgbin);
        MC_recon_reg.fill(HIST("histNClusterTPC"), track.pt() * 2.0, track.tpcNClsCrossedRows(), pdgbin);
        MC_recon_reg.fill(HIST("histNClusterITS"), track.pt() * 2.0, track.itsNCls(), pdgbin);
        MC_recon_reg.fill(HIST("histNClusterITSib"), track.pt() * 2.0, track.itsNClsInnerBarrel(), pdgbin);
        MC_recon_reg.fill(HIST("histTPCnClsFindable"), track.pt() * 2.0, track.tpcNClsFindable(), pdgbin);
        MC_recon_reg.fill(HIST("histTPCnClsFindableMinusFound"), track.pt() * 2.0, track.tpcNClsFindableMinusFound(), pdgbin);
        MC_recon_reg.fill(HIST("histTPCnClsFindableMinusCrossedRows"), track.pt() * 2.0, track.tpcNClsFindableMinusCrossedRows(), pdgbin);
        MC_recon_reg.fill(HIST("histTPCnClsShared"), track.pt() * 2.0, track.tpcNClsShared(), pdgbin);
        MC_recon_reg.fill(HIST("histChi2TPC"), track.pt() * 2.0, track.tpcChi2NCl(), pdgbin);
        MC_recon_reg.fill(HIST("histChi2ITS"), track.pt() * 2.0, track.itsChi2NCl(), pdgbin);
        MC_recon_reg.fill(HIST("histChi2TOF"), track.pt() * 2.0, track.tofChi2(), pdgbin);
      } else {
        MC_recon_reg.fill(HIST("histPt"), track.pt(), pdgbin);
        MC_recon_reg.fill(HIST("histDCA"), track.pt(), track.dcaXY(), pdgbin);
        MC_recon_reg.fill(HIST("histDCAz"), track.pt(), track.dcaZ(), pdgbin);
        MC_recon_reg.fill(HIST("histTpcSignalData"), track.pt() * track.sign(), track.tpcSignal(), pdgbin);
        MC_recon_reg.fill(HIST("histNClusterTPC"), track.pt(), track.tpcNClsCrossedRows(), pdgbin);
        MC_recon_reg.fill(HIST("histNClusterITS"), track.pt(), track.itsNCls(), pdgbin);
        MC_recon_reg.fill(HIST("histNClusterITSib"), track.pt(), track.itsNClsInnerBarrel(), pdgbin);
        MC_recon_reg.fill(HIST("histTPCnClsFindable"), track.pt(), track.tpcNClsFindable(), pdgbin);
        MC_recon_reg.fill(HIST("histTPCnClsFindableMinusFound"), track.pt(), track.tpcNClsFindableMinusFound(), pdgbin);
        MC_recon_reg.fill(HIST("histTPCnClsFindableMinusCrossedRows"), track.pt(), track.tpcNClsFindableMinusCrossedRows(), pdgbin);
        MC_recon_reg.fill(HIST("histTPCnClsShared"), track.pt(), track.tpcNClsShared(), pdgbin);
        MC_recon_reg.fill(HIST("histChi2TPC"), track.pt(), track.tpcChi2NCl(), pdgbin);
        MC_recon_reg.fill(HIST("histChi2ITS"), track.pt(), track.itsChi2NCl(), pdgbin);
        MC_recon_reg.fill(HIST("histChi2TOF"), track.pt(), track.tofChi2(), pdgbin);
      }
      if (track.hasTOF()) {
        Float_t TOFmass2 = ((track.mass()) * (track.mass()));
        if ((particle.pdgCode() == 1000020030) || (particle.pdgCode() == -1000020030) || (particle.pdgCode() == 1000020040) || (particle.pdgCode() == -1000020040)) {
          MC_recon_reg.fill(HIST("histTOFm2"), track.pt() * 2.0, TOFmass2, pdgbin);
          MC_recon_reg.fill(HIST("histTofSignalData"), track.pt() * 2.0 * track.sign(), track.beta(), pdgbin);
        } else {
          MC_recon_reg.fill(HIST("histTOFm2"), track.pt(), TOFmass2, pdgbin);
          MC_recon_reg.fill(HIST("histTofSignalData"), track.pt() * track.sign(), track.beta(), pdgbin);
        }
      }
    }
  }
  PROCESS_SWITCH(QAHistTask, processMCreco, "process MC", false);
};

//****************************************************************************************************

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<QAHistTask>(cfgc, TaskName{"qa-hist"})};
}
