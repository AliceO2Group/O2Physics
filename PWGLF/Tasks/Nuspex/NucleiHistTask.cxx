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
// Date: 30.11.2022

#include <cmath>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <vector>

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

struct NucleiHistTask {

  // Data
  HistogramRegistry spectra_reg{"spectra", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry pion_reg{"pion", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry apion_reg{"apion", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry proton_reg{"proton", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry aproton_reg{"aproton", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry deuteron_reg{"deuteron", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry adeuteron_reg{"adeuteron", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry triton_reg{"triton", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry atriton_reg{"atriton", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry Helium3_reg{"Helium3", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry aHelium3_reg{"aHelium3", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry Helium4_reg{"Helium4", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry aHelium4_reg{"aHelium4", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // MC
  HistogramRegistry MC_recon_reg{"MC_particles_reco", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  OutputObj<TH1I> histPDG{TH1I("PDG", "PDG;PDG code", 20, 0.0, 20.0)};
  OutputObj<TH1I> histTrackcuts_data{TH1I("histTrackcuts_data", "Entires;Track cut", 15, 0, 15)};
  OutputObj<TH1I> histTrackcuts_MC{TH1I("histTrackcuts_MC", "Entires;Track cut", 15, 0, 15)};

  void init(o2::framework::InitContext&)
  {

    if (doprocessData == true && doprocessDataCent == true) {
      LOG(fatal) << "Can't enable processData and processDataCent in the same time, pick one!";
    }

    // +++++++++++++++++++++ Data ++++++++++++++++++++++++
    histTrackcuts_data->GetXaxis()->SetBinLabel(1, "Events read");
    histTrackcuts_data->GetXaxis()->SetBinLabel(2, "Ev. sel. passed");
    histTrackcuts_data->GetXaxis()->SetBinLabel(3, "Rap. cut passed");
    histTrackcuts_data->GetXaxis()->SetBinLabel(4, "DCA cut passed");
    histTrackcuts_data->GetXaxis()->SetBinLabel(5, "TPCnCls cut passed");
    histTrackcuts_data->GetXaxis()->SetBinLabel(6, "TPCCrossedRowsOverFindable cut passed");
    histTrackcuts_data->GetXaxis()->SetBinLabel(7, "Chi2 cut passed");
    histTrackcuts_data->GetXaxis()->SetBinLabel(8, "Passed TPC refit cut");
    histTrackcuts_data->GetXaxis()->SetBinLabel(9, "Passed ITS refit cut");
    histTrackcuts_data->GetXaxis()->SetBinLabel(10, "ITSnCls cut passed");
    histTrackcuts_data->GetXaxis()->SetBinLabel(11, "Momentum cut passed");
    histTrackcuts_data->GetXaxis()->SetBinLabel(12, "hasITS & hasTPC cut passed");
    histTrackcuts_data->GetXaxis()->SetBinLabel(13, "GoldenChi2 cut passed");

    // +++++++++++++++++++++ MC ++++++++++++++++++++++++
    histTrackcuts_MC->GetXaxis()->SetBinLabel(1, "Events read");
    histTrackcuts_MC->GetXaxis()->SetBinLabel(2, "Prim. particle. sel. passed");
    histTrackcuts_MC->GetXaxis()->SetBinLabel(3, "Rap. cut passed");
    histTrackcuts_MC->GetXaxis()->SetBinLabel(4, "DCA cut passed");
    histTrackcuts_MC->GetXaxis()->SetBinLabel(5, "TPCnCls cut passed");
    histTrackcuts_MC->GetXaxis()->SetBinLabel(6, "TPCCrossedRowsOverFindable cut passed");
    histTrackcuts_MC->GetXaxis()->SetBinLabel(7, "Chi2 cut passed");
    histTrackcuts_MC->GetXaxis()->SetBinLabel(8, "Passed TPC refit cut");
    histTrackcuts_MC->GetXaxis()->SetBinLabel(9, "Passed ITS refit cut");
    histTrackcuts_MC->GetXaxis()->SetBinLabel(10, "ITSnCls cut passed");
    histTrackcuts_MC->GetXaxis()->SetBinLabel(11, "Momentum cut passed");
    histTrackcuts_MC->GetXaxis()->SetBinLabel(12, "hasITS & hasTPC cut passed");
    histTrackcuts_MC->GetXaxis()->SetBinLabel(13, "GoldenChi2 cut passed");

    std::vector<double> ptBinning = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 3.6, 4., 5., 6., 8., 10., 12., 14.};
    std::vector<double> etaBinning = {-1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    std::vector<double> PDGBinning = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0};

    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec pAxis = {ptBinning, "#it{p} (GeV/#it{c})"};
    AxisSpec centralityAxis = {100, 0.0, 105.0, "VT0C (%)"};
    AxisSpec etaAxis = {etaBinning, "#eta"};
    AxisSpec PDGBINNING = {PDGBinning, "PDG code"};

    // +++++++++++++++++++++ Data ++++++++++++++++++++++++

    // QA histograms
    spectra_reg.add("histRecVtxZData", "collision z position", HistType::kTH1F, {{200, -20., +20., "z position (cm)"}});
    spectra_reg.add("histTpcSignalData", "Specific energy loss", HistType::kTH2F, {{600, -6., 6., "#it{p*} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}});
    spectra_reg.add("histTofSignalData", "TOF signal", HistType::kTH2F, {{600, -6., 6., "#it{p*} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    spectra_reg.add("histDcaVsPtData_particle", "dcaXY vs Pt (particle)", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    spectra_reg.add("histDcaZVsPtData_particle", "dcaZ vs Pt (particle)", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    spectra_reg.add("histDcaVsPtData_antiparticle", "dcaXY vs Pt (antiparticle)", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    spectra_reg.add("histDcaZVsPtData_antiparticle", "dcaZ vs Pt (antiparticle)", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    spectra_reg.add("histTOFm2", "TOF m^2 vs P", HistType::kTH2F, {pAxis, {400, 0.0, 10.0, "m^2"}});
    spectra_reg.add("histNClusterTPC", "Number of Clusters in TPC vs Pt", HistType::kTH2F, {ptAxis, {160, 0.0, 160.0, "nCluster"}});
    spectra_reg.add("histNClusterITS", "Number of Clusters in ITS vs Pt", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    spectra_reg.add("histNClusterITSib", "Number of Clusters in ib of ITS vs Pt", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    spectra_reg.add("histChi2TPC", "chi^2 TPC vs Pt", HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    spectra_reg.add("histChi2ITS", "chi^2 ITS vs Pt", HistType::kTH2F, {ptAxis, {500, 0.0, 50.0, "chi^2"}});
    spectra_reg.add("histCentrality", "Centrality", HistType::kTH1F, {centralityAxis});
    spectra_reg.add("histEtaWithOverFlow", "Pseudorapidity 0 - 105%% centrality", HistType::kTH1F, {etaAxis});
    spectra_reg.add("histEta", "Pseudorapidity with centrality cut", HistType::kTH1F, {etaAxis});
    spectra_reg.add("histEta_cent", "Pseudorapidity vs Centrality", HistType::kTH2F, {centralityAxis, etaAxis});

    // histograms for pi⁺
    pion_reg.add("histTpcSignalData", "Specific energy loss (#pi^{+})", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}});
    pion_reg.add("histTofSignalData", "TOF signal (#pi^{+})", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    pion_reg.add("histDcaVsPtData", "dcaXY vs Pt (#pi^{+})", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    pion_reg.add("histDcaZVsPtData", "dcaZ vs Pt (#pi^{+})", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    pion_reg.add("histTOFm2", "TOF m^2 vs Pt (#pi^{+})", HistType::kTH2F, {pAxis, {400, 0.0, 10.0, "m^2"}});
    pion_reg.add("histTpcNsigmaData", "n-sigma TPC (#pi^{+})", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{#pi^{+}}"}});
    pion_reg.add("histTofNsigmaData", "n-sigma TOF (#pi^{+})", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{#pi^{+}}"}});
    pion_reg.add("histNClusterTPC", "Number of Clusters in TPC vs Pt (#pi^{+})", HistType::kTH2F, {ptAxis, {160, 0.0, 160.0, "nCluster"}});
    pion_reg.add("histNClusterITS", "Number of Clusters in ITS vs Pt (#pi^{+})", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    pion_reg.add("histNClusterITSib", "Number of Clusters in ib of ITS vs Pt (#pi^{+})", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    pion_reg.add("histChi2TPC", "chi^2 TPC vs Pt (#pi^{+})", HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    pion_reg.add("histChi2ITS", "chi^2 ITS vs Pt (#pi^{+})", HistType::kTH2F, {ptAxis, {500, 0.0, 50.0, "chi^2"}});
    pion_reg.add("histTpcNsigmaData_cent", "n-sigma TPC (#pi^{+}) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{#pi^{+}}"}, centralityAxis});
    pion_reg.add("histTofNsigmaData_cent", "n-sigma TOF (#pi^{+}) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{#pi^{+}}"}, centralityAxis});
    pion_reg.add("histTofm2_cent", "mass^2 TOF (#pi^{+}) centrality", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{#pi^{+}}"}, centralityAxis});
    pion_reg.add("histTpcNsigmaData_eta", "n-sigma TPC (#pi^{+}) vs eta", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{#pi^{+}}"}, etaAxis});
    pion_reg.add("histTofNsigmaData_eta", "n-sigma TOF (#pi^{+}) vs eta", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{#pi^{+}}"}, etaAxis});
    pion_reg.add("histTofm2_eta", "mass^2 TOF (#pi^{+}) vs eta", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{#pi^{+}}"}, etaAxis});

    // histograms for pi⁻
    apion_reg.add("histTpcSignalData", "Specific energy loss (#pi^{-})", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}});
    apion_reg.add("histTofSignalData", "TOF signal (#pi^{-})", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    apion_reg.add("histDcaVsPtData", "dcaXY vs Pt (#pi^{-})", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    apion_reg.add("histDcaZVsPtData", "dcaZ vs Pt (#pi^{-})", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    apion_reg.add("histTOFm2", "TOF m^2 vs Pt (#pi^{-})", HistType::kTH2F, {pAxis, {400, 0.0, 10.0, "m^2"}});
    apion_reg.add("histTpcNsigmaData", "n-sigma TPC (#pi^{-})", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{#pi^{+}}"}});
    apion_reg.add("histTofNsigmaData", "n-sigma TOF (#pi^{-})", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{#pi^{+}}"}});
    apion_reg.add("histNClusterTPC", "Number of Clusters in TPC vs Pt (#pi^{-})", HistType::kTH2F, {ptAxis, {160, 0.0, 160.0, "nCluster"}});
    apion_reg.add("histNClusterITS", "Number of Clusters in ITS vs Pt (#pi^{-})", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    apion_reg.add("histNClusterITSib", "Number of Clusters in ib of ITS vs Pt (#pi^{-})", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    apion_reg.add("histChi2TPC", "chi^2 TPC vs Pt (#pi^{-})", HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    apion_reg.add("histChi2ITS", "chi^2 ITS vs Pt (#pi^{-})", HistType::kTH2F, {ptAxis, {500, 0.0, 50.0, "chi^2"}});
    apion_reg.add("histTpcNsigmaData_cent", "n-sigma TPC (#pi^{-}) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{#pi^{+}}"}, centralityAxis});
    apion_reg.add("histTofNsigmaData_cent", "n-sigma TOF (#pi^{-}) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{#pi^{+}}"}, centralityAxis});
    apion_reg.add("histTofm2_cent", "mass^2 TOF (#pi^{-}) centrality", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{#pi^{+}}"}, centralityAxis});
    apion_reg.add("histTpcNsigmaData_eta", "n-sigma TPC (#pi^{-}) vs eta", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{#pi^{+}}"}, etaAxis});
    apion_reg.add("histTofNsigmaData_eta", "n-sigma TOF (#pi^{-}) vs eta", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{#pi^{+}}"}, etaAxis});
    apion_reg.add("histTofm2_eta", "mass^2 TOF (#pi^{-}) vs eta", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{#pi^{+}}"}, etaAxis});

    // histograms for Proton
    proton_reg.add("histTpcSignalData", "Specific energy loss (p)", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}});
    proton_reg.add("histTofSignalData", "TOF signal (p)", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    proton_reg.add("histDcaVsPtData", "dcaXY vs Pt (p)", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    proton_reg.add("histDcaZVsPtData", "dcaZ vs Pt (p)", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    proton_reg.add("histTOFm2", "TOF m^2 vs Pt (p)", HistType::kTH2F, {pAxis, {400, 0.0, 10.0, "m^2"}});
    proton_reg.add("histTpcNsigmaData", "n-sigma TPC (p)", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{p}"}});
    proton_reg.add("histTofNsigmaData", "n-sigma TOF (p)", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{p}"}});
    proton_reg.add("histNClusterTPC", "Number of Clusters in TPC vs Pt (p)", HistType::kTH2F, {ptAxis, {160, 0.0, 160.0, "nCluster"}});
    proton_reg.add("histNClusterITS", "Number of Clusters in ITS vs Pt (p)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    proton_reg.add("histNClusterITSib", "Number of Clusters in ib of ITS vs Pt (p)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    proton_reg.add("histChi2TPC", "chi^2 TPC vs Pt (p)", HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    proton_reg.add("histChi2ITS", "chi^2 ITS vs Pt (p)", HistType::kTH2F, {ptAxis, {500, 0.0, 50.0, "chi^2"}});
    proton_reg.add("histTpcNsigmaData_cent", "n-sigma TPC (p) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{p}"}, centralityAxis});
    proton_reg.add("histTofNsigmaData_cent", "n-sigma TOF (p) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{p}"}, centralityAxis});
    proton_reg.add("histTofm2_cent", "mass^2 TOF (p) centrality", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{p}"}, centralityAxis});
    proton_reg.add("histTpcNsigmaData_eta", "n-sigma TPC (p) vs eta", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{p}"}, etaAxis});
    proton_reg.add("histTofNsigmaData_eta", "n-sigma TOF (p) vs eta", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{p}"}, etaAxis});
    proton_reg.add("histTofm2_eta", "mass^2 TOF (p) vs eta", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{p}"}, etaAxis});

    // histograms for antiProton
    aproton_reg.add("histTpcSignalData", "Specific energy loss (#bar{p})", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}});
    aproton_reg.add("histTofSignalData", "TOF signal (#bar{p})", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    aproton_reg.add("histDcaVsPtData", "dcaXY vs Pt (#bar{p})", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    aproton_reg.add("histDcaZVsPtData", "dcaZ vs Pt (#bar{p})", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    aproton_reg.add("histTOFm2", "TOF m^2 vs Pt (#bar{p})", HistType::kTH2F, {pAxis, {400, 0.0, 10.0, "m^2"}});
    aproton_reg.add("histTpcNsigmaData", "n-sigma TPC (#bar{p})", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{#bar{p}}"}});
    aproton_reg.add("histTofNsigmaData", "n-sigma TOF (#bar{p})", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{#bar{p}}"}});
    aproton_reg.add("histNClusterTPC", "Number of Clusters in TPC vs Pt (#bar{p})", HistType::kTH2F, {ptAxis, {160, 0.0, 160.0, "nCluster"}});
    aproton_reg.add("histNClusterITS", "Number of Clusters in ITS vs Pt (#bar{p})", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    aproton_reg.add("histNClusterITSib", "Number of Clusters in ib of ITS vs Pt (#bar{p})", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    aproton_reg.add("histChi2TPC", "chi^2 TPC vs Pt (#bar{p})", HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    aproton_reg.add("histChi2ITS", "chi^2 ITS vs Pt (#bar{p})", HistType::kTH2F, {ptAxis, {500, 0.0, 50.0, "chi^2"}});
    aproton_reg.add("histTpcNsigmaData_cent", "n-sigma TPC (#bar{p}) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{#bar{p}}"}, centralityAxis});
    aproton_reg.add("histTofNsigmaData_cent", "n-sigma TOF (#bar{p}) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{#bar{p}}"}, centralityAxis});
    aproton_reg.add("histTofm2_cent", "mass^2 TOF (#bar{p}) centrality", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{#bar{p}}"}, centralityAxis});
    aproton_reg.add("histTpcNsigmaData_eta", "n-sigma TPC (#bar{p}) vs eta", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{#bar{p}}"}, etaAxis});
    aproton_reg.add("histTofNsigmaData_eta", "n-sigma TOF (#bar{p}) vs eta", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{#bar{p}}"}, etaAxis});
    aproton_reg.add("histTofm2_eta", "mass^2 TOF (#bar{p}) vs eta", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{#bar{p}}"}, etaAxis});

    // histograms for Deuterons
    deuteron_reg.add("histTpcSignalData", "Specific energy loss (d)", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}});
    deuteron_reg.add("histTofSignalData", "TOF signal (d)", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    deuteron_reg.add("histDcaVsPtData", "dcaXY vs Pt (d)", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    deuteron_reg.add("histDcaZVsPtData", "dcaZ vs Pt (d)", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    deuteron_reg.add("histTOFm2", "TOF m^2 vs Pt (d)", HistType::kTH2F, {pAxis, {400, 0.0, 10.0, "m^2"}});
    deuteron_reg.add("histTpcNsigmaData", "n-sigma TPC (d)", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{d}"}});
    deuteron_reg.add("histTofNsigmaData", "n-sigma TOF (d)", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{d}"}});
    deuteron_reg.add("histNClusterTPC", "Number of Clusters in TPC vs Pt (d)", HistType::kTH2F, {ptAxis, {160, 0.0, 160.0, "nCluster"}});
    deuteron_reg.add("histNClusterITS", "Number of Clusters in ITS vs Pt (d)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    deuteron_reg.add("histNClusterITSib", "Number of Clusters in ib of ITS vs Pt (d)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    deuteron_reg.add("histChi2TPC", "chi^2 TPC vs Pt (d)", HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    deuteron_reg.add("histChi2ITS", "chi^2 ITS vs Pt (d)", HistType::kTH2F, {ptAxis, {500, 0.0, 50.0, "chi^2"}});
    deuteron_reg.add("histTpcNsigmaData_cent", "n-sigma TPC (d) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{d}"}, centralityAxis});
    deuteron_reg.add("histTofNsigmaData_cent", "n-sigma TOF (d) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{d}"}, centralityAxis});
    deuteron_reg.add("histTofm2_cent", "mass^2 TOF (d) centrality", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{d}"}, centralityAxis});
    deuteron_reg.add("histTpcNsigmaData_eta", "n-sigma TPC (d) vs eta", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{d}"}, etaAxis});
    deuteron_reg.add("histTofNsigmaData_eta", "n-sigma TOF (d) vs eta", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{d}"}, etaAxis});
    deuteron_reg.add("histTofm2_eta", "mass^2 TOF (d) vs eta", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{d}"}, etaAxis});

    // histograms for antiDeuterons
    adeuteron_reg.add("histTpcSignalData", "Specific energy loss (#bar{d})", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}});
    adeuteron_reg.add("histTofSignalData", "TOF signal (#bar{d})", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    adeuteron_reg.add("histDcaVsPtData", "dcaXY vs Pt (#bar{d})", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    adeuteron_reg.add("histDcaZVsPtData", "dcaZ vs Pt (#bar{d})", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    adeuteron_reg.add("histTOFm2", "TOF m^2 vs Pt (#bar{d})", HistType::kTH2F, {pAxis, {400, 0.0, 10.0, "m^2"}});
    adeuteron_reg.add("histTpcNsigmaData", "n-sigma TPC (#bar{d})", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{#bar{d}}"}});
    adeuteron_reg.add("histTofNsigmaData", "n-sigma TOF (#bar{d})", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{#bar{d}}"}});
    adeuteron_reg.add("histNClusterTPC", "Number of Clusters in TPC vs Pt (#bar{d})", HistType::kTH2F, {ptAxis, {160, 0.0, 160.0, "nCluster"}});
    adeuteron_reg.add("histNClusterITS", "Number of Clusters in ITS vs Pt (#bar{d})", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    adeuteron_reg.add("histNClusterITSib", "Number of Clusters in ib of ITS vs Pt (#bar{d})", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    adeuteron_reg.add("histChi2TPC", "chi^2 TPC vs Pt (#bar{d})", HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    adeuteron_reg.add("histChi2ITS", "chi^2 ITS vs Pt (#bar{d})", HistType::kTH2F, {ptAxis, {500, 0.0, 50.0, "chi^2"}});
    adeuteron_reg.add("histTpcNsigmaData_cent", "n-sigma TPC (#bar{d}) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{#bar{d}}"}, centralityAxis});
    adeuteron_reg.add("histTofNsigmaData_cent", "n-sigma TOF (#bar{d}) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{#bar{d}}"}, centralityAxis});
    adeuteron_reg.add("histTofm2_cent", "mass^2 TOF (#bar{d}) centrality", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{#bar{d}}"}, centralityAxis});
    adeuteron_reg.add("histTpcNsigmaData_eta", "n-sigma TPC (#bar{d}) vs eta", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{#bar{d}}"}, etaAxis});
    adeuteron_reg.add("histTofNsigmaData_eta", "n-sigma TOF (#bar{d}) vs eta", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{#bar{d}}"}, etaAxis});
    adeuteron_reg.add("histTofm2_eta", "mass^2 TOF (#bar{d}) vs eta", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{#bar{d}}"}, etaAxis});

    // histograms for Triton
    triton_reg.add("histTpcSignalData", "Specific energy loss (t)", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}});
    triton_reg.add("histTofSignalData", "TOF signal (t)", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    triton_reg.add("histDcaVsPtData", "dcaXY vs Pt (t)", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    triton_reg.add("histDcaZVsPtData", "dcaZ vs Pt (t)", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    triton_reg.add("histTOFm2", "TOF m^2 vs Pt (t)", HistType::kTH2F, {pAxis, {400, 0.0, 10.0, "m^2"}});
    triton_reg.add("histTpcNsigmaData", "n-sigma TPC (t)", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{t}"}});
    triton_reg.add("histTofNsigmaData", "n-sigma TOF (t)", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{t}"}});
    triton_reg.add("histNClusterTPC", "Number of Clusters in TPC vs Pt (t)", HistType::kTH2F, {ptAxis, {160, 0.0, 160.0, "nCluster"}});
    triton_reg.add("histNClusterITS", "Number of Clusters in ITS vs Pt (t)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    triton_reg.add("histNClusterITSib", "Number of Clusters in ib of ITS vs Pt (t)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    triton_reg.add("histChi2TPC", "chi^2 TPC vs Pt (t)", HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    triton_reg.add("histChi2ITS", "chi^2 ITS vs Pt (t)", HistType::kTH2F, {ptAxis, {500, 0.0, 50.0, "chi^2"}});
    triton_reg.add("histTpcNsigmaData_cent", "n-sigma TPC (t) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{t}"}, centralityAxis});
    triton_reg.add("histTofNsigmaData_cent", "n-sigma TOF (t) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{t}"}, centralityAxis});
    triton_reg.add("histTofm2_cent", "mass^2 TOF (t) centrality", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{t}"}, centralityAxis});
    triton_reg.add("histTpcNsigmaData_eta", "n-sigma TPC (t) vs eta", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{t}"}, etaAxis});
    triton_reg.add("histTofNsigmaData_eta", "n-sigma TOF (t) vs eta", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{t}"}, etaAxis});
    triton_reg.add("histTofm2_eta", "mass^2 TOF (t) vs eta", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{t}"}, etaAxis});

    // histograms for antiTriton
    atriton_reg.add("histTpcSignalData", "Specific energy loss (#bar{t})", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}});
    atriton_reg.add("histTofSignalData", "TOF signal (#bar{t})", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    atriton_reg.add("histDcaVsPtData", "dcaXY vs Pt (#bar{t})", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    atriton_reg.add("histDcaZVsPtData", "dcaZ vs Pt (#bar{t})", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    atriton_reg.add("histTOFm2", "TOF m^2 vs Pt (#bar{t})", HistType::kTH2F, {pAxis, {400, 0.0, 10.0, "m^2"}});
    atriton_reg.add("histTpcNsigmaData", "n-sigma TPC (#bar{t})", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{#bar{t}}"}});
    atriton_reg.add("histTofNsigmaData", "n-sigma TOF (#bar{t})", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{#bar{t}}"}});
    atriton_reg.add("histNClusterTPC", "Number of Clusters in TPC vs Pt (#bar{t})", HistType::kTH2F, {ptAxis, {160, 0.0, 160.0, "nCluster"}});
    atriton_reg.add("histNClusterITS", "Number of Clusters in ITS vs Pt (#bar{t})", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    atriton_reg.add("histNClusterITSib", "Number of Clusters in ib of ITS vs Pt (#bar{t})", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    atriton_reg.add("histChi2TPC", "chi^2 TPC vs Pt (#bar{t})", HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    atriton_reg.add("histChi2ITS", "chi^2 ITS vs Pt (#bar{t})", HistType::kTH2F, {ptAxis, {500, 0.0, 50.0, "chi^2"}});
    atriton_reg.add("histTpcNsigmaData_cent", "n-sigma TPC (#bar{t}) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{#bar{t}}"}, centralityAxis});
    atriton_reg.add("histTofNsigmaData_cent", "n-sigma TOF (#bar{t}) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{#bar{t}}"}, centralityAxis});
    atriton_reg.add("histTofm2_cent", "mass^2 TOF (#bar{t}) centrality", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{#bar{t}}"}, centralityAxis});
    atriton_reg.add("histTpcNsigmaData_eta", "n-sigma TPC (#bar{t}) vs eta", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{#bar{t}}"}, etaAxis});
    atriton_reg.add("histTofNsigmaData_eta", "n-sigma TOF (#bar{t}) vs eta", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{#bar{t}}"}, etaAxis});
    atriton_reg.add("histTofm2_eta", "mass^2 TOF (#bar{t}) vs eta", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{#bar{t}}"}, etaAxis});

    // histograms for Helium-3
    Helium3_reg.add("histTpcSignalData", "Specific energy loss (^{3}He)", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}});
    Helium3_reg.add("histTofSignalData", "TOF signal (^{3}He)", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    Helium3_reg.add("histDcaVsPtData", "dcaXY vs Pt (^{3}He)", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    Helium3_reg.add("histDcaZVsPtData", "dcaZ vs Pt (^{3}He)", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    Helium3_reg.add("histTOFm2", "TOF m^2 vs Pt (^{3}He)", HistType::kTH2F, {pAxis, {400, 0.0, 10.0, "m^2"}});
    Helium3_reg.add("histTpcNsigmaData", "n-sigma TPC (^{3}He)", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{^{3}He}"}});
    Helium3_reg.add("histTofNsigmaData", "n-sigma TOF (^{3}He)", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{^{3}He}"}});
    Helium3_reg.add("histNClusterTPC", "Number of Clusters in TPC vs Pt (^{3}He)", HistType::kTH2F, {ptAxis, {160, 0.0, 160.0, "nCluster"}});
    Helium3_reg.add("histNClusterITS", "Number of Clusters in ITS vs Pt (^{3}He)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    Helium3_reg.add("histNClusterITSib", "Number of Clusters in ib of ITS vs Pt (^{3}He)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    Helium3_reg.add("histChi2TPC", "chi^2 TPC vs Pt (^{3}He)", HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    Helium3_reg.add("histChi2ITS", "chi^2 ITS vs Pt (^{3}He)", HistType::kTH2F, {ptAxis, {500, 0.0, 50.0, "chi^2"}});
    Helium3_reg.add("histTpcNsigmaData_cent", "n-sigma TPC (^{3}He) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{^{3}He}"}, centralityAxis});
    Helium3_reg.add("histTofNsigmaData_cent", "n-sigma TOF (^{3}He) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{^{3}He}"}, centralityAxis});
    Helium3_reg.add("histTofm2_cent", "mass^2 TOF (^{3}He) centrality", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{^{3}He}"}, centralityAxis});
    Helium3_reg.add("histTpcNsigmaData_eta", "n-sigma TPC (^{3}He) vs eta", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{^{3}He}"}, etaAxis});
    Helium3_reg.add("histTofNsigmaData_eta", "n-sigma TOF (^{3}He) vs eta", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{^{3}He}"}, etaAxis});
    Helium3_reg.add("histTofm2_eta", "mass^2 TOF (^{3}He) vs eta", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{^{3}He}"}, etaAxis});

    // histograms for antiHelium-3
    aHelium3_reg.add("histTpcSignalData", "Specific energy loss (^{3}#bar{He})", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}});
    aHelium3_reg.add("histTofSignalData", "TOF signal (^{3}#bar{He})", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    aHelium3_reg.add("histDcaVsPtData", "dcaXY vs Pt (^{3}#bar{He})", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    aHelium3_reg.add("histDcaZVsPtData", "dcaZ vs Pt (^{3}#bar{He})", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    aHelium3_reg.add("histTOFm2", "TOF m^2 vs Pt (^{3}#bar{He})", HistType::kTH2F, {pAxis, {400, 0.0, 10.0, "m^2"}});
    aHelium3_reg.add("histTpcNsigmaData", "n-sigma TPC (^{3}#bar{He})", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{^{3}He}"}});
    aHelium3_reg.add("histTofNsigmaData", "n-sigma TOF (^{3}#bar{He})", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{^{3}He}"}});
    aHelium3_reg.add("histNClusterTPC", "Number of Clusters in TPC vs Pt (^{3}#bar{He})", HistType::kTH2F, {ptAxis, {160, 0.0, 160.0, "nCluster"}});
    aHelium3_reg.add("histNClusterITS", "Number of Clusters in ITS vs Pt (^{3}#bar{He})", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    aHelium3_reg.add("histNClusterITSib", "Number of Clusters in ib of ITS vs Pt (^{3}#bar{He})", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    aHelium3_reg.add("histChi2TPC", "chi^2 TPC vs Pt (^{3}#bar{He})", HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    aHelium3_reg.add("histChi2ITS", "chi^2 ITS vs Pt (^{3}#bar{He})", HistType::kTH2F, {ptAxis, {500, 0.0, 50.0, "chi^2"}});
    aHelium3_reg.add("histTpcNsigmaData_cent", "n-sigma TPC (^{3}#bar{He}) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{^{3}He}"}, centralityAxis});
    aHelium3_reg.add("histTofNsigmaData_cent", "n-sigma TOF (^{3}#bar{He}) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{^{3}He}"}, centralityAxis});
    aHelium3_reg.add("histTofm2_cent", "mass^2 TOF (^{3}#bar{He}) centrality", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{^{3}He}"}, centralityAxis});
    aHelium3_reg.add("histTpcNsigmaData_eta", "n-sigma TPC (^{3}#bar{He}) vs eta", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{^{3}He}"}, etaAxis});
    aHelium3_reg.add("histTofNsigmaData_eta", "n-sigma TOF (^{3}#bar{He}) vs eta", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{^{3}He}"}, etaAxis});
    aHelium3_reg.add("histTofm2_eta", "mass^2 TOF (^{3}#bar{He}) vs eta", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{^{3}He}"}, etaAxis});

    // histograms for Helium-4 (alpha)
    Helium4_reg.add("histTpcSignalData", "Specific energy loss (^{4}He)", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}});
    Helium4_reg.add("histTofSignalData", "TOF signal (^{4}He)", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    Helium4_reg.add("histDcaVsPtData", "dcaXY vs Pt (^{4}He)", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    Helium4_reg.add("histDcaZVsPtData", "dcaZ vs Pt (^{4}He)", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    Helium4_reg.add("histTOFm2", "TOF m^2 vs Pt (^{4}He)", HistType::kTH2F, {pAxis, {400, 0.0, 10.0, "m^2"}});
    Helium4_reg.add("histTpcNsigmaData", "n-sigma TPC (^{4}He)", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{^{4}He}"}});
    Helium4_reg.add("histTofNsigmaData", "n-sigma TOF (^{4}He)", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{^{4}He}"}});
    Helium4_reg.add("histNClusterTPC", "Number of Clusters in TPC vs Pt (^{4}He)", HistType::kTH2F, {ptAxis, {160, 0.0, 160.0, "nCluster"}});
    Helium4_reg.add("histNClusterITS", "Number of Clusters in ITS vs Pt (^{4}He)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    Helium4_reg.add("histNClusterITSib", "Number of Clusters in ib of ITS vs Pt (^{4}He)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    Helium4_reg.add("histChi2TPC", "chi^2 TPC vs Pt (^{4}He)", HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    Helium4_reg.add("histChi2ITS", "chi^2 ITS vs Pt (^{4}He)", HistType::kTH2F, {ptAxis, {500, 0.0, 50.0, "chi^2"}});
    Helium4_reg.add("histTpcNsigmaData_cent", "n-sigma TPC (^{4}He) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{^{4}He}"}, centralityAxis});
    Helium4_reg.add("histTofNsigmaData_cent", "n-sigma TOF (^{4}He) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{^{4}He}"}, centralityAxis});
    Helium4_reg.add("histTofm2_cent", "mass^2 TOF (^{4}He) centrality", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{^{4}He}"}, centralityAxis});
    Helium4_reg.add("histTpcNsigmaData_eta", "n-sigma TPC (^{4}He) vs eta", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{^{4}He}"}, etaAxis});
    Helium4_reg.add("histTofNsigmaData_eta", "n-sigma TOF (^{4}He) vs eta", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{^{4}He}"}, etaAxis});
    Helium4_reg.add("histTofm2_eta", "mass^2 TOF (^{4}He) vs eta", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{^{4}He}"}, etaAxis});

    // histograms for antiHelium-4 (alpha)
    aHelium4_reg.add("histTpcSignalData", "Specific energy loss (^{4}#bar{He})", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}});
    aHelium4_reg.add("histTofSignalData", "TOF signal (^{4}#bar{He})", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    aHelium4_reg.add("histDcaVsPtData", "dcaXY vs Pt (^{4}#bar{He})", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    aHelium4_reg.add("histDcaZVsPtData", "dcaZ vs Pt (^{4}#bar{He})", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    aHelium4_reg.add("histTOFm2", "TOF m^2 vs Pt (^{4}#bar{He})", HistType::kTH2F, {pAxis, {400, 0.0, 10.0, "m^2"}});
    aHelium4_reg.add("histTpcNsigmaData", "n-sigma TPC (^{4}#bar{He})", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{antiHe-4}"}});
    aHelium4_reg.add("histTofNsigmaData", "n-sigma TOF (^{4}#bar{He})", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{antiHe-4}"}});
    aHelium4_reg.add("histNClusterTPC", "Number of Clusters in TPC vs Pt (^{4}#bar{He})", HistType::kTH2F, {ptAxis, {160, 0.0, 160.0, "nCluster"}});
    aHelium4_reg.add("histNClusterITS", "Number of Clusters in ITS vs Pt (^{4}#bar{He})", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    aHelium4_reg.add("histNClusterITSib", "Number of Clusters in ib of ITS vs Pt (^{4}#bar{He})", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    aHelium4_reg.add("histChi2TPC", "chi^2 TPC vs Pt (^{4}#bar{He})", HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    aHelium4_reg.add("histChi2ITS", "chi^2 ITS vs Pt (^{4}#bar{He})", HistType::kTH2F, {ptAxis, {500, 0.0, 50.0, "chi^2"}});
    aHelium4_reg.add("histTpcNsigmaData_cent", "n-sigma TPC (^{4}#bar{He}) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{antiHe-4}"}, centralityAxis});
    aHelium4_reg.add("histTofNsigmaData_cent", "n-sigma TOF (^{4}#bar{He}) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{antiHe-4}"}, centralityAxis});
    aHelium4_reg.add("histTofm2_cent", "mass^2 TOF (^{4}#bar{He}) centrality", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{antiHe-4}"}, centralityAxis});
    aHelium4_reg.add("histTpcNsigmaData_eta", "n-sigma TPC (^{4}#bar{He}) vs eta", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{p}"}, etaAxis});
    aHelium4_reg.add("histTofNsigmaData_eta", "n-sigma TOF (^{4}#bar{He}) vs eta", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{p}"}, etaAxis});
    aHelium4_reg.add("histTofm2_eta", "mass^2 TOF (^{4}#bar{He}) vs eta", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{p}"}, etaAxis});

    // +++++++++++++++++++++ MC ++++++++++++++++++++++++++

    // MC reconstructed
    MC_recon_reg.add("histRecVtxMC", "MC reconstructed vertex z position", HistType::kTH1F, {{400, -40., +40., "z position (cm)"}});
    MC_recon_reg.add("histCentrality", "Centrality", HistType::kTH1F, {centralityAxis});
    MC_recon_reg.add("histEta", "#eta", HistType::kTH2F, {{102, -2.01, 2.01}, PDGBINNING});
    MC_recon_reg.add("histPt", "p_{t}", HistType::kTH2F, {ptAxis, PDGBINNING});
    MC_recon_reg.add("histDCA", "DCA xy", HistType::kTH3F, {ptAxis, {250, -0.5, 0.5, "dca"}, PDGBINNING});
    MC_recon_reg.add("histDCAz", "DCA z", HistType::kTH3F, {ptAxis, {1000, -2.0, 2.0, "dca"}, PDGBINNING});
    MC_recon_reg.add("histTpcSignalData", "Specific energy loss", HistType::kTH3F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}, PDGBINNING});
    MC_recon_reg.add("histTofSignalData", "TOF signal", HistType::kTH3F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}, PDGBINNING});
    MC_recon_reg.add("histTpcSignalData_all_species", "Specific energy loss", HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}});
    MC_recon_reg.add("histTofSignalData_all_species", "TOF signal", HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    MC_recon_reg.add("histTOFm2", "TOF m^2 vs Pt", HistType::kTH3F, {pAxis, {400, 0.0, 10.0, "m^2"}, PDGBINNING});
    MC_recon_reg.add("histNClusterTPC", "Number of Clusters in TPC vs Pt", HistType::kTH3F, {ptAxis, {80, 0.0, 160.0, "nCluster"}, PDGBINNING});
    MC_recon_reg.add("histNClusterITS", "Number of Clusters in ITS vs Pt", HistType::kTH3F, {ptAxis, {10, 0.0, 10.0, "ITS nCls"}, PDGBINNING});
    MC_recon_reg.add("histNClusterITSib", "Number of Clusters in ib of ITS vs Pt", HistType::kTH3F, {ptAxis, {10, 0.0, 10.0, "ITS ib nCls"}, PDGBINNING});
    MC_recon_reg.add("histChi2TPC", "chi^2 TPC vs Pt", HistType::kTH3F, {ptAxis, {100, 0.0, 5.0, "chi^2"}, PDGBINNING});
    MC_recon_reg.add("histChi2ITS", "chi^2 ITS vs Pt", HistType::kTH3F, {ptAxis, {500, 0.0, 50.0, "chi^2"}, PDGBINNING});
    MC_recon_reg.add("histChi2TOF", "chi^2 TOF vs Pt", HistType::kTH3F, {ptAxis, {75, 0.0, 15.0, "chi^2"}, PDGBINNING});
    MC_recon_reg.add("histTpcNsigmaDataPi", "TPC nSigma (#pi^{+})", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{#pi^{+}}"}});
    MC_recon_reg.add("histTpcNsigmaDataaPi", "TPC nSigma (#pi^{-})", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{#pi^{+}}"}});
    MC_recon_reg.add("histTpcNsigmaDataPr", "TPC nSigma (p)", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{p}"}});
    MC_recon_reg.add("histTpcNsigmaDataaPr", "TPC nSigma (#bar{p})", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{p}"}});
    MC_recon_reg.add("histTpcNsigmaDataDe", "TPC nSigma (d)", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{d}"}});
    MC_recon_reg.add("histTpcNsigmaDataaDe", "TPC nSigma (#bar{d})", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{d}"}});
    MC_recon_reg.add("histTpcNsigmaDataTr", "TPC nSigma (t)", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{t}"}});
    MC_recon_reg.add("histTpcNsigmaDataaTr", "TPC nSigma (#bar{t})", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{t}"}});
    MC_recon_reg.add("histTpcNsigmaDataHe", "TPC nSigma (^{3}He)", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{He-3}"}});
    MC_recon_reg.add("histTpcNsigmaDataaHe", "TPC nSigma (^{3}#bar{He})", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{He-3}"}});
    MC_recon_reg.add("histTpcNsigmaDataAl", "TPC nSigma (^{4}He)", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{He-4}"}});
    MC_recon_reg.add("histTpcNsigmaDataaAl", "TPC nSigma (^{4}#bar{He})", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{He-4}"}});
    MC_recon_reg.add("histTofNsigmaDataPi", "TOF nSigma (#pi^{+})", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{#pi^{+}}"}});
    MC_recon_reg.add("histTofNsigmaDataaPi", "TOF nSigma (#pi^{-})", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{#pi^{+}}"}});
    MC_recon_reg.add("histTofNsigmaDataPr", "TOF nSigma (p)", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{p}"}});
    MC_recon_reg.add("histTofNsigmaDataaPr", "TOF nSigma (#bar{p})", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{p}"}});
    MC_recon_reg.add("histTofNsigmaDataDe", "TOF nSigma (d)", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{d}"}});
    MC_recon_reg.add("histTofNsigmaDataaDe", "TOF nSigma (#bar{d})", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{d}"}});
    MC_recon_reg.add("histTofNsigmaDataTr", "TOF nSigma (t)", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{t}"}});
    MC_recon_reg.add("histTofNsigmaDataaTr", "TOF nSigma (#bar{t})", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{t}"}});
    MC_recon_reg.add("histTofNsigmaDataHe", "TOF nSigma (^{3}He)", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{^{3}He}"}});
    MC_recon_reg.add("histTofNsigmaDataaHe", "TOF nSigma (^{3}#bar{He})", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{^{3}He}"}});
    MC_recon_reg.add("histTofNsigmaDataAl", "TOF nSigma (^{4}He)", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{^{4}He}"}});
    MC_recon_reg.add("histTofNsigmaDataaAl", "TOF nSigma (^{4}#bar{He})", HistType::kTH2F, {pAxis, {160, -20., +20., "n#sigma_{^{4}He}"}});
  }

  // Configurables
  Configurable<int> use_momentum_getter{"use_momentum_getter", 0, "0: track.p(), 1: track.pt(), 2: track.tpcInnerParam()"};
  Configurable<int> momentum_He3{"momentum_He3", 0, "0: momentum * 1.0, 1: momentum * 2.0"};
  Configurable<float> yMin{"yMin", -0.5, "Maximum rapidity"};
  Configurable<float> yMax{"yMax", 0.5, "Minimum rapidity"};
  Configurable<float> p_min{"p_min", 0.1f, "min momentum"};
  Configurable<float> p_max{"p_max", 1e+10f, "max momentum"};
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.9f, "Eta range for tracks"};
  Configurable<float> nsigmacutLow{"nsigmacutLow", -3.0, "Value of the Nsigma cut"};
  Configurable<float> nsigmacutHigh{"nsigmacutHigh", +3.0, "Value of the Nsigma cut"};
  Configurable<float> minCentrality{"minCentrality", 0.0, "min Centrality used"};
  Configurable<float> maxCentrality{"maxCentrality", 90.0, "max Centrality used"};
  Configurable<bool> enable_Centrality_cut_global{"enable_Centrality_cut_global", true, "use Centrality cut"};

  // Replacement for globalTrack filter
  Configurable<bool> requireITS{"requireITS", true, "Additional cut on the ITS requirement"};
  Configurable<bool> requireTPC{"requireTPC", true, "Additional cut on the TPC requirement"};
  Configurable<bool> passedITSRefit{"passedITSRefit", true, "Additional cut on the ITS refit requirement"};
  Configurable<bool> passedTPCRefit{"passedTPCRefit", true, "Additional cut on the TPC refit requirement"};
  Configurable<float> minReqClusterITS{"minReqClusterITS", 1.0, "min number of clusters required in ITS"};
  Configurable<float> minReqClusterITSib{"minReqClusterITSib", 1.0, "min number of clusters required in ITS inner barrel"};
  Configurable<float> minTPCnClsFound{"minTPCnClsFound", 0.0f, "min number of crossed rows TPC"};
  Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 70.0f, "min number of crossed rows TPC"};
  Configurable<float> minRatioCrossedRowsTPC{"minRatioCrossedRowsTPC", 0.8f, "min ratio of crossed rows over findable clusters TPC"};
  Configurable<float> maxRatioCrossedRowsTPC{"maxRatioCrossedRowsTPC", 2.0f, "max ratio of crossed rows over findable clusters TPC"};
  Configurable<float> maxChi2PerClusterTPC{"maxChi2PerClusterTPC", 4.f, "Cut on the maximum value of the chi2 per cluster in the TPC"};
  Configurable<float> minChi2PerClusterTPC{"minChi2PerClusterTPC", 0.5f, "Cut on the minimum value of the chi2 per cluster in the TPC"};
  Configurable<float> maxChi2PerClusterITS{"maxChi2PerClusterITS", 36.f, "Cut on the maximum value of the chi2 per cluster in the ITS"};
  Configurable<float> maxDCA_XY{"maxDCA_XY", 0.5f, "max DCA to vertex xy"};
  Configurable<float> maxDCA_Z{"maxDCA_Z", 2.0f, "max DCA to vertex z"};
  Configurable<int> lastRequiredTrdCluster{"lastRequiredTrdCluster", -1, "Last cluster to required in TRD for track selection. -1 does not require any TRD cluster"};
  Configurable<bool> requireGoldenChi2{"requireGoldenChi2", false, "Enable the requirement of GoldenChi2"};
  Configurable<bool> event_selection_sel8{"event_selection_sel8", true, "Enable sel8 event selection"};
  Configurable<bool> event_selection_MC_sel8{"event_selection_MC_sel8", true, "Enable sel8 event selection in MC processing"};

  template <typename CollisionType, typename TracksType>
  void fillHistograms(const CollisionType& event, const TracksType& tracks)
  {

    if (event_selection_sel8 && event.sel8()) {
      spectra_reg.fill(HIST("histRecVtxZData"), event.posZ());
    }
    if (!event_selection_sel8) {
      spectra_reg.fill(HIST("histRecVtxZData"), event.posZ());
    }

    for (auto track : tracks) {

      histTrackcuts_data->AddBinContent(1);
      if (event_selection_sel8 && !event.sel8())
        continue;

      histTrackcuts_data->AddBinContent(2);
      double momentum = 0.0;
      double momentum_Z_equal_2 = 0.0;
      int momentum_getter = 0;
      int momentum_getter_Z_equal_2 = 0;
      switch (use_momentum_getter) {
        case 0:
          momentum = track.p();
          switch (momentum_He3) {
            case 0:
              momentum_Z_equal_2 = track.p();
              break;
            case 1:
              momentum_Z_equal_2 = track.p() * 2.0;
              break;
            default:
              momentum_getter_Z_equal_2 = -1;
              break;
          }
          break;
        case 1:
          momentum = track.pt();
          switch (momentum_He3) {
            case 0:
              momentum_Z_equal_2 = track.pt();
              break;
            case 1:
              momentum_Z_equal_2 = track.pt() * 2.0;
              break;
            default:
              momentum_getter_Z_equal_2 = -1;
              break;
          }
          break;
        case 2:
          momentum = track.tpcInnerParam();
          switch (momentum_He3) {
            case 0:
              momentum_Z_equal_2 = track.tpcInnerParam();
              break;
            case 1:
              momentum_Z_equal_2 = track.tpcInnerParam() * 2.0;
              break;
            default:
              momentum_getter_Z_equal_2 = -1;
              break;
          }
          break;
        default:
          momentum_getter = -1;
          break;
      }
      if (momentum_getter == -1 || momentum_getter_Z_equal_2 == -1)
        break;

      float TPCnumberClsFound = track.tpcNClsFound();
      float TPC_nCls_Crossed_Rows = track.tpcNClsCrossedRows();
      float RatioCrossedRowsOverFindableTPC = track.tpcCrossedRowsOverFindableCls();
      float Chi2perClusterTPC = track.tpcChi2NCl();
      float Chi2perClusterITS = track.itsChi2NCl();

      if (track.sign() > 0) {
        spectra_reg.fill(HIST("histDcaVsPtData_particle"), track.pt(), track.dcaXY());
        spectra_reg.fill(HIST("histDcaZVsPtData_particle"), track.pt(), track.dcaZ());
      }
      if (track.sign() < 0) {
        spectra_reg.fill(HIST("histDcaVsPtData_antiparticle"), track.pt(), track.dcaXY());
        spectra_reg.fill(HIST("histDcaZVsPtData_antiparticle"), track.pt(), track.dcaZ());
      }

      TLorentzVector lorentzVector_pion{};
      TLorentzVector lorentzVector_proton{};
      TLorentzVector lorentzVector_deuteron{};
      TLorentzVector lorentzVector_triton{};
      TLorentzVector lorentzVector_He3{};
      TLorentzVector lorentzVector_He4{};

      lorentzVector_pion.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassPiPlus);
      lorentzVector_proton.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassProton);
      lorentzVector_deuteron.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassDeuteron);
      lorentzVector_triton.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassTriton);
      lorentzVector_He3.SetPtEtaPhiM(track.pt() * 2.0, track.eta(), track.phi(), constants::physics::MassHelium3);
      lorentzVector_He4.SetPtEtaPhiM(track.pt() * 2.0, track.eta(), track.phi(), constants::physics::MassAlpha);

      if (lorentzVector_pion.Rapidity() < yMin || lorentzVector_pion.Rapidity() > yMax ||
          lorentzVector_proton.Rapidity() < yMin || lorentzVector_proton.Rapidity() > yMax ||
          lorentzVector_deuteron.Rapidity() < yMin || lorentzVector_deuteron.Rapidity() > yMax ||
          lorentzVector_triton.Rapidity() < yMin || lorentzVector_triton.Rapidity() > yMax ||
          lorentzVector_He3.Rapidity() < yMin || lorentzVector_He3.Rapidity() > yMax ||
          lorentzVector_He4.Rapidity() < yMin || lorentzVector_He4.Rapidity() > yMax)
        continue;
      histTrackcuts_data->AddBinContent(3);

      float nSigmaPion = track.tpcNSigmaPi();
      float nSigmaProton = track.tpcNSigmaPr();
      float nSigmaDeut = track.tpcNSigmaDe();
      float nSigmaTriton = track.tpcNSigmaTr();
      float nSigmaHe3 = track.tpcNSigmaHe();
      float nSigmaHe4 = track.tpcNSigmaAl();

      spectra_reg.fill(HIST("histNClusterTPC"), track.pt(), track.tpcNClsCrossedRows());
      spectra_reg.fill(HIST("histNClusterITS"), track.pt(), track.itsNCls());
      spectra_reg.fill(HIST("histNClusterITSib"), track.pt(), track.itsNClsInnerBarrel());
      spectra_reg.fill(HIST("histChi2TPC"), track.pt(), track.tpcChi2NCl());
      spectra_reg.fill(HIST("histChi2ITS"), track.pt(), track.itsChi2NCl());

      if (TMath::Abs(track.dcaXY()) > maxDCA_XY || TMath::Abs(track.dcaZ()) > maxDCA_Z)
        continue;
      histTrackcuts_data->AddBinContent(4);
      if (TPCnumberClsFound < minTPCnClsFound || TPC_nCls_Crossed_Rows < minNCrossedRowsTPC)
        continue;
      histTrackcuts_data->AddBinContent(5);
      if (RatioCrossedRowsOverFindableTPC < minRatioCrossedRowsTPC || RatioCrossedRowsOverFindableTPC > maxRatioCrossedRowsTPC)
        continue;
      histTrackcuts_data->AddBinContent(6);
      if (Chi2perClusterTPC > maxChi2PerClusterTPC || Chi2perClusterTPC < minChi2PerClusterTPC || Chi2perClusterITS > maxChi2PerClusterITS)
        continue;
      histTrackcuts_data->AddBinContent(7);
      if (!(track.passedTPCRefit()))
        continue;
      histTrackcuts_data->AddBinContent(8);
      if (!(track.passedITSRefit()))
        continue;
      histTrackcuts_data->AddBinContent(9);
      if ((track.itsNClsInnerBarrel()) < minReqClusterITSib || (track.itsNCls()) < minReqClusterITS)
        continue;
      histTrackcuts_data->AddBinContent(10);
      if (track.pt() < p_min || track.pt() > p_max)
        continue;
      histTrackcuts_data->AddBinContent(11);
      if ((requireITS && !(track.hasITS())) || (requireTPC && !(track.hasTPC())))
        continue;
      histTrackcuts_data->AddBinContent(12);
      if (requireGoldenChi2 && !(track.passedGoldenChi2()))
        continue;
      histTrackcuts_data->AddBinContent(13);

      spectra_reg.fill(HIST("histTpcSignalData"), momentum * track.sign(), track.tpcSignal());

      if (track.sign() > 0) {
        pion_reg.fill(HIST("histTpcNsigmaData"), momentum, nSigmaPion);
        proton_reg.fill(HIST("histTpcNsigmaData"), momentum, nSigmaProton);
        deuteron_reg.fill(HIST("histTpcNsigmaData"), momentum, nSigmaDeut);
        triton_reg.fill(HIST("histTpcNsigmaData"), momentum, nSigmaTriton);
        Helium3_reg.fill(HIST("histTpcNsigmaData"), momentum_Z_equal_2, nSigmaHe3);
        Helium4_reg.fill(HIST("histTpcNsigmaData"), momentum_Z_equal_2, nSigmaHe4);

        if (track.hasTOF()) {
          if (track.hasTRD() && (lastRequiredTrdCluster > 0)) {
            int lastLayer = 0;
            for (int l = 7; l >= 0; l--) {
              if (track.trdPattern() & (1 << l)) {
                lastLayer = l;
                break;
              }
            }
            if (lastLayer < lastRequiredTrdCluster)
              continue;
          }
          Float_t TOFmass2 = ((track.mass()) * (track.mass()));
          spectra_reg.fill(HIST("histTOFm2"), momentum, TOFmass2);
        }
      }

      if (track.sign() < 0) {
        apion_reg.fill(HIST("histTpcNsigmaData"), momentum, nSigmaPion);
        aproton_reg.fill(HIST("histTpcNsigmaData"), momentum, nSigmaProton);
        adeuteron_reg.fill(HIST("histTpcNsigmaData"), momentum, nSigmaDeut);
        atriton_reg.fill(HIST("histTpcNsigmaData"), momentum, nSigmaTriton);
        aHelium3_reg.fill(HIST("histTpcNsigmaData"), momentum_Z_equal_2, nSigmaHe3);
        aHelium4_reg.fill(HIST("histTpcNsigmaData"), momentum_Z_equal_2, nSigmaHe4);

        if (track.hasTOF()) {
          if (track.hasTRD() && (lastRequiredTrdCluster > 0)) {
            int lastLayer = 0;
            for (int l = 7; l >= 0; l--) {
              if (track.trdPattern() & (1 << l)) {
                lastLayer = l;
                break;
              }
            }
            if (lastLayer < lastRequiredTrdCluster)
              continue;
          }
          Float_t TOFmass2 = ((track.mass()) * (track.mass()));
          spectra_reg.fill(HIST("histTOFm2"), momentum, TOFmass2);
        }
      }

      //**************   Pion   *******************

      if (nSigmaPion > nsigmacutLow && nSigmaPion < nsigmacutHigh) {

        if (track.sign() > 0) {
          pion_reg.fill(HIST("histDcaVsPtData"), track.pt(), track.dcaXY());
          pion_reg.fill(HIST("histDcaZVsPtData"), track.pt(), track.dcaZ());
          pion_reg.fill(HIST("histTpcSignalData"), momentum, track.tpcSignal());
          pion_reg.fill(HIST("histNClusterTPC"), track.pt(), track.tpcNClsFound());
          pion_reg.fill(HIST("histNClusterITS"), track.pt(), track.itsNCls());
          pion_reg.fill(HIST("histNClusterITSib"), track.pt(), track.itsNClsInnerBarrel());
          pion_reg.fill(HIST("histChi2TPC"), track.pt(), track.tpcChi2NCl());
          pion_reg.fill(HIST("histChi2ITS"), track.pt(), track.itsChi2NCl());

          if (track.hasTOF()) {
            if (track.hasTRD() && (lastRequiredTrdCluster > 0)) {
              int lastLayer = 0;
              for (int l = 7; l >= 0; l--) {
                if (track.trdPattern() & (1 << l)) {
                  lastLayer = l;
                  break;
                }
              }
              if (lastLayer < lastRequiredTrdCluster)
                continue;
            }
            Float_t TOFmass2 = ((track.mass()) * (track.mass()));
            Float_t beta = track.beta();
            pion_reg.fill(HIST("histTOFm2"), momentum, TOFmass2);
            pion_reg.fill(HIST("histTofSignalData"), momentum, beta);
            pion_reg.fill(HIST("histTofNsigmaData"), momentum, track.tofNSigmaPi());
          }
        }

        if (track.sign() < 0) {
          apion_reg.fill(HIST("histDcaVsPtData"), track.pt(), track.dcaXY());
          apion_reg.fill(HIST("histDcaZVsPtData"), track.pt(), track.dcaZ());
          apion_reg.fill(HIST("histTpcSignalData"), momentum, track.tpcSignal());
          apion_reg.fill(HIST("histNClusterTPC"), track.pt(), track.tpcNClsFound());
          apion_reg.fill(HIST("histNClusterITS"), track.pt(), track.itsNCls());
          apion_reg.fill(HIST("histNClusterITSib"), track.pt(), track.itsNClsInnerBarrel());
          apion_reg.fill(HIST("histChi2TPC"), track.pt(), track.tpcChi2NCl());
          apion_reg.fill(HIST("histChi2ITS"), track.pt(), track.itsChi2NCl());

          if (track.hasTOF()) {
            if (track.hasTRD() && (lastRequiredTrdCluster > 0)) {
              int lastLayer = 0;
              for (int l = 7; l >= 0; l--) {
                if (track.trdPattern() & (1 << l)) {
                  lastLayer = l;
                  break;
                }
              }
              if (lastLayer < lastRequiredTrdCluster)
                continue;
            }
            Float_t TOFmass2 = ((track.mass()) * (track.mass()));
            Float_t beta = track.beta();
            apion_reg.fill(HIST("histTOFm2"), momentum, TOFmass2);
            apion_reg.fill(HIST("histTofSignalData"), momentum, beta);
            apion_reg.fill(HIST("histTofNsigmaData"), momentum, track.tofNSigmaPi());
          }
        }

        if (track.hasTOF()) {
          if (track.hasTRD() && (lastRequiredTrdCluster > 0)) {
            int lastLayer = 0;
            for (int l = 7; l >= 0; l--) {
              if (track.trdPattern() & (1 << l)) {
                lastLayer = l;
                break;
              }
            }
            if (lastLayer < lastRequiredTrdCluster)
              continue;
          }
          spectra_reg.fill(HIST("histTofSignalData"), momentum * track.sign(), track.beta());
        }
      }

      //**************   Proton   *******************

      if (nSigmaProton > nsigmacutLow && nSigmaProton < nsigmacutHigh) {

        if (track.sign() > 0) {
          proton_reg.fill(HIST("histDcaVsPtData"), track.pt(), track.dcaXY());
          proton_reg.fill(HIST("histDcaZVsPtData"), track.pt(), track.dcaZ());
          proton_reg.fill(HIST("histTpcSignalData"), momentum, track.tpcSignal());
          proton_reg.fill(HIST("histNClusterTPC"), track.pt(), track.tpcNClsFound());
          proton_reg.fill(HIST("histNClusterITS"), track.pt(), track.itsNCls());
          proton_reg.fill(HIST("histNClusterITSib"), track.pt(), track.itsNClsInnerBarrel());
          proton_reg.fill(HIST("histChi2TPC"), track.pt(), track.tpcChi2NCl());
          proton_reg.fill(HIST("histChi2ITS"), track.pt(), track.itsChi2NCl());

          if (track.hasTOF()) {
            if (track.hasTRD() && (lastRequiredTrdCluster > 0)) {
              int lastLayer = 0;
              for (int l = 7; l >= 0; l--) {
                if (track.trdPattern() & (1 << l)) {
                  lastLayer = l;
                  break;
                }
              }
              if (lastLayer < lastRequiredTrdCluster)
                continue;
            }
            Float_t TOFmass2 = ((track.mass()) * (track.mass()));
            Float_t beta = track.beta();
            proton_reg.fill(HIST("histTOFm2"), momentum, TOFmass2);
            proton_reg.fill(HIST("histTofSignalData"), momentum, beta);
            proton_reg.fill(HIST("histTofNsigmaData"), momentum, track.tofNSigmaPr());
          }
        }

        if (track.sign() < 0) {
          aproton_reg.fill(HIST("histDcaVsPtData"), track.pt(), track.dcaXY());
          aproton_reg.fill(HIST("histDcaZVsPtData"), track.pt(), track.dcaZ());
          aproton_reg.fill(HIST("histTpcSignalData"), momentum, track.tpcSignal());
          aproton_reg.fill(HIST("histNClusterTPC"), track.pt(), track.tpcNClsFound());
          aproton_reg.fill(HIST("histNClusterITS"), track.pt(), track.itsNCls());
          aproton_reg.fill(HIST("histNClusterITSib"), track.pt(), track.itsNClsInnerBarrel());
          aproton_reg.fill(HIST("histChi2TPC"), track.pt(), track.tpcChi2NCl());
          aproton_reg.fill(HIST("histChi2ITS"), track.pt(), track.itsChi2NCl());

          if (track.hasTOF()) {
            if (track.hasTRD() && (lastRequiredTrdCluster > 0)) {
              int lastLayer = 0;
              for (int l = 7; l >= 0; l--) {
                if (track.trdPattern() & (1 << l)) {
                  lastLayer = l;
                  break;
                }
              }
              if (lastLayer < lastRequiredTrdCluster)
                continue;
            }
            Float_t TOFmass2 = ((track.mass()) * (track.mass()));
            Float_t beta = track.beta();
            aproton_reg.fill(HIST("histTOFm2"), momentum, TOFmass2);
            aproton_reg.fill(HIST("histTofSignalData"), momentum, beta);
            aproton_reg.fill(HIST("histTofNsigmaData"), momentum, track.tofNSigmaPr());
          }
        }

        if (track.hasTOF()) {
          if (track.hasTRD() && (lastRequiredTrdCluster > 0)) {
            int lastLayer = 0;
            for (int l = 7; l >= 0; l--) {
              if (track.trdPattern() & (1 << l)) {
                lastLayer = l;
                break;
              }
            }
            if (lastLayer < lastRequiredTrdCluster)
              continue;
          }
          spectra_reg.fill(HIST("histTofSignalData"), momentum * track.sign(), track.beta());
        }
      }

      //**************   Deuteron   *******************

      if (nSigmaDeut > nsigmacutLow && nSigmaDeut < nsigmacutHigh) {

        if (track.sign() > 0) {
          deuteron_reg.fill(HIST("histDcaVsPtData"), track.pt(), track.dcaXY());
          deuteron_reg.fill(HIST("histDcaZVsPtData"), track.pt(), track.dcaZ());
          deuteron_reg.fill(HIST("histTpcSignalData"), momentum, track.tpcSignal());
          deuteron_reg.fill(HIST("histNClusterTPC"), track.pt(), track.tpcNClsFound());
          deuteron_reg.fill(HIST("histNClusterITS"), track.pt(), track.itsNCls());
          deuteron_reg.fill(HIST("histNClusterITSib"), track.pt(), track.itsNClsInnerBarrel());
          deuteron_reg.fill(HIST("histChi2TPC"), track.pt(), track.tpcChi2NCl());
          deuteron_reg.fill(HIST("histChi2ITS"), track.pt(), track.itsChi2NCl());

          if (track.hasTOF()) {
            if (track.hasTRD() && (lastRequiredTrdCluster > 0)) {
              int lastLayer = 0;
              for (int l = 7; l >= 0; l--) {
                if (track.trdPattern() & (1 << l)) {
                  lastLayer = l;
                  break;
                }
              }
              if (lastLayer < lastRequiredTrdCluster)
                continue;
            }
            Float_t TOFmass2 = ((track.mass()) * (track.mass()));
            Float_t beta = track.beta();
            deuteron_reg.fill(HIST("histTOFm2"), momentum, TOFmass2);
            deuteron_reg.fill(HIST("histTofSignalData"), momentum, beta);
            deuteron_reg.fill(HIST("histTofNsigmaData"), momentum, track.tofNSigmaDe());
          }
        }

        if (track.sign() < 0) {
          adeuteron_reg.fill(HIST("histDcaVsPtData"), track.pt(), track.dcaXY());
          adeuteron_reg.fill(HIST("histDcaZVsPtData"), track.pt(), track.dcaZ());
          adeuteron_reg.fill(HIST("histTpcSignalData"), momentum, track.tpcSignal());
          adeuteron_reg.fill(HIST("histNClusterTPC"), track.pt(), track.tpcNClsFound());
          adeuteron_reg.fill(HIST("histNClusterITS"), track.pt(), track.itsNCls());
          adeuteron_reg.fill(HIST("histNClusterITSib"), track.pt(), track.itsNClsInnerBarrel());
          adeuteron_reg.fill(HIST("histChi2TPC"), track.pt(), track.tpcChi2NCl());
          adeuteron_reg.fill(HIST("histChi2ITS"), track.pt(), track.itsChi2NCl());

          if (track.hasTOF()) {
            if (track.hasTRD() && (lastRequiredTrdCluster > 0)) {
              int lastLayer = 0;
              for (int l = 7; l >= 0; l--) {
                if (track.trdPattern() & (1 << l)) {
                  lastLayer = l;
                  break;
                }
              }
              if (lastLayer < lastRequiredTrdCluster)
                continue;
            }
            Float_t TOFmass2 = ((track.mass()) * (track.mass()));
            Float_t beta = track.beta();
            adeuteron_reg.fill(HIST("histTOFm2"), momentum, TOFmass2);
            adeuteron_reg.fill(HIST("histTofSignalData"), momentum, beta);
            adeuteron_reg.fill(HIST("histTofNsigmaData"), momentum, track.tofNSigmaDe());
          }
        }

        if (track.hasTOF()) {
          if (track.hasTRD() && (lastRequiredTrdCluster > 0)) {
            int lastLayer = 0;
            for (int l = 7; l >= 0; l--) {
              if (track.trdPattern() & (1 << l)) {
                lastLayer = l;
                break;
              }
            }
            if (lastLayer < lastRequiredTrdCluster)
              continue;
          }
          spectra_reg.fill(HIST("histTofSignalData"), momentum * track.sign(), track.beta());
        }
      }

      //**************   Triton   *******************

      if (nSigmaTriton > nsigmacutLow && nSigmaTriton < nsigmacutHigh) {

        if (track.sign() > 0) {
          triton_reg.fill(HIST("histDcaVsPtData"), track.pt(), track.dcaXY());
          triton_reg.fill(HIST("histDcaZVsPtData"), track.pt(), track.dcaZ());
          triton_reg.fill(HIST("histTpcSignalData"), momentum, track.tpcSignal());
          triton_reg.fill(HIST("histNClusterTPC"), track.pt(), track.tpcNClsFound());
          triton_reg.fill(HIST("histNClusterITS"), track.pt(), track.itsNCls());
          triton_reg.fill(HIST("histNClusterITSib"), track.pt(), track.itsNClsInnerBarrel());
          triton_reg.fill(HIST("histChi2TPC"), track.pt(), track.tpcChi2NCl());
          triton_reg.fill(HIST("histChi2ITS"), track.pt(), track.itsChi2NCl());

          if (track.hasTOF()) {
            if (track.hasTRD() && (lastRequiredTrdCluster > 0)) {
              int lastLayer = 0;
              for (int l = 7; l >= 0; l--) {
                if (track.trdPattern() & (1 << l)) {
                  lastLayer = l;
                  break;
                }
              }
              if (lastLayer < lastRequiredTrdCluster)
                continue;
            }
            Float_t TOFmass2 = ((track.mass()) * (track.mass()));
            Float_t beta = track.beta();
            triton_reg.fill(HIST("histTOFm2"), momentum, TOFmass2);
            triton_reg.fill(HIST("histTofSignalData"), momentum, beta);
            triton_reg.fill(HIST("histTofNsigmaData"), momentum, track.tofNSigmaTr());
          }
        }

        if (track.sign() < 0) {
          atriton_reg.fill(HIST("histDcaVsPtData"), track.pt(), track.dcaXY());
          atriton_reg.fill(HIST("histDcaZVsPtData"), track.pt(), track.dcaZ());
          atriton_reg.fill(HIST("histTpcSignalData"), momentum, track.tpcSignal());
          atriton_reg.fill(HIST("histNClusterTPC"), track.pt(), track.tpcNClsFound());
          atriton_reg.fill(HIST("histNClusterITS"), track.pt(), track.itsNCls());
          atriton_reg.fill(HIST("histNClusterITSib"), track.pt(), track.itsNClsInnerBarrel());
          atriton_reg.fill(HIST("histChi2TPC"), track.pt(), track.tpcChi2NCl());
          atriton_reg.fill(HIST("histChi2ITS"), track.pt(), track.itsChi2NCl());

          if (track.hasTOF()) {
            if (track.hasTRD() && (lastRequiredTrdCluster > 0)) {
              int lastLayer = 0;
              for (int l = 7; l >= 0; l--) {
                if (track.trdPattern() & (1 << l)) {
                  lastLayer = l;
                  break;
                }
              }
              if (lastLayer < lastRequiredTrdCluster)
                continue;
            }
            Float_t TOFmass2 = ((track.mass()) * (track.mass()));
            Float_t beta = track.beta();
            atriton_reg.fill(HIST("histTOFm2"), momentum, TOFmass2);
            atriton_reg.fill(HIST("histTofSignalData"), momentum, beta);
            atriton_reg.fill(HIST("histTofNsigmaData"), momentum, track.tofNSigmaTr());
          }
        }

        if (track.hasTOF()) {
          if (track.hasTRD() && (lastRequiredTrdCluster > 0)) {
            int lastLayer = 0;
            for (int l = 7; l >= 0; l--) {
              if (track.trdPattern() & (1 << l)) {
                lastLayer = l;
                break;
              }
            }
            if (lastLayer < lastRequiredTrdCluster)
              continue;
          }
          spectra_reg.fill(HIST("histTofSignalData"), momentum * track.sign(), track.beta());
        }
      }

      //**************   Helium-3   *******************

      if (nSigmaHe3 > nsigmacutLow && nSigmaHe3 < nsigmacutHigh) {

        if (track.sign() > 0) {
          Helium3_reg.fill(HIST("histDcaVsPtData"), track.pt() * 2.0, track.dcaXY());
          Helium3_reg.fill(HIST("histDcaZVsPtData"), track.pt() * 2.0, track.dcaZ());
          Helium3_reg.fill(HIST("histTpcSignalData"), momentum_Z_equal_2, track.tpcSignal());
          Helium3_reg.fill(HIST("histNClusterTPC"), track.pt() * 2.0, track.tpcNClsFound());
          Helium3_reg.fill(HIST("histNClusterITS"), track.pt() * 2.0, track.itsNCls());
          Helium3_reg.fill(HIST("histNClusterITSib"), track.pt() * 2.0, track.itsNClsInnerBarrel());
          Helium3_reg.fill(HIST("histChi2TPC"), track.pt() * 2.0, track.tpcChi2NCl());
          Helium3_reg.fill(HIST("histChi2ITS"), track.pt() * 2.0, track.itsChi2NCl());

          if (track.hasTOF()) {
            if (track.hasTRD() && (lastRequiredTrdCluster > 0)) {
              int lastLayer = 0;
              for (int l = 7; l >= 0; l--) {
                if (track.trdPattern() & (1 << l)) {
                  lastLayer = l;
                  break;
                }
              }
              if (lastLayer < lastRequiredTrdCluster)
                continue;
            }
            Float_t TOFmass2 = ((track.mass()) * (track.mass()));
            Float_t beta = track.beta();
            Helium3_reg.fill(HIST("histTOFm2"), momentum_Z_equal_2, TOFmass2);
            Helium3_reg.fill(HIST("histTofSignalData"), momentum_Z_equal_2, beta);
            Helium3_reg.fill(HIST("histTofNsigmaData"), momentum_Z_equal_2, track.tofNSigmaHe());
          }
        }

        if (track.sign() < 0) {
          aHelium3_reg.fill(HIST("histDcaVsPtData"), track.pt() * 2.0, track.dcaXY());
          aHelium3_reg.fill(HIST("histDcaZVsPtData"), track.pt() * 2.0, track.dcaZ());
          aHelium3_reg.fill(HIST("histTpcSignalData"), momentum_Z_equal_2, track.tpcSignal());
          aHelium3_reg.fill(HIST("histNClusterTPC"), track.pt() * 2.0, track.tpcNClsFound());
          aHelium3_reg.fill(HIST("histNClusterITS"), track.pt() * 2.0, track.itsNCls());
          aHelium3_reg.fill(HIST("histNClusterITSib"), track.pt() * 2.0, track.itsNClsInnerBarrel());
          aHelium3_reg.fill(HIST("histChi2TPC"), track.pt() * 2.0, track.tpcChi2NCl());
          aHelium3_reg.fill(HIST("histChi2ITS"), track.pt() * 2.0, track.itsChi2NCl());

          if (track.hasTOF()) {
            if (track.hasTRD() && (lastRequiredTrdCluster > 0)) {
              int lastLayer = 0;
              for (int l = 7; l >= 0; l--) {
                if (track.trdPattern() & (1 << l)) {
                  lastLayer = l;
                  break;
                }
              }
              if (lastLayer < lastRequiredTrdCluster)
                continue;
            }
            Float_t TOFmass2 = ((track.mass()) * (track.mass()));
            Float_t beta = track.beta();
            aHelium3_reg.fill(HIST("histTOFm2"), momentum_Z_equal_2, TOFmass2);
            aHelium3_reg.fill(HIST("histTofSignalData"), momentum_Z_equal_2, beta);
            aHelium3_reg.fill(HIST("histTofNsigmaData"), momentum_Z_equal_2, track.tofNSigmaHe());
          }
        }

        if (track.hasTOF()) {
          if (track.hasTRD() && (lastRequiredTrdCluster > 0)) {
            int lastLayer = 0;
            for (int l = 7; l >= 0; l--) {
              if (track.trdPattern() & (1 << l)) {
                lastLayer = l;
                break;
              }
            }
            if (lastLayer < lastRequiredTrdCluster)
              continue;
          }
          spectra_reg.fill(HIST("histTofSignalData"), momentum_Z_equal_2 * track.sign(), track.beta());
        }
      }

      //**************   Helium-4   *******************

      if (nSigmaHe4 > nsigmacutLow && nSigmaHe4 < nsigmacutHigh) {

        if (track.sign() > 0) {
          Helium4_reg.fill(HIST("histDcaVsPtData"), track.pt() * 2.0, track.dcaXY());
          Helium4_reg.fill(HIST("histDcaZVsPtData"), track.pt() * 2.0, track.dcaZ());
          Helium4_reg.fill(HIST("histTpcSignalData"), momentum_Z_equal_2, track.tpcSignal());
          Helium4_reg.fill(HIST("histNClusterTPC"), track.pt() * 2.0, track.tpcNClsFound());
          Helium4_reg.fill(HIST("histNClusterITS"), track.pt() * 2.0, track.itsNCls());
          Helium4_reg.fill(HIST("histNClusterITSib"), track.pt() * 2.0, track.itsNClsInnerBarrel());
          Helium4_reg.fill(HIST("histChi2TPC"), track.pt() * 2.0, track.tpcChi2NCl());
          Helium4_reg.fill(HIST("histChi2ITS"), track.pt() * 2.0, track.itsChi2NCl());

          if (track.hasTOF()) {
            if (track.hasTRD() && (lastRequiredTrdCluster > 0)) {
              int lastLayer = 0;
              for (int l = 7; l >= 0; l--) {
                if (track.trdPattern() & (1 << l)) {
                  lastLayer = l;
                  break;
                }
              }
              if (lastLayer < lastRequiredTrdCluster)
                continue;
            }
            Float_t TOFmass2 = ((track.mass()) * (track.mass()));
            Float_t beta = track.beta();
            Helium4_reg.fill(HIST("histTOFm2"), momentum_Z_equal_2, TOFmass2);
            Helium4_reg.fill(HIST("histTofSignalData"), momentum_Z_equal_2, beta);
            Helium4_reg.fill(HIST("histTofNsigmaData"), momentum_Z_equal_2, track.tofNSigmaAl());
          }
        }

        if (track.sign() < 0) {
          aHelium4_reg.fill(HIST("histDcaVsPtData"), track.pt() * 2.0, track.dcaXY());
          aHelium4_reg.fill(HIST("histDcaZVsPtData"), track.pt() * 2.0, track.dcaZ());
          aHelium4_reg.fill(HIST("histTpcSignalData"), momentum_Z_equal_2, track.tpcSignal());
          aHelium4_reg.fill(HIST("histNClusterTPC"), track.pt() * 2.0, track.tpcNClsFound());
          aHelium4_reg.fill(HIST("histNClusterITS"), track.pt() * 2.0, track.itsNCls());
          aHelium4_reg.fill(HIST("histNClusterITSib"), track.pt() * 2.0, track.itsNClsInnerBarrel());
          aHelium4_reg.fill(HIST("histChi2TPC"), track.pt() * 2.0, track.tpcChi2NCl());
          aHelium4_reg.fill(HIST("histChi2ITS"), track.pt() * 2.0, track.itsChi2NCl());

          if (track.hasTOF()) {
            if (track.hasTRD() && (lastRequiredTrdCluster > 0)) {
              int lastLayer = 0;
              for (int l = 7; l >= 0; l--) {
                if (track.trdPattern() & (1 << l)) {
                  lastLayer = l;
                  break;
                }
              }
              if (lastLayer < lastRequiredTrdCluster)
                continue;
            }
            Float_t TOFmass2 = ((track.mass()) * (track.mass()));
            Float_t beta = track.beta();
            aHelium4_reg.fill(HIST("histTOFm2"), momentum_Z_equal_2, TOFmass2);
            aHelium4_reg.fill(HIST("histTofSignalData"), momentum_Z_equal_2, beta);
            aHelium4_reg.fill(HIST("histTofNsigmaData"), momentum_Z_equal_2, track.tofNSigmaAl());
          }
        }

        if (track.hasTOF()) {
          if (track.hasTRD() && (lastRequiredTrdCluster > 0)) {
            int lastLayer = 0;
            for (int l = 7; l >= 0; l--) {
              if (track.trdPattern() & (1 << l)) {
                lastLayer = l;
                break;
              }
            }
            if (lastLayer < lastRequiredTrdCluster)
              continue;
          }
          spectra_reg.fill(HIST("histTofSignalData"), momentum_Z_equal_2 * track.sign(), track.beta());
        }
      }
    }
  }

  //****************************************************************************************************

  template <typename CollisionType, typename TracksType>
  void fillCentHistorgrams(const CollisionType& event, const TracksType& tracks)
  {

    if (event_selection_sel8 && event.sel8())
      spectra_reg.fill(HIST("histCentrality"), event.centFT0C());
    if (!event_selection_sel8)
      spectra_reg.fill(HIST("histCentrality"), event.centFT0C());

    for (auto track : tracks) {

      if ((event_selection_sel8 && !event.sel8()) || (enable_Centrality_cut_global && (event.centFT0C() < minCentrality) && (event.centFT0C() > maxCentrality)))
        continue;

      double momentum = 0.0;
      int momentum_getter = 0;
      switch (use_momentum_getter) {
        case 0:
          momentum = track.p();
          break;
        case 1:
          momentum = track.pt();
          break;
        case 2:
          momentum = track.tpcInnerParam();
          break;
        default:
          momentum_getter = -1;
          break;
      }
      if (momentum_getter == -1)
        break;

      spectra_reg.fill(HIST("histEtaWithOverFlow"), track.eta());
      spectra_reg.fill(HIST("histEta_cent"), event.centFT0C(), track.eta());

      if ((event.centFT0C() > minCentrality) && (event.centFT0C() < maxCentrality)) {
        spectra_reg.fill(HIST("histEta"), track.eta());
      }

      float TPCnumberClsFound = track.tpcNClsFound();
      float TPC_nCls_Crossed_Rows = track.tpcNClsCrossedRows();
      float RatioCrossedRowsOverFindableTPC = track.tpcCrossedRowsOverFindableCls();
      float Chi2perClusterTPC = track.tpcChi2NCl();
      float Chi2perClusterITS = track.itsChi2NCl();

      if (TMath::Abs(track.dcaXY()) > maxDCA_XY || TMath::Abs(track.dcaZ()) > maxDCA_Z)
        continue;
      if (TPCnumberClsFound < minTPCnClsFound || TPC_nCls_Crossed_Rows < minNCrossedRowsTPC)
        continue;
      if (RatioCrossedRowsOverFindableTPC < minRatioCrossedRowsTPC || RatioCrossedRowsOverFindableTPC > maxRatioCrossedRowsTPC)
        continue;
      if (Chi2perClusterTPC > maxChi2PerClusterTPC || Chi2perClusterTPC < minChi2PerClusterTPC || Chi2perClusterITS > maxChi2PerClusterITS)
        continue;
      if (!(track.passedTPCRefit()))
        continue;
      if (!(track.passedITSRefit()))
        continue;
      if ((track.itsNClsInnerBarrel()) < minReqClusterITSib || (track.itsNCls()) < minReqClusterITS)
        continue;
      if (track.pt() < p_min || track.pt() > p_max)
        continue;
      if ((requireITS && !(track.hasITS())) || (requireTPC && !(track.hasTPC())))
        continue;
      if (requireGoldenChi2 && !(track.passedGoldenChi2()))
        continue;

      TLorentzVector lorentzVector_pion{};
      TLorentzVector lorentzVector_proton{};
      TLorentzVector lorentzVector_deuteron{};
      TLorentzVector lorentzVector_triton{};
      TLorentzVector lorentzVector_He3{};
      TLorentzVector lorentzVector_He4{};

      lorentzVector_pion.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassPiPlus);
      lorentzVector_proton.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassProton);
      lorentzVector_deuteron.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassDeuteron);
      lorentzVector_triton.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassTriton);
      lorentzVector_He3.SetPtEtaPhiM(track.pt() * 2.0, track.eta(), track.phi(), constants::physics::MassHelium3);
      lorentzVector_He4.SetPtEtaPhiM(track.pt() * 2.0, track.eta(), track.phi(), constants::physics::MassAlpha);

      if (lorentzVector_pion.Rapidity() < yMin || lorentzVector_pion.Rapidity() > yMax ||
          lorentzVector_proton.Rapidity() < yMin || lorentzVector_proton.Rapidity() > yMax ||
          lorentzVector_deuteron.Rapidity() < yMin || lorentzVector_deuteron.Rapidity() > yMax ||
          lorentzVector_triton.Rapidity() < yMin || lorentzVector_triton.Rapidity() > yMax ||
          lorentzVector_He3.Rapidity() < yMin || lorentzVector_He3.Rapidity() > yMax ||
          lorentzVector_He4.Rapidity() < yMin || lorentzVector_He4.Rapidity() > yMax)
        continue;

      if (track.sign() > 0) {

        pion_reg.fill(HIST("histTpcNsigmaData_cent"), momentum, track.tpcNSigmaPi(), event.centFT0C());
        proton_reg.fill(HIST("histTpcNsigmaData_cent"), momentum, track.tpcNSigmaPr(), event.centFT0C());
        deuteron_reg.fill(HIST("histTpcNsigmaData_cent"), momentum, track.tpcNSigmaDe(), event.centFT0C());
        triton_reg.fill(HIST("histTpcNsigmaData_cent"), momentum, track.tpcNSigmaTr(), event.centFT0C());
        Helium3_reg.fill(HIST("histTpcNsigmaData_cent"), momentum * 2.0, track.tpcNSigmaHe(), event.centFT0C());
        Helium4_reg.fill(HIST("histTpcNsigmaData_cent"), momentum * 2.0, track.tpcNSigmaAl(), event.centFT0C());

        if ((event.centFT0C() > minCentrality) && (event.centFT0C() < maxCentrality)) {
          pion_reg.fill(HIST("histTpcNsigmaData_eta"), momentum, track.tpcNSigmaPi(), track.eta());
          proton_reg.fill(HIST("histTpcNsigmaData_eta"), momentum, track.tpcNSigmaPr(), track.eta());
          deuteron_reg.fill(HIST("histTpcNsigmaData_eta"), momentum, track.tpcNSigmaDe(), track.eta());
          triton_reg.fill(HIST("histTpcNsigmaData_eta"), momentum, track.tpcNSigmaTr(), track.eta());
          Helium3_reg.fill(HIST("histTpcNsigmaData_eta"), momentum * 2.0, track.tpcNSigmaHe(), track.eta());
          Helium4_reg.fill(HIST("histTpcNsigmaData_eta"), momentum * 2.0, track.tpcNSigmaAl(), track.eta());
        }

        if (track.hasTOF()) {
          if (track.hasTRD() && (lastRequiredTrdCluster > 0)) {
            int lastLayer = 0;
            for (int l = 7; l >= 0; l--) {
              if (track.trdPattern() & (1 << l)) {
                lastLayer = l;
                break;
              }
            }
            if (lastLayer < lastRequiredTrdCluster)
              continue;
          }
          pion_reg.fill(HIST("histTofNsigmaData_cent"), momentum, track.tofNSigmaPi(), event.centFT0C());
          proton_reg.fill(HIST("histTofNsigmaData_cent"), momentum, track.tofNSigmaPr(), event.centFT0C());
          deuteron_reg.fill(HIST("histTofNsigmaData_cent"), momentum, track.tofNSigmaDe(), event.centFT0C());
          triton_reg.fill(HIST("histTofNsigmaData_cent"), momentum, track.tofNSigmaTr(), event.centFT0C());
          Helium3_reg.fill(HIST("histTofNsigmaData_cent"), momentum * 2.0, track.tofNSigmaHe(), event.centFT0C());
          Helium4_reg.fill(HIST("histTofNsigmaData_cent"), momentum * 2.0, track.tofNSigmaAl(), event.centFT0C());

          pion_reg.fill(HIST("histTofm2_cent"), momentum, track.mass() * track.mass(), event.centFT0C());
          proton_reg.fill(HIST("histTofm2_cent"), momentum, track.mass() * track.mass(), event.centFT0C());
          deuteron_reg.fill(HIST("histTofm2_cent"), momentum, track.mass() * track.mass(), event.centFT0C());
          triton_reg.fill(HIST("histTofm2_cent"), momentum, track.mass() * track.mass(), event.centFT0C());
          Helium3_reg.fill(HIST("histTofm2_cent"), momentum * 2.0, track.mass() * track.mass(), event.centFT0C());
          Helium4_reg.fill(HIST("histTofm2_cent"), momentum * 2.0, track.mass() * track.mass(), event.centFT0C());

          if ((event.centFT0C() > minCentrality) && (event.centFT0C() < maxCentrality)) {
            pion_reg.fill(HIST("histTofm2_eta"), momentum, track.mass() * track.mass(), track.eta());
            pion_reg.fill(HIST("histTofNsigmaData_eta"), momentum, track.tofNSigmaPi(), track.eta());
            proton_reg.fill(HIST("histTofm2_eta"), momentum, track.mass() * track.mass(), track.eta());
            proton_reg.fill(HIST("histTofNsigmaData_eta"), momentum, track.tofNSigmaPr(), track.eta());
            deuteron_reg.fill(HIST("histTofm2_eta"), momentum, track.mass() * track.mass(), track.eta());
            deuteron_reg.fill(HIST("histTofNsigmaData_eta"), momentum, track.tofNSigmaDe(), track.eta());
            triton_reg.fill(HIST("histTofm2_eta"), momentum, track.mass() * track.mass(), track.eta());
            triton_reg.fill(HIST("histTofNsigmaData_eta"), momentum, track.tofNSigmaTr(), track.eta());
            Helium3_reg.fill(HIST("histTofm2_eta"), momentum * 2.0, track.mass() * track.mass(), track.eta());
            Helium3_reg.fill(HIST("histTofNsigmaData_eta"), momentum * 2.0, track.tofNSigmaHe(), track.eta());
            Helium4_reg.fill(HIST("histTofm2_eta"), momentum * 2.0, track.mass() * track.mass(), track.eta());
            Helium4_reg.fill(HIST("histTofNsigmaData_eta"), momentum * 2.0, track.tofNSigmaAl(), track.eta());
          }
        }
      }

      if (track.sign() < 0) {

        apion_reg.fill(HIST("histTpcNsigmaData_cent"), momentum, track.tpcNSigmaPi(), event.centFT0C());
        aproton_reg.fill(HIST("histTpcNsigmaData_cent"), momentum, track.tpcNSigmaPr(), event.centFT0C());
        adeuteron_reg.fill(HIST("histTpcNsigmaData_cent"), momentum, track.tpcNSigmaDe(), event.centFT0C());
        atriton_reg.fill(HIST("histTpcNsigmaData_cent"), momentum, track.tpcNSigmaTr(), event.centFT0C());
        aHelium3_reg.fill(HIST("histTpcNsigmaData_cent"), momentum * 2.0, track.tpcNSigmaHe(), event.centFT0C());
        aHelium4_reg.fill(HIST("histTpcNsigmaData_cent"), momentum * 2.0, track.tpcNSigmaAl(), event.centFT0C());

        if ((event.centFT0C() > minCentrality) && (event.centFT0C() < maxCentrality)) {
          apion_reg.fill(HIST("histTpcNsigmaData_eta"), momentum, track.tpcNSigmaPi(), track.eta());
          aproton_reg.fill(HIST("histTpcNsigmaData_eta"), momentum, track.tpcNSigmaPr(), track.eta());
          adeuteron_reg.fill(HIST("histTpcNsigmaData_eta"), momentum, track.tpcNSigmaDe(), track.eta());
          atriton_reg.fill(HIST("histTpcNsigmaData_eta"), momentum, track.tpcNSigmaTr(), track.eta());
          aHelium3_reg.fill(HIST("histTpcNsigmaData_eta"), momentum * 2.0, track.tpcNSigmaHe(), track.eta());
          aHelium4_reg.fill(HIST("histTpcNsigmaData_eta"), momentum * 2.0, track.tpcNSigmaAl(), track.eta());
        }

        if (track.hasTOF()) {
          if (track.hasTRD() && (lastRequiredTrdCluster > 0)) {
            int lastLayer = 0;
            for (int l = 7; l >= 0; l--) {
              if (track.trdPattern() & (1 << l)) {
                lastLayer = l;
                break;
              }
            }
            if (lastLayer < lastRequiredTrdCluster)
              continue;
          }
          apion_reg.fill(HIST("histTofNsigmaData_cent"), momentum, track.tofNSigmaPi(), event.centFT0C());
          aproton_reg.fill(HIST("histTofNsigmaData_cent"), momentum, track.tofNSigmaPr(), event.centFT0C());
          adeuteron_reg.fill(HIST("histTofNsigmaData_cent"), momentum, track.tofNSigmaDe(), event.centFT0C());
          atriton_reg.fill(HIST("histTofNsigmaData_cent"), momentum, track.tofNSigmaTr(), event.centFT0C());
          aHelium3_reg.fill(HIST("histTofNsigmaData_cent"), momentum * 2.0, track.tofNSigmaHe(), event.centFT0C());
          aHelium4_reg.fill(HIST("histTofNsigmaData_cent"), momentum * 2.0, track.tofNSigmaAl(), event.centFT0C());

          apion_reg.fill(HIST("histTofm2_cent"), momentum, track.mass() * track.mass(), event.centFT0C());
          aproton_reg.fill(HIST("histTofm2_cent"), momentum, track.mass() * track.mass(), event.centFT0C());
          adeuteron_reg.fill(HIST("histTofm2_cent"), momentum, track.mass() * track.mass(), event.centFT0C());
          atriton_reg.fill(HIST("histTofm2_cent"), momentum, track.mass() * track.mass(), event.centFT0C());
          aHelium3_reg.fill(HIST("histTofm2_cent"), momentum * 2.0, track.mass() * track.mass(), event.centFT0C());
          aHelium4_reg.fill(HIST("histTofm2_cent"), momentum * 2.0, track.mass() * track.mass(), event.centFT0C());

          if ((event.centFT0C() > minCentrality) && (event.centFT0C() < maxCentrality)) {
            apion_reg.fill(HIST("histTofm2_eta"), momentum, track.mass() * track.mass(), track.eta());
            apion_reg.fill(HIST("histTofNsigmaData_eta"), momentum, track.tofNSigmaPi(), track.eta());
            aproton_reg.fill(HIST("histTofm2_eta"), momentum, track.mass() * track.mass(), track.eta());
            aproton_reg.fill(HIST("histTofNsigmaData_eta"), momentum, track.tofNSigmaPr(), track.eta());
            adeuteron_reg.fill(HIST("histTofm2_eta"), momentum, track.mass() * track.mass(), track.eta());
            adeuteron_reg.fill(HIST("histTofNsigmaData_eta"), momentum, track.tofNSigmaDe(), track.eta());
            atriton_reg.fill(HIST("histTofm2_eta"), momentum, track.mass() * track.mass(), track.eta());
            atriton_reg.fill(HIST("histTofNsigmaData_eta"), momentum, track.tofNSigmaTr(), track.eta());
            aHelium3_reg.fill(HIST("histTofm2_eta"), momentum * 2.0, track.mass() * track.mass(), track.eta());
            aHelium3_reg.fill(HIST("histTofNsigmaData_eta"), momentum * 2.0, track.tofNSigmaHe(), track.eta());
            aHelium4_reg.fill(HIST("histTofm2_eta"), momentum * 2.0, track.mass() * track.mass(), track.eta());
            aHelium4_reg.fill(HIST("histTofNsigmaData_eta"), momentum * 2.0, track.tofNSigmaAl(), track.eta());
          }
        }
      }
    }
  }

  //****************************************************************************************************

  Filter collisionFilter = (nabs(aod::collision::posZ) < cfgCutVertex);
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta && requireGlobalTrackWoDCAInFilter());

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>;

  using EventCandidatesCent = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>>;

  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCLfFullPi, aod::pidTOFFullPi, aod::pidTPCLfFullKa, aod::pidTOFFullKa, aod::pidTPCLfFullPr, aod::pidTOFFullPr, aod::pidTPCLfFullDe, aod::pidTOFFullDe, aod::pidTPCLfFullTr, aod::pidTOFFullTr, aod::pidTPCLfFullHe, aod::pidTOFFullHe, aod::pidTPCLfFullAl, aod::pidTOFFullAl, aod::TrackSelection, aod::TrackSelectionExtension, aod::TOFSignal, aod::pidTOFmass, aod::pidTOFbeta>>;

  void processData(EventCandidates::iterator const& event, TrackCandidates const& tracks)
  {
    fillHistograms(event, tracks);
  }
  PROCESS_SWITCH(NucleiHistTask, processData, "process data", true);

  void processDataCent(EventCandidatesCent::iterator const& event, TrackCandidates const& tracks)
  {
    fillHistograms(event, tracks);
    fillCentHistorgrams(event, tracks);
  }
  PROCESS_SWITCH(NucleiHistTask, processDataCent, "process data with centralities", false);

  void processMC(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0Cs>::iterator const& collisions, soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCLfFullPi, aod::pidTOFFullPi, aod::pidTPCLfFullKa, aod::pidTOFFullKa, aod::pidTPCLfFullPr, aod::pidTOFFullPr, aod::pidTPCLfFullDe, aod::pidTOFFullDe, aod::pidTPCLfFullTr, aod::pidTOFFullTr, aod::pidTPCLfFullHe, aod::pidTOFFullHe, aod::pidTPCLfFullAl, aod::pidTOFFullAl, aod::McTrackLabels, aod::TrackSelection, aod::TrackSelectionExtension, aod::TOFSignal, aod::pidTOFmass, aod::pidTOFbeta>> const& tracks,
                 aod::McParticles& /*mcParticles*/, aod::McCollisions const& /*mcCollisions*/)
  {

    if (event_selection_MC_sel8 && !collisions.sel8())
      return;
    MC_recon_reg.fill(HIST("histRecVtxMC"), collisions.posZ());
    MC_recon_reg.fill(HIST("histCentrality"), collisions.centFT0C());

    for (auto& track : tracks) {
      histTrackcuts_MC->AddBinContent(1);
      const auto particle = track.mcParticle();
      if (!particle.isPhysicalPrimary())
        continue;
      histTrackcuts_MC->AddBinContent(2);

      int pdgbin = 0;
      switch (particle.pdgCode()) {
        case +211:
          histPDG->AddBinContent(1);
          pdgbin = 0;
          break;
        case -211:
          histPDG->AddBinContent(2);
          pdgbin = 1;
          break;
        case +321:
          histPDG->AddBinContent(3);
          pdgbin = 2;
          break;
        case -321:
          histPDG->AddBinContent(4);
          pdgbin = 3;
          break;
        case +2212:
          histPDG->AddBinContent(5);
          pdgbin = 4;
          break;
        case -2212:
          histPDG->AddBinContent(6);
          pdgbin = 5;
          break;
        case +1000010020:
          histPDG->AddBinContent(7);
          pdgbin = 6;
          break;
        case -1000010020:
          histPDG->AddBinContent(8);
          pdgbin = 7;
          break;
        case +1000010030:
          histPDG->AddBinContent(9);
          pdgbin = 8;
          break;
        case -1000010030:
          histPDG->AddBinContent(10);
          pdgbin = 9;
          break;
        case +1000020030:
          histPDG->AddBinContent(11);
          pdgbin = 10;
          break;
        case -1000020030:
          histPDG->AddBinContent(12);
          pdgbin = 11;
          break;
        case +1000020040:
          histPDG->AddBinContent(13);
          pdgbin = 12;
          break;
        case -1000020040:
          histPDG->AddBinContent(14);
          pdgbin = 13;
          break;
        default:
          pdgbin = -1;
          break;
      }

      double momentum = 0.0;
      int momentum_getter = 0;
      switch (use_momentum_getter) {
        case 0:
          momentum = track.p();
          break;
        case 1:
          momentum = track.pt();
          break;
        case 2:
          momentum = track.tpcInnerParam();
          break;
        default:
          momentum_getter = -1;
          break;
      }
      if (momentum_getter == -1)
        break;

      TLorentzVector lorentzVector_pion{};
      TLorentzVector lorentzVector_kaon{};
      TLorentzVector lorentzVector_proton{};
      TLorentzVector lorentzVector_deuteron{};
      TLorentzVector lorentzVector_triton{};
      TLorentzVector lorentzVector_He3{};
      TLorentzVector lorentzVector_He4{};

      lorentzVector_pion.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassPiPlus);
      lorentzVector_kaon.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassKPlus);
      lorentzVector_proton.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassProton);
      lorentzVector_deuteron.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassDeuteron);
      lorentzVector_triton.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassTriton);
      lorentzVector_He3.SetPtEtaPhiM(track.pt() * 2.0, track.eta(), track.phi(), constants::physics::MassHelium3);
      lorentzVector_He4.SetPtEtaPhiM(track.pt() * 2.0, track.eta(), track.phi(), constants::physics::MassAlpha);

      if (lorentzVector_pion.Rapidity() < yMin || lorentzVector_pion.Rapidity() > yMax ||
          lorentzVector_kaon.Rapidity() < yMin || lorentzVector_kaon.Rapidity() > yMax ||
          lorentzVector_proton.Rapidity() < yMin || lorentzVector_proton.Rapidity() > yMax ||
          lorentzVector_deuteron.Rapidity() < yMin || lorentzVector_deuteron.Rapidity() > yMax ||
          lorentzVector_triton.Rapidity() < yMin || lorentzVector_triton.Rapidity() > yMax ||
          lorentzVector_He3.Rapidity() < yMin || lorentzVector_He3.Rapidity() > yMax ||
          lorentzVector_He4.Rapidity() < yMin || lorentzVector_He4.Rapidity() > yMax)
        continue;
      histTrackcuts_MC->AddBinContent(3);

      MC_recon_reg.fill(HIST("histEta"), track.eta(), pdgbin);

      if ((particle.pdgCode() == 1000020030) || (particle.pdgCode() == -1000020030) || (particle.pdgCode() == 1000020040) || (particle.pdgCode() == -1000020040)) {
        MC_recon_reg.fill(HIST("histPt"), track.pt() * 2.0, pdgbin);
        MC_recon_reg.fill(HIST("histDCA"), track.pt() * 2.0, track.dcaXY(), pdgbin);
        MC_recon_reg.fill(HIST("histDCAz"), track.pt() * 2.0, track.dcaZ(), pdgbin);
        MC_recon_reg.fill(HIST("histTpcSignalData"), momentum * 2.0 * track.sign(), track.tpcSignal(), pdgbin);
        MC_recon_reg.fill(HIST("histTpcSignalData_all_species"), momentum * 2.0 * track.sign(), track.tpcSignal());
        MC_recon_reg.fill(HIST("histNClusterTPC"), track.pt() * 2.0, track.tpcNClsCrossedRows(), pdgbin);
        MC_recon_reg.fill(HIST("histNClusterITS"), track.pt() * 2.0, track.itsNCls(), pdgbin);
        MC_recon_reg.fill(HIST("histNClusterITSib"), track.pt() * 2.0, track.itsNClsInnerBarrel(), pdgbin);
        MC_recon_reg.fill(HIST("histChi2TPC"), track.pt() * 2.0, track.tpcChi2NCl(), pdgbin);
        MC_recon_reg.fill(HIST("histChi2ITS"), track.pt() * 2.0, track.itsChi2NCl(), pdgbin);
        MC_recon_reg.fill(HIST("histChi2TOF"), track.pt() * 2.0, track.tofChi2(), pdgbin);
      } else {
        MC_recon_reg.fill(HIST("histPt"), track.pt(), pdgbin);
        MC_recon_reg.fill(HIST("histDCA"), track.pt(), track.dcaXY(), pdgbin);
        MC_recon_reg.fill(HIST("histDCAz"), track.pt(), track.dcaZ(), pdgbin);
        MC_recon_reg.fill(HIST("histTpcSignalData"), momentum * track.sign(), track.tpcSignal(), pdgbin);
        MC_recon_reg.fill(HIST("histTpcSignalData_all_species"), momentum * track.sign(), track.tpcSignal());
        MC_recon_reg.fill(HIST("histNClusterTPC"), track.pt(), track.tpcNClsCrossedRows(), pdgbin);
        MC_recon_reg.fill(HIST("histNClusterITS"), track.pt(), track.itsNCls(), pdgbin);
        MC_recon_reg.fill(HIST("histNClusterITSib"), track.pt(), track.itsNClsInnerBarrel(), pdgbin);
        MC_recon_reg.fill(HIST("histChi2TPC"), track.pt(), track.tpcChi2NCl(), pdgbin);
        MC_recon_reg.fill(HIST("histChi2ITS"), track.pt(), track.itsChi2NCl(), pdgbin);
        MC_recon_reg.fill(HIST("histChi2TOF"), track.pt(), track.tofChi2(), pdgbin);
      }

      float TPCnumberClsFound = track.tpcNClsFound();
      float TPC_nCls_Crossed_Rows = track.tpcNClsCrossedRows();
      float RatioCrossedRowsOverFindableTPC = track.tpcCrossedRowsOverFindableCls();
      float Chi2perClusterTPC = track.tpcChi2NCl();
      float Chi2perClusterITS = track.itsChi2NCl();

      if (TMath::Abs(track.dcaXY()) > maxDCA_XY || TMath::Abs(track.dcaZ()) > maxDCA_Z)
        continue;
      histTrackcuts_MC->AddBinContent(4);
      if (TPCnumberClsFound < minTPCnClsFound || TPC_nCls_Crossed_Rows < minNCrossedRowsTPC)
        continue;
      histTrackcuts_MC->AddBinContent(5);
      if (RatioCrossedRowsOverFindableTPC < minRatioCrossedRowsTPC || RatioCrossedRowsOverFindableTPC > maxRatioCrossedRowsTPC)
        continue;
      histTrackcuts_MC->AddBinContent(6);
      if (Chi2perClusterTPC > maxChi2PerClusterTPC || Chi2perClusterTPC < minChi2PerClusterTPC || Chi2perClusterITS > maxChi2PerClusterITS)
        continue;
      histTrackcuts_MC->AddBinContent(7);
      if (!(track.passedTPCRefit()))
        continue;
      histTrackcuts_MC->AddBinContent(8);
      if (!(track.passedITSRefit()))
        continue;
      histTrackcuts_MC->AddBinContent(9);
      if ((track.itsNClsInnerBarrel()) < minReqClusterITSib || (track.itsNCls()) < minReqClusterITS)
        continue;
      histTrackcuts_MC->AddBinContent(10);
      if (track.pt() < p_min || track.pt() > p_max)
        continue;
      histTrackcuts_MC->AddBinContent(11);

      if ((requireITS && !(track.hasITS())) || (requireTPC && !(track.hasTPC())))
        continue;
      histTrackcuts_MC->AddBinContent(12);
      if (requireGoldenChi2 && !(track.passedGoldenChi2()))
        continue;
      histTrackcuts_MC->AddBinContent(13);

      float nSigmaPion = track.tpcNSigmaPi();
      float nSigmaProton = track.tpcNSigmaPr();
      float nSigmaDeuteron = track.tpcNSigmaDe();
      float nSigmaTriton = track.tpcNSigmaTr();
      float nSigmaHe3 = track.tpcNSigmaHe();
      float nSigmaHe4 = track.tpcNSigmaAl();

      if (track.sign() > 0) {
        MC_recon_reg.fill(HIST("histTpcNsigmaDataPi"), momentum, nSigmaPion);
        MC_recon_reg.fill(HIST("histTpcNsigmaDataPr"), momentum, nSigmaProton);
        MC_recon_reg.fill(HIST("histTpcNsigmaDataDe"), momentum, nSigmaDeuteron);
        MC_recon_reg.fill(HIST("histTpcNsigmaDataTr"), momentum, nSigmaTriton);
        MC_recon_reg.fill(HIST("histTpcNsigmaDataHe"), momentum * 2.0, nSigmaHe3);
        MC_recon_reg.fill(HIST("histTpcNsigmaDataAl"), momentum * 2.0, nSigmaHe4);
      }
      if (track.sign() < 0) {
        MC_recon_reg.fill(HIST("histTpcNsigmaDataaPi"), momentum, nSigmaPion);
        MC_recon_reg.fill(HIST("histTpcNsigmaDataaPr"), momentum, nSigmaProton);
        MC_recon_reg.fill(HIST("histTpcNsigmaDataaDe"), momentum, nSigmaDeuteron);
        MC_recon_reg.fill(HIST("histTpcNsigmaDataaTr"), momentum, nSigmaTriton);
        MC_recon_reg.fill(HIST("histTpcNsigmaDataaHe"), momentum * 2.0, nSigmaHe3);
        MC_recon_reg.fill(HIST("histTpcNsigmaDataaAl"), momentum * 2.0, nSigmaHe4);
      }

      if (track.hasTOF()) {
        Float_t TOFmass2 = ((track.mass()) * (track.mass()));

        MC_recon_reg.fill(HIST("histTOFm2"), momentum, TOFmass2, pdgbin);
        MC_recon_reg.fill(HIST("histTofSignalData"), momentum * track.sign(), track.beta(), pdgbin);
        MC_recon_reg.fill(HIST("histTofSignalData_all_species"), momentum * track.sign(), track.beta());

        if (track.sign() > 0) {
          if (nSigmaPion > nsigmacutLow && nSigmaPion < nsigmacutHigh)
            MC_recon_reg.fill(HIST("histTofNsigmaDataPi"), momentum, track.tofNSigmaPi());
          if (nSigmaProton > nsigmacutLow && nSigmaProton < nsigmacutHigh)
            MC_recon_reg.fill(HIST("histTofNsigmaDataPr"), momentum, track.tofNSigmaPr());
          if (nSigmaDeuteron > nsigmacutLow && nSigmaDeuteron < nsigmacutHigh)
            MC_recon_reg.fill(HIST("histTofNsigmaDataDe"), momentum, track.tofNSigmaDe());
          if (nSigmaTriton > nsigmacutLow && nSigmaTriton < nsigmacutHigh)
            MC_recon_reg.fill(HIST("histTofNsigmaDataTr"), momentum, track.tofNSigmaTr());
          if (nSigmaHe3 > nsigmacutLow && nSigmaHe3 < nsigmacutHigh)
            MC_recon_reg.fill(HIST("histTofNsigmaDataHe"), momentum * 2.0, track.tofNSigmaHe());
          if (nSigmaHe4 > nsigmacutLow && nSigmaHe4 < nsigmacutHigh)
            MC_recon_reg.fill(HIST("histTofNsigmaDataAl"), momentum * 2.0, track.tofNSigmaAl());
        }
        if (track.sign() < 0) {
          if (nSigmaPion > nsigmacutLow && nSigmaPion < nsigmacutHigh)
            MC_recon_reg.fill(HIST("histTofNsigmaDataaPi"), momentum, track.tofNSigmaPi());
          if (nSigmaProton > nsigmacutLow && nSigmaProton < nsigmacutHigh)
            MC_recon_reg.fill(HIST("histTofNsigmaDataaPr"), momentum, track.tofNSigmaPr());
          if (nSigmaDeuteron > nsigmacutLow && nSigmaDeuteron < nsigmacutHigh)
            MC_recon_reg.fill(HIST("histTofNsigmaDataaDe"), momentum, track.tofNSigmaDe());
          if (nSigmaTriton > nsigmacutLow && nSigmaTriton < nsigmacutHigh)
            MC_recon_reg.fill(HIST("histTofNsigmaDataaTr"), momentum, track.tofNSigmaTr());
          if (nSigmaHe3 > nsigmacutLow && nSigmaHe3 < nsigmacutHigh)
            MC_recon_reg.fill(HIST("histTofNsigmaDataaHe"), momentum * 2.0, track.tofNSigmaHe());
          if (nSigmaHe4 > nsigmacutLow && nSigmaHe4 < nsigmacutHigh)
            MC_recon_reg.fill(HIST("histTofNsigmaDataaAl"), momentum * 2.0, track.tofNSigmaAl());
        }
      }
    }
  }
  PROCESS_SWITCH(NucleiHistTask, processMC, "process MC", false);
};

//****************************************************************************************************

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<NucleiHistTask>(cfgc, TaskName{"nuclei-hist"})};
}
