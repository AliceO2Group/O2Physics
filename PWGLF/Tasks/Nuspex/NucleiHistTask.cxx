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
#include <TF1.h>
#include <string>

#include "ReconstructionDataFormats/Track.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/McCollisionExtra.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Framework/HistogramRegistry.h"
#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "TPDGCode.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/Core/TrackSelection.h"
#include "Framework/StaticFor.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "PWGLF/DataModel/spectraTOF.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "PWGLF/Utils/inelGt.h"
#include "PWGLF/DataModel/mcCentrality.h"
#include "Common/Core/RecoDecay.h"

using namespace o2;
using namespace o2::track;
using namespace o2::framework;
using namespace o2::framework::expressions;

std::vector<float> ptBinning = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 3.6, 4., 5., 6., 8., 10., 12., 14.};

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
  OutputObj<TH1I> histTrackcuts_data_spectra{TH1I("histTrackcuts_data_spectra", "Entires;Track cut", 18, 0, 18)};
  OutputObj<TH1I> histTrackcuts_data_particle{TH1I("histTrackcuts_data_particle", "Entires;Track cut", 18, 0, 18)};

  // MC
  HistogramRegistry MC_recon_reg{"MC_particles_reco", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry MC_gen_reg{"MC_particles_gen", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry MC_DCA{"MC_DCA", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  OutputObj<TH1I> histPDG_reco{TH1I("PDG reconstructed", "PDG;PDG code", 18, 0.0, 18)};
  OutputObj<TH1I> histPDG_gen{TH1I("PDG generated", "PDG;PDG code", 18, 0.0, 18)};
  OutputObj<TH1I> histTrackcuts_MC{TH1I("histTrackcuts_MC", "Entires;Track cut", 18, 0, 18)};

  void init(o2::framework::InitContext&)
  {

    // +++++++++++++++++++++ Spectra ++++++++++++++++++++++++
    histTrackcuts_data_spectra->GetXaxis()->SetBinLabel(1, "Events read");
    histTrackcuts_data_spectra->GetXaxis()->SetBinLabel(2, "Ev. sel. passed");
    histTrackcuts_data_spectra->GetXaxis()->SetBinLabel(3, "DCA cut passed");
    histTrackcuts_data_spectra->GetXaxis()->SetBinLabel(4, "TPCnCls cut passed");
    histTrackcuts_data_spectra->GetXaxis()->SetBinLabel(5, "TPCCrossedRowsOverFindable cut passed");
    histTrackcuts_data_spectra->GetXaxis()->SetBinLabel(6, "Chi2 cut passed");
    histTrackcuts_data_spectra->GetXaxis()->SetBinLabel(7, "Passed TPC refit cut");
    histTrackcuts_data_spectra->GetXaxis()->SetBinLabel(8, "Passed ITS refit cut");
    histTrackcuts_data_spectra->GetXaxis()->SetBinLabel(9, "ITSnCls cut passed");
    histTrackcuts_data_spectra->GetXaxis()->SetBinLabel(10, "track.pt() cut passed");
    histTrackcuts_data_spectra->GetXaxis()->SetBinLabel(11, "hasITS & hasTPC cut passed");
    histTrackcuts_data_spectra->GetXaxis()->SetBinLabel(12, "GoldenChi2 cut passed");

    // +++++++++++++++++++++ Data particle ++++++++++++++++++++++++
    histTrackcuts_data_particle->GetXaxis()->SetBinLabel(1, "Events read");
    histTrackcuts_data_particle->GetXaxis()->SetBinLabel(2, "Ev. sel. passed");
    histTrackcuts_data_particle->GetXaxis()->SetBinLabel(3, "Rap. cut passed");
    histTrackcuts_data_particle->GetXaxis()->SetBinLabel(4, "TPCnCls cut passed");
    histTrackcuts_data_particle->GetXaxis()->SetBinLabel(5, "TPCCrossedRowsOverFindable cut passed");
    histTrackcuts_data_particle->GetXaxis()->SetBinLabel(6, "Chi2 cut passed");
    histTrackcuts_data_particle->GetXaxis()->SetBinLabel(7, "Passed TPC refit cut");
    histTrackcuts_data_particle->GetXaxis()->SetBinLabel(8, "Passed ITS refit cut");
    histTrackcuts_data_particle->GetXaxis()->SetBinLabel(9, "ITSnCls cut passed");
    histTrackcuts_data_particle->GetXaxis()->SetBinLabel(10, "track.pt() cut passed");
    histTrackcuts_data_particle->GetXaxis()->SetBinLabel(11, "hasITS & hasTPC cut passed");
    histTrackcuts_data_particle->GetXaxis()->SetBinLabel(12, "GoldenChi2 cut passed");
    histTrackcuts_data_particle->GetXaxis()->SetBinLabel(13, "DCA cut passed");

    // +++++++++++++++++++ reconstructed MC ++++++++++++++++++++++
    histTrackcuts_MC->GetXaxis()->SetBinLabel(1, "Events read");
    histTrackcuts_MC->GetXaxis()->SetBinLabel(2, "Is Physical Primary");
    histTrackcuts_MC->GetXaxis()->SetBinLabel(3, "Rap. cut passed");
    histTrackcuts_MC->GetXaxis()->SetBinLabel(4, "DCA cut passed");
    histTrackcuts_MC->GetXaxis()->SetBinLabel(5, "TPCnCls cut passed");
    histTrackcuts_MC->GetXaxis()->SetBinLabel(6, "TPCCrossedRowsOverFindable cut passed");
    histTrackcuts_MC->GetXaxis()->SetBinLabel(7, "Chi2 cut passed");
    histTrackcuts_MC->GetXaxis()->SetBinLabel(8, "Passed TPC refit cut");
    histTrackcuts_MC->GetXaxis()->SetBinLabel(9, "Passed ITS refit cut");
    histTrackcuts_MC->GetXaxis()->SetBinLabel(10, "ITSnCls cut passed");
    histTrackcuts_MC->GetXaxis()->SetBinLabel(11, "track.pt() cut passed");
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
    pion_reg.add("histDcaVsPtData_after_cut", "dcaXY vs Pt (#pi^{+}) after cut", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    pion_reg.add("histDcaZVsPtData_after_cut", "dcaZ vs Pt (#pi^{+}) after cut", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
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
    apion_reg.add("histDcaVsPtData_after_cut", "dcaXY vs Pt (#pi^{-}) after cut", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    apion_reg.add("histDcaZVsPtData_after_cut", "dcaZ vs Pt (#pi^{-}) after cut", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
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
    proton_reg.add("histDcaVsPtData_after_cut", "dcaXY vs Pt (p) after cut", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    proton_reg.add("histDcaZVsPtData_after_cut", "dcaZ vs Pt (p) after cut", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
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
    aproton_reg.add("histDcaVsPtData_after_cut", "dcaXY vs Pt (#bar{p}) after cut", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    aproton_reg.add("histDcaZVsPtData_after_cut", "dcaZ vs Pt (#bar{p}) after cut", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
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
    deuteron_reg.add("histDcaVsPtData_after_cut", "dcaXY vs Pt (d) after cut", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    deuteron_reg.add("histDcaZVsPtData_after_cut", "dcaZ vs Pt (d) after cut", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
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
    adeuteron_reg.add("histDcaVsPtData_after_cut", "dcaXY vs Pt (#bar{d}) after cut", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    adeuteron_reg.add("histDcaZVsPtData_after_cut", "dcaZ vs Pt (#bar{d}) after cut", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
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
    triton_reg.add("histDcaVsPtData_after_cut", "dcaXY vs Pt (t) after cut", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    triton_reg.add("histDcaZVsPtData_after_cut", "dcaZ vs Pt (t) after cut", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
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
    atriton_reg.add("histDcaVsPtData_after_cut", "dcaXY vs Pt (#bar{t}) after cut", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    atriton_reg.add("histDcaZVsPtData_after_cut", "dcaZ vs Pt (#bar{t}) after cut", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
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
    Helium3_reg.add("histDcaVsPtData_after_cut", "dcaXY vs Pt (^{3}He) after cut", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    Helium3_reg.add("histDcaZVsPtData_after_cut", "dcaZ vs Pt (^{3}He) after cut", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
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
    aHelium3_reg.add("histDcaVsPtData_after_cut", "dcaXY vs Pt (^{3}#bar{He}) after cut", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    aHelium3_reg.add("histDcaZVsPtData_after_cut", "dcaZ vs Pt (^{3}#bar{He}) after cut", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
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
    Helium4_reg.add("histDcaVsPtData_after_cut", "dcaXY vs Pt (^{4}He) after cut", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    Helium4_reg.add("histDcaZVsPtData_after_cut", "dcaZ vs Pt (^{4}He) after cut", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
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
    aHelium4_reg.add("histDcaVsPtData_after_cut", "dcaXY vs Pt (^{4}#bar{He}) after cut", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    aHelium4_reg.add("histDcaZVsPtData_after_cut", "dcaZ vs Pt (^{4}#bar{He}) after cut", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
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

    // MC generated
    MC_gen_reg.add("histGenVtxMC", "MC generated vertex z position", HistType::kTH1F, {{400, -40., +40., "z position (cm)"}});
    MC_gen_reg.add("histCentrality", "Centrality", HistType::kTH1F, {centralityAxis});
    MC_gen_reg.add("histEta", "#eta", HistType::kTH2F, {{102, -2.01, 2.01}, PDGBINNING});
    MC_gen_reg.add("histPt", "p_{t}", HistType::kTH2F, {ptAxis, PDGBINNING});

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

    MC_DCA.add("histEta", "#eta", HistType::kTH2F, {{204, -2.01, 2.01}, PDGBINNING});
    MC_DCA.add("histDCA_prim", "DCA xy", HistType::kTH3F, {ptAxis, {250, -0.5, 0.5, "dca"}, PDGBINNING});
    MC_DCA.add("histDCAz_prim", "DCA z", HistType::kTH3F, {ptAxis, {1000, -2.0, 2.0, "dca"}, PDGBINNING});
    MC_DCA.add("histDCA_weak", "DCA xy", HistType::kTH3F, {ptAxis, {250, -0.5, 0.5, "dca"}, PDGBINNING});
    MC_DCA.add("histDCAz_weak", "DCA z", HistType::kTH3F, {ptAxis, {1000, -2.0, 2.0, "dca"}, PDGBINNING});
    MC_DCA.add("histDCA_mat", "DCA xy", HistType::kTH3F, {ptAxis, {250, -0.5, 0.5, "dca"}, PDGBINNING});
    MC_DCA.add("histDCAz_mat", "DCA z", HistType::kTH3F, {ptAxis, {1000, -2.0, 2.0, "dca"}, PDGBINNING});
  }

  // particle configurables
  Configurable<bool> do_pion{"do_pion", false, "fill pion histograms"};
  Configurable<bool> do_proton{"do_proton", false, "fill proton histograms"};
  Configurable<bool> do_deuteron{"do_deuteron", false, "fill deuteron histograms"};
  Configurable<bool> do_triton{"do_triton", false, "fill triton histograms"};
  Configurable<bool> do_He3{"do_He3", false, "fill Helium-3 histograms"};
  Configurable<bool> do_He4{"do_He4", false, "fill Helium-4 histograms"};

  // Configurables
  Configurable<float> yMin{"yMin", -0.5, "Maximum rapidity"};
  Configurable<float> yMax{"yMax", 0.5, "Minimum rapidity"};
  Configurable<float> p_min{"p_min", 0.1f, "min track.pt()"};
  Configurable<float> p_max{"p_max", 1e+10f, "max track.pt()"};
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.9f, "Eta range for tracks"};
  Configurable<float> nsigmacutLow{"nsigmacutLow", -3.0, "Value of the Nsigma cut"};
  Configurable<float> nsigmacutHigh{"nsigmacutHigh", +3.0, "Value of the Nsigma cut"};
  Configurable<float> minCentrality{"minCentrality", 0.0, "min Centrality used"};
  Configurable<float> maxCentrality{"maxCentrality", 80.0, "max Centrality used"};
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
  // Configurable<float> maxDCA_XY{"maxDCA_XY", 0.5f, "max DCA to vertex xy"};
  Configurable<float> maxDcaXYFactor{"maxDcaXYFactor", 0.5f, "DCA xy factor"};
  Configurable<float> maxDCA_Z{"maxDCA_Z", 2.0f, "max DCA to vertex z"};
  Configurable<int> lastRequiredTrdCluster{"lastRequiredTrdCluster", -1, "Last cluster to required in TRD for track selection. -1 does not require any TRD cluster"};
  Configurable<bool> requireGoldenChi2{"requireGoldenChi2", false, "Enable the requirement of GoldenChi2"};
  Configurable<bool> event_selection_sel8{"event_selection_sel8", true, "Enable sel8 event selection"};
  Configurable<bool> event_selection_MC_sel8{"event_selection_MC_sel8", true, "Enable sel8 event selection in MC processing"};
  Configurable<bool> require_PhysicalPrimary_MC_reco{"require_PhysicalPrimary_MC_reco", true, "Enable PhysicalPrimary selection in reconstructed MC processing"};
  Configurable<bool> require_PhysicalPrimary_MC_gen{"require_PhysicalPrimary_MC_gen", true, "Enable PhysicalPrimary selection in generated MC processing"};
  Configurable<bool> removeITSROFrameBorder{"removeITSROFrameBorder", false, "Remove TF border"};
  Configurable<bool> removeNoSameBunchPileup{"removeNoSameBunchPileup", false, "Remove TF border"};
  Configurable<bool> requireIsGoodZvtxFT0vsPV{"requireIsGoodZvtxFT0vsPV", false, "Remove TF border"};
  Configurable<bool> requireIsVertexITSTPC{"requireIsVertexITSTPC", false, "Remove TF border"};
  Configurable<bool> removeNoTimeFrameBorder{"removeNoTimeFrameBorder", false, "Remove TF border"};

  TF1* Particle_Tpc_nSigma_shift = 0;
  Configurable<bool> enable_pT_shift_tpc_nSigma{"enable_pT_shift_tpc_nSigma", false, "Enable Data TPC nSigma recentering by TF1"};
  Configurable<std::vector<float>> parShiftPt{"parShiftPt", {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f}, "Parameters for pT shift."};

  TF1* Pion_Tpc_nSigma_shift = 0;
  TF1* Proton_Tpc_nSigma_shift = 0;
  TF1* Deuteron_Tpc_nSigma_shift = 0;
  TF1* Triton_Tpc_nSigma_shift = 0;
  TF1* He3_Tpc_nSigma_shift = 0;
  TF1* He4_Tpc_nSigma_shift = 0;

  Configurable<bool> enable_pT_shift_pion_tpc_nSigma{"enable_pT_shift_pion_tpc_nSigma", false, "Enable MC Pi plus TPC nSigma recentering by TF1"};
  Configurable<bool> enable_pT_shift_proton_tpc_nSigma{"enable_pT_shift_proton_tpc_nSigma", false, "Enable MC Proton TPC nSigma recentering by TF1"};
  Configurable<bool> enable_pT_shift_deuteron_tpc_nSigma{"enable_pT_shift_deuteron_tpc_nSigma", false, "Enable MC Deuteron TPC nSigma recentering by TF1"};
  Configurable<bool> enable_pT_shift_triton_tpc_nSigma{"enable_pT_shift_triton_tpc_nSigma", false, "Enable MC Triton TPC nSigma recentering by TF1"};
  Configurable<bool> enable_pT_shift_He3_tpc_nSigma{"enable_pT_shift_He3_tpc_nSigma", false, "Enable MC Helium-3 TPC nSigma recentering by TF1"};
  Configurable<bool> enable_pT_shift_He4_tpc_nSigma{"enable_pT_shift_He4_tpc_nSigma", false, "Enable MC Helium-4 TPC nSigma recentering by TF1"};
  Configurable<std::vector<float>> parShiftPtPion{"parShiftPtPion", {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f}, "MC Parameters for Pi plus pT shift."};
  Configurable<std::vector<float>> parShiftPtProton{"parShiftPtProton", {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f}, "MC Parameters for Proton pT shift."};
  Configurable<std::vector<float>> parShiftPtDeuteron{"parShiftPtDeuteron", {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f}, "MC Parameters for Deuteron pT shift."};
  Configurable<std::vector<float>> parShiftPtTriton{"parShiftPtTriton", {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f}, "MC Parameters for Triton pT shift."};
  Configurable<std::vector<float>> parShiftPtHe3{"parShiftPtHe3", {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f}, "MC Parameters for Helium-3 pT shift."};
  Configurable<std::vector<float>> parShiftPtHe4{"parShiftPtHe4", {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f}, "MC Parameters for Alpha pT shift."};

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
  void fillHistograms_spectra(const CollisionType& event, const TracksType& tracks)
  {
    if (event_selection_sel8 && event.sel8()) {
      spectra_reg.fill(HIST("histRecVtxZData"), event.posZ());
    }
    if (!event_selection_sel8) {
      spectra_reg.fill(HIST("histRecVtxZData"), event.posZ());
    }
    if (!isEventSelected(event))
      return;

    for (auto track : tracks) { // start loop over all tracks

      histTrackcuts_data_spectra->AddBinContent(1);
      if (event_selection_sel8 && !event.sel8())
        continue;
      histTrackcuts_data_spectra->AddBinContent(2);

      if (track.sign() > 0) {
        spectra_reg.fill(HIST("histDcaVsPtData_particle"), track.pt(), track.dcaXY());
        spectra_reg.fill(HIST("histDcaZVsPtData_particle"), track.pt(), track.dcaZ());
      }
      if (track.sign() < 0) {
        spectra_reg.fill(HIST("histDcaVsPtData_antiparticle"), track.pt(), track.dcaXY());
        spectra_reg.fill(HIST("histDcaZVsPtData_antiparticle"), track.pt(), track.dcaZ());
      }
      spectra_reg.fill(HIST("histNClusterTPC"), track.pt(), track.tpcNClsCrossedRows());
      spectra_reg.fill(HIST("histNClusterITS"), track.pt(), track.itsNCls());
      spectra_reg.fill(HIST("histNClusterITSib"), track.pt(), track.itsNClsInnerBarrel());
      spectra_reg.fill(HIST("histChi2TPC"), track.pt(), track.tpcChi2NCl());
      spectra_reg.fill(HIST("histChi2ITS"), track.pt(), track.itsChi2NCl());

      float TPCnumberClsFound = track.tpcNClsFound();
      float TPC_nCls_Crossed_Rows = track.tpcNClsCrossedRows();
      float RatioCrossedRowsOverFindableTPC = track.tpcCrossedRowsOverFindableCls();
      float Chi2perClusterTPC = track.tpcChi2NCl();
      float Chi2perClusterITS = track.itsChi2NCl();

      bool insideDCAxy = (std::abs(track.dcaXY()) <= (maxDcaXYFactor.value * (0.0105f + 0.0350f / pow(track.pt(), 1.1f))));

      if (!(insideDCAxy) || TMath::Abs(track.dcaZ()) > maxDCA_Z)
        continue;
      histTrackcuts_data_spectra->AddBinContent(3);
      if (TPCnumberClsFound < minTPCnClsFound || TPC_nCls_Crossed_Rows < minNCrossedRowsTPC)
        continue;
      histTrackcuts_data_spectra->AddBinContent(4);
      if (RatioCrossedRowsOverFindableTPC < minRatioCrossedRowsTPC || RatioCrossedRowsOverFindableTPC > maxRatioCrossedRowsTPC)
        continue;
      histTrackcuts_data_spectra->AddBinContent(5);
      if (Chi2perClusterTPC > maxChi2PerClusterTPC || Chi2perClusterTPC < minChi2PerClusterTPC || Chi2perClusterITS > maxChi2PerClusterITS)
        continue;
      histTrackcuts_data_spectra->AddBinContent(6);
      if (!(track.passedTPCRefit()))
        continue;
      histTrackcuts_data_spectra->AddBinContent(7);
      if (!(track.passedITSRefit()))
        continue;
      histTrackcuts_data_spectra->AddBinContent(8);
      if ((track.itsNClsInnerBarrel()) < minReqClusterITSib || (track.itsNCls()) < minReqClusterITS)
        continue;
      histTrackcuts_data_spectra->AddBinContent(9);
      if (track.pt() < p_min || track.pt() > p_max)
        continue;
      histTrackcuts_data_spectra->AddBinContent(10);
      if ((requireITS && !(track.hasITS())) || (requireTPC && !(track.hasTPC())))
        continue;
      histTrackcuts_data_spectra->AddBinContent(11);
      if (requireGoldenChi2 && !(track.passedGoldenChi2()))
        continue;
      histTrackcuts_data_spectra->AddBinContent(12);

      spectra_reg.fill(HIST("histTpcSignalData"), track.pt() * track.sign(), track.tpcSignal());
      if (track.hasTOF()) {
        spectra_reg.fill(HIST("histTOFm2"), track.pt(), track.mass() * track.mass());
        spectra_reg.fill(HIST("histTofSignalData"), track.pt() * track.sign(), track.beta());
      }
    }
  }

  //***********************************************************************************

  template <typename CollisionType, typename TracksType>
  void fillHistograms_particle(const CollisionType& event, const TracksType& tracks, const int Partilce_type)
  {

    if (!isEventSelected(event))
      return;

    if (enable_pT_shift_tpc_nSigma) {
      Particle_Tpc_nSigma_shift = new TF1("Particle_Tpc_nSigma_shift", "[0] * TMath::Exp([1] + [2] * x) + [3] + [4] * x + [5] * x * x", 0.f, 14.f);
      auto par = (std::vector<float>)parShiftPt;
      Particle_Tpc_nSigma_shift->SetParameters(par[0], par[1], par[2], par[3], par[4], par[5]);
    }
    for (auto track : tracks) {

      float TPCnSigma_particle = -100;
      float TOFnSigma_particle = -100;
      float momentum;
      TLorentzVector lorentzVector_particle{};
      HistogramRegistry particle_reg;
      HistogramRegistry aparticle_reg;

      switch (Partilce_type) {
        case 0: // pi plus/minus
          lorentzVector_particle.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassPiPlus);
          TPCnSigma_particle = track.tpcNSigmaPi();
          TOFnSigma_particle = track.tofNSigmaPi();
          particle_reg = pion_reg;
          aparticle_reg = apion_reg;
          momentum = track.pt();
          break;
        case 1: // (anti)proton
          lorentzVector_particle.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassProton);
          TPCnSigma_particle = track.tpcNSigmaPr();
          TOFnSigma_particle = track.tofNSigmaPr();
          particle_reg = proton_reg;
          aparticle_reg = aproton_reg;
          momentum = track.pt();
          break;
        case 2: // (anti)deuteron
          lorentzVector_particle.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassDeuteron);
          TPCnSigma_particle = track.tpcNSigmaDe();
          TOFnSigma_particle = track.tofNSigmaDe();
          particle_reg = deuteron_reg;
          aparticle_reg = adeuteron_reg;
          momentum = track.pt();
          break;
        case 3: // (anti)triton
          lorentzVector_particle.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassTriton);
          TPCnSigma_particle = track.tpcNSigmaTr();
          TOFnSigma_particle = track.tofNSigmaTr();
          particle_reg = triton_reg;
          aparticle_reg = atriton_reg;
          momentum = track.pt();
          break;
        case 4: // (anti)Helium-3
          lorentzVector_particle.SetPtEtaPhiM(track.pt() * 2.0, track.eta(), track.phi(), constants::physics::MassHelium3);
          TPCnSigma_particle = track.tpcNSigmaHe();
          TOFnSigma_particle = track.tofNSigmaHe();
          particle_reg = Helium3_reg;
          aparticle_reg = aHelium3_reg;
          momentum = track.pt() * 2.0;
          break;
        case 5: // (anti)Helium-4
          lorentzVector_particle.SetPtEtaPhiM(track.pt() * 2.0, track.eta(), track.phi(), constants::physics::MassAlpha);
          TPCnSigma_particle = track.tpcNSigmaAl();
          TOFnSigma_particle = track.tofNSigmaAl();
          particle_reg = Helium4_reg;
          aparticle_reg = aHelium4_reg;
          momentum = track.pt() * 2.0;
          break;
        default:
          continue;
          break;
      }

      if (enable_pT_shift_tpc_nSigma) {
        float nSigma_shift = Particle_Tpc_nSigma_shift->Eval(momentum);
        TPCnSigma_particle -= nSigma_shift;
      }

      histTrackcuts_data_particle->AddBinContent(1);
      if (event_selection_sel8 && !event.sel8())
        continue;
      histTrackcuts_data_particle->AddBinContent(2);

      float TPCnumberClsFound = track.tpcNClsFound();
      float TPC_nCls_Crossed_Rows = track.tpcNClsCrossedRows();
      float RatioCrossedRowsOverFindableTPC = track.tpcCrossedRowsOverFindableCls();
      float Chi2perClusterTPC = track.tpcChi2NCl();
      float Chi2perClusterITS = track.itsChi2NCl();

      if (lorentzVector_particle.Rapidity() < yMin || lorentzVector_particle.Rapidity() > yMax)
        continue;
      histTrackcuts_data_particle->AddBinContent(3);
      if (TPCnumberClsFound < minTPCnClsFound || TPC_nCls_Crossed_Rows < minNCrossedRowsTPC)
        continue;
      histTrackcuts_data_particle->AddBinContent(4);
      if (RatioCrossedRowsOverFindableTPC < minRatioCrossedRowsTPC || RatioCrossedRowsOverFindableTPC > maxRatioCrossedRowsTPC)
        continue;
      histTrackcuts_data_particle->AddBinContent(5);
      if (Chi2perClusterTPC > maxChi2PerClusterTPC || Chi2perClusterTPC < minChi2PerClusterTPC || Chi2perClusterITS > maxChi2PerClusterITS)
        continue;
      histTrackcuts_data_particle->AddBinContent(6);
      if (!(track.passedTPCRefit()))
        continue;
      histTrackcuts_data_particle->AddBinContent(7);
      if (!(track.passedITSRefit()))
        continue;
      histTrackcuts_data_particle->AddBinContent(8);
      if ((track.itsNClsInnerBarrel()) < minReqClusterITSib || (track.itsNCls()) < minReqClusterITS)
        continue;
      histTrackcuts_data_particle->AddBinContent(9);
      if (momentum < p_min || momentum > p_max)
        continue;
      histTrackcuts_data_particle->AddBinContent(10);
      if ((requireITS && !(track.hasITS())) || (requireTPC && !(track.hasTPC())))
        continue;
      histTrackcuts_data_particle->AddBinContent(11);
      if (requireGoldenChi2 && !(track.passedGoldenChi2()))
        continue;

      histTrackcuts_data_particle->AddBinContent(12);

      if (TPCnSigma_particle > nsigmacutLow && TPCnSigma_particle < nsigmacutHigh && TMath::Abs(track.dcaZ()) < 2.0 && TMath::Abs(track.dcaXY()) < 0.5) {
        if (track.sign() > 0) {
          particle_reg.fill(HIST("histDcaVsPtData"), momentum, track.dcaXY());
          particle_reg.fill(HIST("histDcaZVsPtData"), momentum, track.dcaZ());
        }
        if (track.sign() < 0) {
          aparticle_reg.fill(HIST("histDcaVsPtData"), momentum, track.dcaXY());
          aparticle_reg.fill(HIST("histDcaZVsPtData"), momentum, track.dcaZ());
        }
      }

      bool insideDCAxy = (std::abs(track.dcaXY()) <= (maxDcaXYFactor.value * (0.0105f + 0.0350f / pow(momentum, 1.1f))));

      if (!(insideDCAxy) || TMath::Abs(track.dcaZ()) > maxDCA_Z)
        continue;
      histTrackcuts_data_particle->AddBinContent(13);

      //******* Fill particle histograms ***********
      if (track.sign() > 0) {
        particle_reg.fill(HIST("histTpcNsigmaData"), momentum, TPCnSigma_particle);
      }

      if (track.sign() < 0) {
        aparticle_reg.fill(HIST("histTpcNsigmaData"), momentum, TPCnSigma_particle);
      }

      if (TPCnSigma_particle > nsigmacutLow && TPCnSigma_particle < nsigmacutHigh) {

        if (track.sign() > 0) {
          particle_reg.fill(HIST("histDcaVsPtData_after_cut"), momentum, track.dcaXY());
          particle_reg.fill(HIST("histDcaZVsPtData_after_cut"), momentum, track.dcaZ());
          particle_reg.fill(HIST("histTpcSignalData"), momentum, track.tpcSignal());
          particle_reg.fill(HIST("histNClusterTPC"), momentum, track.tpcNClsFound());
          particle_reg.fill(HIST("histNClusterITS"), momentum, track.itsNCls());
          particle_reg.fill(HIST("histNClusterITSib"), momentum, track.itsNClsInnerBarrel());
          particle_reg.fill(HIST("histChi2TPC"), momentum, track.tpcChi2NCl());
          particle_reg.fill(HIST("histChi2ITS"), momentum, track.itsChi2NCl());

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
            particle_reg.fill(HIST("histTOFm2"), momentum, TOFmass2);
            particle_reg.fill(HIST("histTofSignalData"), momentum, beta);
            particle_reg.fill(HIST("histTofNsigmaData"), momentum, TOFnSigma_particle);
          }
        }

        if (track.sign() < 0) {
          aparticle_reg.fill(HIST("histDcaVsPtData_after_cut"), momentum, track.dcaXY());
          aparticle_reg.fill(HIST("histDcaZVsPtData_after_cut"), momentum, track.dcaZ());
          aparticle_reg.fill(HIST("histTpcSignalData"), momentum, track.tpcSignal());
          aparticle_reg.fill(HIST("histNClusterTPC"), momentum, track.tpcNClsFound());
          aparticle_reg.fill(HIST("histNClusterITS"), momentum, track.itsNCls());
          aparticle_reg.fill(HIST("histNClusterITSib"), momentum, track.itsNClsInnerBarrel());
          aparticle_reg.fill(HIST("histChi2TPC"), momentum, track.tpcChi2NCl());
          aparticle_reg.fill(HIST("histChi2ITS"), momentum, track.itsChi2NCl());

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
            aparticle_reg.fill(HIST("histTOFm2"), momentum, TOFmass2);
            aparticle_reg.fill(HIST("histTofSignalData"), momentum, beta);
            aparticle_reg.fill(HIST("histTofNsigmaData"), momentum, TOFnSigma_particle);
          }
        }
      }
    }
  }

  //***********************************************************************************

  template <typename CollisionType, typename TracksType>
  void fillHistograms_particle_cent(const CollisionType& event, const TracksType& tracks, const int Partilce_type)
  {
    if (!isEventSelected(event))
      return;

    if (enable_pT_shift_tpc_nSigma) {
      Particle_Tpc_nSigma_shift = new TF1("Particle_Tpc_nSigma_shift", "[0] * TMath::Exp([1] + [2] * x) + [3] + [4] * x + [5] * x * x", 0.f, 14.f);
      auto par = (std::vector<float>)parShiftPt;
      Particle_Tpc_nSigma_shift->SetParameters(par[0], par[1], par[2], par[3], par[4], par[5]);
    }

    if (event_selection_sel8 && event.sel8())
      spectra_reg.fill(HIST("histCentrality"), event.centFT0C());
    if (!event_selection_sel8)
      spectra_reg.fill(HIST("histCentrality"), event.centFT0C());

    for (auto track : tracks) {
      if ((event_selection_sel8 && !event.sel8()) || (enable_Centrality_cut_global && (event.centFT0C() < minCentrality) && (event.centFT0C() > maxCentrality)))
        continue;

      float TPCnSigma_particle = -100;
      float TOFnSigma_particle = -100;
      float momentum;
      TLorentzVector lorentzVector_particle{};
      HistogramRegistry particle_reg;
      HistogramRegistry aparticle_reg;

      switch (Partilce_type) {
        case 0: // pi plus/minus
          lorentzVector_particle.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassPiPlus);
          TPCnSigma_particle = track.tpcNSigmaPi();
          TOFnSigma_particle = track.tofNSigmaPi();
          particle_reg = pion_reg;
          aparticle_reg = apion_reg;
          momentum = track.pt();
          break;
        case 1: // (anti)proton
          lorentzVector_particle.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassProton);
          TPCnSigma_particle = track.tpcNSigmaPr();
          TOFnSigma_particle = track.tofNSigmaPr();
          particle_reg = proton_reg;
          aparticle_reg = aproton_reg;
          momentum = track.pt();
          break;
        case 2: // (anti)deuteron
          lorentzVector_particle.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassDeuteron);
          TPCnSigma_particle = track.tpcNSigmaDe();
          TOFnSigma_particle = track.tofNSigmaDe();
          particle_reg = deuteron_reg;
          aparticle_reg = adeuteron_reg;
          momentum = track.pt();
          break;
        case 3: // (anti)triton
          lorentzVector_particle.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassTriton);
          TPCnSigma_particle = track.tpcNSigmaTr();
          TOFnSigma_particle = track.tofNSigmaTr();
          particle_reg = triton_reg;
          aparticle_reg = atriton_reg;
          momentum = track.pt();
          break;
        case 4: // (anti)Helium-3
          lorentzVector_particle.SetPtEtaPhiM(track.pt() * 2.0, track.eta(), track.phi(), constants::physics::MassHelium3);
          TPCnSigma_particle = track.tpcNSigmaHe();
          TOFnSigma_particle = track.tofNSigmaHe();
          particle_reg = Helium3_reg;
          aparticle_reg = aHelium3_reg;
          momentum = track.pt() * 2.0;
          break;
        case 5: // (anti)Helium-4
          lorentzVector_particle.SetPtEtaPhiM(track.pt() * 2.0, track.eta(), track.phi(), constants::physics::MassAlpha);
          TPCnSigma_particle = track.tpcNSigmaAl();
          TOFnSigma_particle = track.tofNSigmaAl();
          particle_reg = Helium4_reg;
          aparticle_reg = aHelium4_reg;
          momentum = track.pt() * 2.0;
          break;
        default:
          continue;
          break;
      }

      if (enable_pT_shift_tpc_nSigma) {
        float nSigma_shift = Particle_Tpc_nSigma_shift->Eval(momentum);
        TPCnSigma_particle -= nSigma_shift;
      }

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

      if (lorentzVector_particle.Rapidity() < yMin || lorentzVector_particle.Rapidity() > yMax)
        continue;

      bool insideDCAxy = (std::abs(track.dcaXY()) <= (maxDcaXYFactor.value * (0.0105f + 0.0350f / pow(momentum, 1.1f))));

      if (!(insideDCAxy) || TMath::Abs(track.dcaZ()) > maxDCA_Z)
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
      if (momentum < p_min || momentum > p_max)
        continue;
      if ((requireITS && !(track.hasITS())) || (requireTPC && !(track.hasTPC())))
        continue;
      if (requireGoldenChi2 && !(track.passedGoldenChi2()))
        continue;

      if (track.sign() > 0) {
        particle_reg.fill(HIST("histTpcNsigmaData_cent"), momentum, TPCnSigma_particle, event.centFT0C());
        if ((event.centFT0C() > minCentrality) && (event.centFT0C() < maxCentrality)) {
          particle_reg.fill(HIST("histTpcNsigmaData_eta"), momentum, TPCnSigma_particle, track.eta());
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
          particle_reg.fill(HIST("histTofNsigmaData_cent"), momentum, TOFnSigma_particle, event.centFT0C());
          particle_reg.fill(HIST("histTofm2_cent"), momentum, track.mass() * track.mass(), event.centFT0C());

          if ((event.centFT0C() > minCentrality) && (event.centFT0C() < maxCentrality)) {
            particle_reg.fill(HIST("histTofm2_eta"), momentum, track.mass() * track.mass(), track.eta());
            particle_reg.fill(HIST("histTofNsigmaData_eta"), momentum, TOFnSigma_particle, track.eta());
          }
        }
      }
      if (track.sign() < 0) {
        aparticle_reg.fill(HIST("histTpcNsigmaData_cent"), momentum, TPCnSigma_particle, event.centFT0C());
        if ((event.centFT0C() > minCentrality) && (event.centFT0C() < maxCentrality)) {
          aparticle_reg.fill(HIST("histTpcNsigmaData_eta"), momentum, TPCnSigma_particle, track.eta());
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
          aparticle_reg.fill(HIST("histTofNsigmaData_cent"), momentum, TOFnSigma_particle, event.centFT0C());
          aparticle_reg.fill(HIST("histTofm2_cent"), momentum, track.mass() * track.mass(), event.centFT0C());

          if ((event.centFT0C() > minCentrality) && (event.centFT0C() < maxCentrality)) {
            aparticle_reg.fill(HIST("histTofm2_eta"), momentum, track.mass() * track.mass(), track.eta());
            aparticle_reg.fill(HIST("histTofNsigmaData_eta"), momentum, TOFnSigma_particle, track.eta());
          }
        }
      }
    }
  }

  //***********************************************************************************

  Filter collisionFilter = (nabs(aod::collision::posZ) < cfgCutVertex);
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta && requireGlobalTrackWoDCAInFilter());

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>;

  using EventCandidatesCent = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>>;

  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCLfFullPi, aod::pidTOFFullPi, aod::pidTPCLfFullKa, aod::pidTOFFullKa, aod::pidTPCLfFullPr, aod::pidTOFFullPr, aod::pidTPCLfFullDe, aod::pidTOFFullDe, aod::pidTPCLfFullTr, aod::pidTOFFullTr, aod::pidTPCLfFullHe, aod::pidTOFFullHe, aod::pidTPCLfFullAl, aod::pidTOFFullAl, aod::TrackSelection, aod::TrackSelectionExtension, aod::TOFSignal, aod::pidTOFmass, aod::pidTOFbeta>>;

  void processData(EventCandidates::iterator const& event, TrackCandidates const& tracks)
  {
    fillHistograms_spectra(event, tracks);
    if (do_pion)
      fillHistograms_particle(event, tracks, 0); // pion
    if (do_proton)
      fillHistograms_particle(event, tracks, 1); // proton
    if (do_deuteron)
      fillHistograms_particle(event, tracks, 2); // deuteron
    if (do_triton)
      fillHistograms_particle(event, tracks, 3); // triton
    if (do_He3)
      fillHistograms_particle(event, tracks, 4); // He3
    if (do_He4)
      fillHistograms_particle(event, tracks, 5); // He4
  }
  PROCESS_SWITCH(NucleiHistTask, processData, "process data", true);

  void processDataCent(EventCandidatesCent::iterator const& event, TrackCandidates const& tracks)
  {
    fillHistograms_spectra(event, tracks);
    if (do_pion)
      fillHistograms_particle(event, tracks, 0); // pion
    if (do_proton)
      fillHistograms_particle(event, tracks, 1); // proton
    if (do_deuteron)
      fillHistograms_particle(event, tracks, 2); // deuteron
    if (do_triton)
      fillHistograms_particle(event, tracks, 3); // triton
    if (do_He3)
      fillHistograms_particle(event, tracks, 4); // He3
    if (do_He4)
      fillHistograms_particle(event, tracks, 5); // He4
    if (do_pion)
      fillHistograms_particle_cent(event, tracks, 0); // pion
    if (do_proton)
      fillHistograms_particle_cent(event, tracks, 1); // proton
    if (do_deuteron)
      fillHistograms_particle_cent(event, tracks, 2); // deuteron
    if (do_triton)
      fillHistograms_particle_cent(event, tracks, 3); // triton
    if (do_He3)
      fillHistograms_particle_cent(event, tracks, 4); // He3
    if (do_He4)
      fillHistograms_particle_cent(event, tracks, 5); // He4
  }
  PROCESS_SWITCH(NucleiHistTask, processDataCent, "process data with centralities", false);

  //********************** MC ****************************************

  void processMCgen(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles)
  {
    MC_gen_reg.fill(HIST("histGenVtxMC"), mcCollision.posZ());
    MC_gen_reg.fill(HIST("histCentrality"), mcCollision.impactParameter());

    for (const auto& mcParticleGen : mcParticles) {
      if (require_PhysicalPrimary_MC_gen && !mcParticleGen.isPhysicalPrimary())
        continue;
      int pdgCode = mcParticleGen.pdgCode();

      if (mcParticleGen.y() > yMax || mcParticleGen.y() < yMin)
        continue;
      if (mcParticleGen.eta() > cfgCutEta || mcParticleGen.eta() < -cfgCutEta)
        continue;

      int pdgbin = 0;
      switch (pdgCode) {
        case +211:
          histPDG_gen->AddBinContent(1);
          pdgbin = 0;
          break;
        case -211:
          histPDG_gen->AddBinContent(2);
          pdgbin = 1;
          break;
        case +321:
          histPDG_gen->AddBinContent(3);
          pdgbin = 2;
          break;
        case -321:
          histPDG_gen->AddBinContent(4);
          pdgbin = 3;
          break;
        case +2212:
          histPDG_gen->AddBinContent(5);
          pdgbin = 4;
          break;
        case -2212:
          histPDG_gen->AddBinContent(6);
          pdgbin = 5;
          break;
        case +1000010020:
          histPDG_gen->AddBinContent(7);
          pdgbin = 6;
          break;
        case -1000010020:
          histPDG_gen->AddBinContent(8);
          pdgbin = 7;
          break;
        case +1000010030:
          histPDG_gen->AddBinContent(9);
          pdgbin = 8;
          break;
        case -1000010030:
          histPDG_gen->AddBinContent(10);
          pdgbin = 9;
          break;
        case +1000020030:
          histPDG_gen->AddBinContent(11);
          pdgbin = 10;
          break;
        case -1000020030:
          histPDG_gen->AddBinContent(12);
          pdgbin = 11;
          break;
        case +1000020040:
          histPDG_gen->AddBinContent(13);
          pdgbin = 12;
          break;
        case -1000020040:
          histPDG_gen->AddBinContent(14);
          pdgbin = 13;
          break;
        default:
          break;
      }
      MC_gen_reg.fill(HIST("histEta"), mcParticleGen.eta(), pdgbin);
      if ((pdgCode == 1000020030) || (pdgCode == -1000020030) || (pdgCode == 1000020040) || (pdgCode == -1000020040)) {
        MC_gen_reg.fill(HIST("histPt"), mcParticleGen.pt() * 2.0, pdgbin);
      } else {
        MC_gen_reg.fill(HIST("histPt"), mcParticleGen.pt(), pdgbin);
      }
    }
  }
  PROCESS_SWITCH(NucleiHistTask, processMCgen, "process generated MC", false);

  void processMCreco(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0Cs>::iterator const& collisions, soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCLfFullPi, aod::pidTOFFullPi, aod::pidTPCLfFullKa, aod::pidTOFFullKa, aod::pidTPCLfFullPr, aod::pidTOFFullPr, aod::pidTPCLfFullDe, aod::pidTOFFullDe, aod::pidTPCLfFullTr, aod::pidTOFFullTr, aod::pidTPCLfFullHe, aod::pidTOFFullHe, aod::pidTPCLfFullAl, aod::pidTOFFullAl, aod::McTrackLabels, aod::TrackSelection, aod::TrackSelectionExtension, aod::TOFSignal, aod::pidTOFmass, aod::pidTOFbeta>> const& tracks,
                     aod::McParticles& /*mcParticles*/, aod::McCollisions const& /*mcCollisions*/)
  {

    if (event_selection_MC_sel8 && !collisions.sel8())
      return;
    MC_recon_reg.fill(HIST("histRecVtxMC"), collisions.posZ());
    MC_recon_reg.fill(HIST("histCentrality"), collisions.centFT0C());
    if (!isEventSelected(collisions))
      return;

    if (enable_pT_shift_pion_tpc_nSigma) {
      Pion_Tpc_nSigma_shift = new TF1("Pion_Tpc_nSigma_shift", "[0] * TMath::Exp([1] + [2] * x) + [3] + [4] * x + [5] * x * x", 0.f, 14.f);
      auto par = (std::vector<float>)parShiftPtPion;
      Pion_Tpc_nSigma_shift->SetParameters(par[0], par[1], par[2], par[3], par[4], par[5]);
    }
    if (enable_pT_shift_proton_tpc_nSigma) {
      Proton_Tpc_nSigma_shift = new TF1("Proton_Tpc_nSigma_shift", "[0] * TMath::Exp([1] + [2] * x) + [3] + [4] * x + [5] * x * x", 0.f, 14.f);
      auto par = (std::vector<float>)parShiftPtProton;
      Proton_Tpc_nSigma_shift->SetParameters(par[0], par[1], par[2], par[3], par[4], par[5]);
    }
    if (enable_pT_shift_deuteron_tpc_nSigma) {
      Deuteron_Tpc_nSigma_shift = new TF1("Deuteron_Tpc_nSigma_shift", "[0] * TMath::Exp([1] + [2] * x) + [3] + [4] * x + [5] * x * x", 0.f, 14.f);
      auto par = (std::vector<float>)parShiftPtDeuteron;
      Deuteron_Tpc_nSigma_shift->SetParameters(par[0], par[1], par[2], par[3], par[4], par[5]);
    }
    if (enable_pT_shift_triton_tpc_nSigma) {
      Triton_Tpc_nSigma_shift = new TF1("Triton_Tpc_nSigma_shift", "[0] * TMath::Exp([1] + [2] * x) + [3] + [4] * x + [5] * x * x", 0.f, 14.f);
      auto par = (std::vector<float>)parShiftPtTriton;
      Triton_Tpc_nSigma_shift->SetParameters(par[0], par[1], par[2], par[3], par[4], par[5]);
    }
    if (enable_pT_shift_He3_tpc_nSigma) {
      He3_Tpc_nSigma_shift = new TF1("He3_Tpc_nSigma_shift", "[0] * TMath::Exp([1] + [2] * x) + [3] + [4] * x + [5] * x * x", 0.f, 14.f);
      auto par = (std::vector<float>)parShiftPtHe3;
      He3_Tpc_nSigma_shift->SetParameters(par[0], par[1], par[2], par[3], par[4], par[5]);
    }
    if (enable_pT_shift_He4_tpc_nSigma) {
      He4_Tpc_nSigma_shift = new TF1("He4_Tpc_nSigma_shift", "[0] * TMath::Exp([1] + [2] * x) + [3] + [4] * x + [5] * x * x", 0.f, 14.f);
      auto par = (std::vector<float>)parShiftPtHe4;
      He4_Tpc_nSigma_shift->SetParameters(par[0], par[1], par[2], par[3], par[4], par[5]);
    }

    for (auto& track : tracks) {
      histTrackcuts_MC->AddBinContent(1);
      const auto particle = track.mcParticle();

      int pdgbin = 0;
      TLorentzVector lorentzVector_particle_MC{};
      switch (particle.pdgCode()) {
        case +211:
          pdgbin = 0;
          histPDG_reco->AddBinContent(1);
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassPiPlus);
          break;
        case -211:
          pdgbin = 1;
          histPDG_reco->AddBinContent(2);
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassPiPlus);
          break;
        case +321:
          pdgbin = 2;
          histPDG_reco->AddBinContent(3);
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassKPlus);
          break;
        case -321:
          pdgbin = 3;
          histPDG_reco->AddBinContent(4);
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassKPlus);
          break;
        case +2212:
          pdgbin = 4;
          histPDG_reco->AddBinContent(5);
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassProton);
          break;
        case -2212:
          pdgbin = 5;
          histPDG_reco->AddBinContent(6);
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassProton);
          break;
        case +1000010020:
          pdgbin = 6;
          histPDG_reco->AddBinContent(7);
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassDeuteron);
          break;
        case -1000010020:
          pdgbin = 7;
          histPDG_reco->AddBinContent(8);
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassDeuteron);
          break;
        case +1000010030:
          pdgbin = 8;
          histPDG_reco->AddBinContent(9);
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassTriton);
          break;
        case -1000010030:
          pdgbin = 9;
          histPDG_reco->AddBinContent(10);
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassTriton);
          break;
        case +1000020030:
          pdgbin = 10;
          histPDG_reco->AddBinContent(11);
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt() * 2.0, track.eta(), track.phi(), constants::physics::MassHelium3);
          break;
        case -1000020030:
          pdgbin = 11;
          histPDG_reco->AddBinContent(12);
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt() * 2.0, track.eta(), track.phi(), constants::physics::MassHelium3);
          break;
        case +1000020040:
          pdgbin = 12;
          histPDG_reco->AddBinContent(13);
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt() * 2.0, track.eta(), track.phi(), constants::physics::MassAlpha);
          break;
        case -1000020040:
          pdgbin = 13;
          histPDG_reco->AddBinContent(14);
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt() * 2.0, track.eta(), track.phi(), constants::physics::MassAlpha);
          break;
        default:
          pdgbin = -1;
          break;
      }

      if (require_PhysicalPrimary_MC_reco && !particle.isPhysicalPrimary())
        continue;
      histTrackcuts_MC->AddBinContent(2);

      if (lorentzVector_particle_MC.Rapidity() < yMin || lorentzVector_particle_MC.Rapidity() > yMax)
        continue;
      histTrackcuts_MC->AddBinContent(3);

      MC_recon_reg.fill(HIST("histEta"), track.eta(), pdgbin);

      if ((particle.pdgCode() == 1000020030) || (particle.pdgCode() == -1000020030) || (particle.pdgCode() == 1000020040) || (particle.pdgCode() == -1000020040)) {
        MC_recon_reg.fill(HIST("histPt"), track.pt() * 2.0, pdgbin);
        MC_recon_reg.fill(HIST("histDCA"), track.pt() * 2.0, track.dcaXY(), pdgbin);
        MC_recon_reg.fill(HIST("histDCAz"), track.pt() * 2.0, track.dcaZ(), pdgbin);
        MC_recon_reg.fill(HIST("histTpcSignalData"), track.pt() * 2.0 * track.sign(), track.tpcSignal(), pdgbin);
        MC_recon_reg.fill(HIST("histTpcSignalData_all_species"), track.pt() * 2.0 * track.sign(), track.tpcSignal());
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
        MC_recon_reg.fill(HIST("histTpcSignalData"), track.pt() * track.sign(), track.tpcSignal(), pdgbin);
        MC_recon_reg.fill(HIST("histTpcSignalData_all_species"), track.pt() * track.sign(), track.tpcSignal());
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

      bool insideDCAxy = (std::abs(track.dcaXY()) <= (maxDcaXYFactor.value * (0.0105f + 0.0350f / pow(track.pt(), 1.1f))));

      if (!(insideDCAxy) || TMath::Abs(track.dcaZ()) > maxDCA_Z)
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

      if (enable_pT_shift_pion_tpc_nSigma) {
        float nSigmaPion_shift = Pion_Tpc_nSigma_shift->Eval(track.pt());
        nSigmaPion -= nSigmaPion_shift;
      }
      if (enable_pT_shift_proton_tpc_nSigma) {
        float nSigmaProton_shift = Proton_Tpc_nSigma_shift->Eval(track.pt());
        nSigmaProton -= nSigmaProton_shift;
      }
      if (enable_pT_shift_deuteron_tpc_nSigma) {
        float nSigmaDeuteron_shift = Deuteron_Tpc_nSigma_shift->Eval(track.pt());
        nSigmaDeuteron -= nSigmaDeuteron_shift;
      }
      if (enable_pT_shift_triton_tpc_nSigma) {
        float nSigmaTriton_shift = Triton_Tpc_nSigma_shift->Eval(track.pt());
        nSigmaTriton -= nSigmaTriton_shift;
      }
      if (enable_pT_shift_He3_tpc_nSigma) {
        float nSigmaHe3_shift = He3_Tpc_nSigma_shift->Eval(track.pt() * 2.0);
        nSigmaHe3 -= nSigmaHe3_shift;
      }
      if (enable_pT_shift_He4_tpc_nSigma) {
        float nSigmaHe4_shift = He4_Tpc_nSigma_shift->Eval(track.pt() * 2.0);
        nSigmaHe4 -= nSigmaHe4_shift;
      }

      if (track.sign() > 0) {
        MC_recon_reg.fill(HIST("histTpcNsigmaDataPi"), track.pt(), nSigmaPion);
        MC_recon_reg.fill(HIST("histTpcNsigmaDataPr"), track.pt(), nSigmaProton);
        MC_recon_reg.fill(HIST("histTpcNsigmaDataDe"), track.pt(), nSigmaDeuteron);
        MC_recon_reg.fill(HIST("histTpcNsigmaDataTr"), track.pt(), nSigmaTriton);
        MC_recon_reg.fill(HIST("histTpcNsigmaDataHe"), track.pt() * 2.0, nSigmaHe3);
        MC_recon_reg.fill(HIST("histTpcNsigmaDataAl"), track.pt() * 2.0, nSigmaHe4);
      }
      if (track.sign() < 0) {
        MC_recon_reg.fill(HIST("histTpcNsigmaDataaPi"), track.pt(), nSigmaPion);
        MC_recon_reg.fill(HIST("histTpcNsigmaDataaPr"), track.pt(), nSigmaProton);
        MC_recon_reg.fill(HIST("histTpcNsigmaDataaDe"), track.pt(), nSigmaDeuteron);
        MC_recon_reg.fill(HIST("histTpcNsigmaDataaTr"), track.pt(), nSigmaTriton);
        MC_recon_reg.fill(HIST("histTpcNsigmaDataaHe"), track.pt() * 2.0, nSigmaHe3);
        MC_recon_reg.fill(HIST("histTpcNsigmaDataaAl"), track.pt() * 2.0, nSigmaHe4);
      }

      if (track.hasTOF()) {
        Float_t TOFmass2 = ((track.mass()) * (track.mass()));

        MC_recon_reg.fill(HIST("histTOFm2"), track.pt(), TOFmass2, pdgbin);
        MC_recon_reg.fill(HIST("histTofSignalData"), track.pt() * track.sign(), track.beta(), pdgbin);
        MC_recon_reg.fill(HIST("histTofSignalData_all_species"), track.pt() * track.sign(), track.beta());

        if (track.sign() > 0) {
          if (nSigmaPion > nsigmacutLow && nSigmaPion < nsigmacutHigh)
            MC_recon_reg.fill(HIST("histTofNsigmaDataPi"), track.pt(), track.tofNSigmaPi());
          if (nSigmaProton > nsigmacutLow && nSigmaProton < nsigmacutHigh)
            MC_recon_reg.fill(HIST("histTofNsigmaDataPr"), track.pt(), track.tofNSigmaPr());
          if (nSigmaDeuteron > nsigmacutLow && nSigmaDeuteron < nsigmacutHigh)
            MC_recon_reg.fill(HIST("histTofNsigmaDataDe"), track.pt(), track.tofNSigmaDe());
          if (nSigmaTriton > nsigmacutLow && nSigmaTriton < nsigmacutHigh)
            MC_recon_reg.fill(HIST("histTofNsigmaDataTr"), track.pt(), track.tofNSigmaTr());
          if (nSigmaHe3 > nsigmacutLow && nSigmaHe3 < nsigmacutHigh)
            MC_recon_reg.fill(HIST("histTofNsigmaDataHe"), track.pt() * 2.0, track.tofNSigmaHe());
          if (nSigmaHe4 > nsigmacutLow && nSigmaHe4 < nsigmacutHigh)
            MC_recon_reg.fill(HIST("histTofNsigmaDataAl"), track.pt() * 2.0, track.tofNSigmaAl());
        }
        if (track.sign() < 0) {
          if (nSigmaPion > nsigmacutLow && nSigmaPion < nsigmacutHigh)
            MC_recon_reg.fill(HIST("histTofNsigmaDataaPi"), track.pt(), track.tofNSigmaPi());
          if (nSigmaProton > nsigmacutLow && nSigmaProton < nsigmacutHigh)
            MC_recon_reg.fill(HIST("histTofNsigmaDataaPr"), track.pt(), track.tofNSigmaPr());
          if (nSigmaDeuteron > nsigmacutLow && nSigmaDeuteron < nsigmacutHigh)
            MC_recon_reg.fill(HIST("histTofNsigmaDataaDe"), track.pt(), track.tofNSigmaDe());
          if (nSigmaTriton > nsigmacutLow && nSigmaTriton < nsigmacutHigh)
            MC_recon_reg.fill(HIST("histTofNsigmaDataaTr"), track.pt(), track.tofNSigmaTr());
          if (nSigmaHe3 > nsigmacutLow && nSigmaHe3 < nsigmacutHigh)
            MC_recon_reg.fill(HIST("histTofNsigmaDataaHe"), track.pt() * 2.0, track.tofNSigmaHe());
          if (nSigmaHe4 > nsigmacutLow && nSigmaHe4 < nsigmacutHigh)
            MC_recon_reg.fill(HIST("histTofNsigmaDataaAl"), track.pt() * 2.0, track.tofNSigmaAl());
        }
      }
    }
  }
  PROCESS_SWITCH(NucleiHistTask, processMCreco, "process reconstructed MC", false);

  void processMCdca(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0Cs>::iterator const& collisions, soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCLfFullPi, aod::pidTOFFullPi, aod::pidTPCLfFullKa, aod::pidTOFFullKa, aod::pidTPCLfFullPr, aod::pidTOFFullPr, aod::pidTPCLfFullDe, aod::pidTOFFullDe, aod::pidTPCLfFullTr, aod::pidTOFFullTr, aod::pidTPCLfFullHe, aod::pidTOFFullHe, aod::pidTPCLfFullAl, aod::pidTOFFullAl, aod::McTrackLabels, aod::TrackSelection, aod::TrackSelectionExtension, aod::TOFSignal, aod::pidTOFmass, aod::pidTOFbeta>> const& tracks,
                    aod::McParticles& /*mcParticles*/, aod::McCollisions const& /*mcCollisions*/)
  {

    if (event_selection_MC_sel8 && !collisions.sel8())
      return;
    if (!isEventSelected(collisions))
      return;

    for (auto& track : tracks) {
      histTrackcuts_MC->AddBinContent(1);
      const auto particle = track.mcParticle();

      int pdgbin = 0;
      TLorentzVector lorentzVector_particle_MC{};
      switch (particle.pdgCode()) {
        case +211:
          pdgbin = 0;
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassPiPlus);
          break;
        case -211:
          pdgbin = 1;
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassPiPlus);
          break;
        case +321:
          pdgbin = 2;
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassKPlus);
          break;
        case -321:
          pdgbin = 3;
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassKPlus);
          break;
        case +2212:
          pdgbin = 4;
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassProton);
          break;
        case -2212:
          pdgbin = 5;
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassProton);
          break;
        case +1000010020:
          pdgbin = 6;
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassDeuteron);
          break;
        case -1000010020:
          pdgbin = 7;
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassDeuteron);
          break;
        case +1000010030:
          pdgbin = 8;
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassTriton);
          break;
        case -1000010030:
          pdgbin = 9;
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassTriton);
          break;
        case +1000020030:
          pdgbin = 10;
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt() * 2.0, track.eta(), track.phi(), constants::physics::MassHelium3);
          break;
        case -1000020030:
          pdgbin = 11;
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt() * 2.0, track.eta(), track.phi(), constants::physics::MassHelium3);
          break;
        case +1000020040:
          pdgbin = 12;
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt() * 2.0, track.eta(), track.phi(), constants::physics::MassAlpha);
          break;
        case -1000020040:
          pdgbin = 13;
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt() * 2.0, track.eta(), track.phi(), constants::physics::MassAlpha);
          break;
        default:
          pdgbin = -1;
          break;
      }

      if (lorentzVector_particle_MC.Rapidity() < yMin || lorentzVector_particle_MC.Rapidity() > yMax)
        continue;
      float TPCnumberClsFound = track.tpcNClsFound();
      float TPC_nCls_Crossed_Rows = track.tpcNClsCrossedRows();
      float RatioCrossedRowsOverFindableTPC = track.tpcCrossedRowsOverFindableCls();
      float Chi2perClusterTPC = track.tpcChi2NCl();
      float Chi2perClusterITS = track.itsChi2NCl();

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

      MC_DCA.fill(HIST("histEta"), track.eta(), pdgbin);

      if (particle.isPhysicalPrimary()) {
        if ((particle.pdgCode() == 1000020030) || (particle.pdgCode() == -1000020030) || (particle.pdgCode() == 1000020040) || (particle.pdgCode() == -1000020040)) {
          MC_DCA.fill(HIST("histDCA_prim"), track.pt() * 2.0, track.dcaXY(), pdgbin);
          MC_DCA.fill(HIST("histDCAz_prim"), track.pt() * 2.0, track.dcaZ(), pdgbin);
        } else {
          MC_DCA.fill(HIST("histDCA_prim"), track.pt(), track.dcaXY(), pdgbin);
          MC_DCA.fill(HIST("histDCAz_prim"), track.pt(), track.dcaZ(), pdgbin);
        }
      } else if (particle.getProcess() == 4) {
        if ((particle.pdgCode() == 1000020030) || (particle.pdgCode() == -1000020030) || (particle.pdgCode() == 1000020040) || (particle.pdgCode() == -1000020040)) {
          MC_DCA.fill(HIST("histDCA_weak"), track.pt() * 2.0, track.dcaXY(), pdgbin);
          MC_DCA.fill(HIST("histDCAz_weak"), track.pt() * 2.0, track.dcaZ(), pdgbin);
        } else {
          MC_DCA.fill(HIST("histDCA_weak"), track.pt(), track.dcaXY(), pdgbin);
          MC_DCA.fill(HIST("histDCAz_weak"), track.pt(), track.dcaZ(), pdgbin);
        }
      } else if (particle.getProcess() == 23) {
        if ((particle.pdgCode() == 1000020030) || (particle.pdgCode() == -1000020030) || (particle.pdgCode() == 1000020040) || (particle.pdgCode() == -1000020040)) {
          MC_DCA.fill(HIST("histDCA_mat"), track.pt() * 2.0, track.dcaXY(), pdgbin);
          MC_DCA.fill(HIST("histDCAz_mat"), track.pt() * 2.0, track.dcaZ(), pdgbin);
        } else {
          MC_DCA.fill(HIST("histDCA_mat"), track.pt(), track.dcaXY(), pdgbin);
          MC_DCA.fill(HIST("histDCAz_mat"), track.pt(), track.dcaZ(), pdgbin);
        }
      }
    }
  }
  PROCESS_SWITCH(NucleiHistTask, processMCdca, "process MC DCA", false);
};

//***********************************************************************************

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<NucleiHistTask>(cfgc, TaskName{"nuclei-hist"})};
}
