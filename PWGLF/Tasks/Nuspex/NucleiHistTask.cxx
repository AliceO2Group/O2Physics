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
  HistogramRegistry MC_truth_reg{"MC_particles", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry MC_recon_reg{"MC_particles_reconstructed", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry MC_recon_diff_reg{"MC_reconstructed_diff", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  OutputObj<TH1F> histPDG{TH1F("PDG", "PDG;PDG code", 100, 0.0, 100.0)};

  void init(o2::framework::InitContext&)
  {

    if (doprocessData == true && doprocessDataCent == true) {
      LOG(fatal) << "Can't enable processData and processDataCent in the same time, pick one!";
    }

    std::vector<double> ptBinning = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 3.6, 4., 5., 6., 8., 10., 12., 14.};
    std::vector<double> ptBinning_short = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 3.6, 4};
    std::vector<double> ptBinning_diff = {-14.0, -12.0, -10.0, -8.0, -6.0, -5.0, -4.0, -3.6, -3.2, -2.8, -2.4, -2.2, -2.0, -1.8, -1.6, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 3.6, 4., 5., 6., 8., 10., 12., 14.};
    std::vector<double> centBinning = {0., 1., 5., 10., 20., 30., 40., 50., 70., 100.};
    std::vector<double> etaBinning = {-1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec ptAxis_reduced = {ptBinning_short, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec centAxis = {centBinning, "V0M (%)"};
    AxisSpec centralityAxis = {100, 0.0, 100.0, "VT0C (%)"};
    AxisSpec centralityAxis_extended = {105, 0.0, 105.0, "VT0C (%)"};
    AxisSpec etaAxis = {etaBinning, "#eta"};
    AxisSpec PDGBINNING = {20, 0.0, 20.0, "PDG code"};
    AxisSpec ptAxis_diff = {ptBinning_diff, "#it{p}_{T,diff} (GeV/#it{c})"};

    // +++++++++++++++++++++ Data ++++++++++++++++++++++++

    // QA histograms
    spectra_reg.add("histRecVtxZData", "collision z position", HistType::kTH1F, {{200, -20., +20., "z position (cm)"}});
    spectra_reg.add("histTpcSignalData", "Specific energy loss", HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}});
    spectra_reg.add("histTofSignalData", "TOF signal", HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    spectra_reg.add("histDcaVsPtData_particle", "dcaXY vs Pt (particle)", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    spectra_reg.add("histDcaZVsPtData_particle", "dcaZ vs Pt (particle)", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    spectra_reg.add("histDcaVsPtData_antiparticle", "dcaXY vs Pt (antiparticle)", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    spectra_reg.add("histDcaZVsPtData_antiparticle", "dcaZ vs Pt (antiparticle)", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    spectra_reg.add("histTOFm2", "TOF m^2 vs Pt", HistType::kTH2F, {ptAxis, {400, 0.0, 10.0, "m^2"}});
    spectra_reg.add("histNClusterTPC", "Number of Clusters in TPC vs Pt", HistType::kTH2F, {ptAxis, {160, 0.0, 160.0, "nCluster"}});
    spectra_reg.add("histNClusterITS", "Number of Clusters in ITS vs Pt", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    spectra_reg.add("histNClusterITSib", "Number of Clusters in ib of ITS vs Pt", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    spectra_reg.add("histChi2TPC", "chi^2 TPC vs Pt", HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    spectra_reg.add("histChi2ITS", "chi^2 ITS vs Pt", HistType::kTH2F, {ptAxis, {500, 0.0, 50.0, "chi^2"}});
    spectra_reg.add("histCentrality", "Centrality", HistType::kTH1F, {centralityAxis_extended});
    spectra_reg.add("histEtaWithOverFlow", "Pseudorapidity 0 - 105%% centrality", HistType::kTH1F, {etaAxis});
    spectra_reg.add("histEta", "Pseudorapidity with centrality cut", HistType::kTH1F, {etaAxis});
    spectra_reg.add("histEta_cent", "Pseudorapidity vs Centrality", HistType::kTH2F, {centralityAxis_extended, etaAxis});

    // histograms for Proton
    proton_reg.add("histKeepEventData", "skimming histogram (p)", HistType::kTH1F, {{2, -0.5, +1.5, "true: keep event, false: reject event"}});
    proton_reg.add("histTpcSignalData", "Specific energy loss (p)", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}});
    proton_reg.add("histTofSignalData", "TOF signal (p)", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    proton_reg.add("histDcaVsPtData", "dcaXY vs Pt (p)", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    proton_reg.add("histDcaZVsPtData", "dcaZ vs Pt (p)", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    proton_reg.add("histTOFm2", "TOF m^2 vs Pt (p)", HistType::kTH2F, {ptAxis, {400, 0.0, 10.0, "m^2"}});
    proton_reg.add("histTpcNsigmaData", "n-sigma TPC (p)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{p}"}});
    proton_reg.add("histTofNsigmaData", "n-sigma TOF (p)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{p}"}});
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
    aproton_reg.add("histKeepEventData", "skimming histogram (antip)", HistType::kTH1F, {{2, -0.5, +1.5, "true: keep event, false: reject event"}});
    aproton_reg.add("histTpcSignalData", "Specific energy loss (antip)", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}});
    aproton_reg.add("histTofSignalData", "TOF signal (antip)", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    aproton_reg.add("histDcaVsPtData", "dcaXY vs Pt (antip)", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    aproton_reg.add("histDcaZVsPtData", "dcaZ vs Pt (antip)", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    aproton_reg.add("histTOFm2", "TOF m^2 vs Pt (antip)", HistType::kTH2F, {ptAxis, {400, 0.0, 10.0, "m^2"}});
    aproton_reg.add("histTpcNsigmaData", "n-sigma TPC (antip)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{antip}"}});
    aproton_reg.add("histTofNsigmaData", "n-sigma TOF (antip)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{antip}"}});
    aproton_reg.add("histNClusterTPC", "Number of Clusters in TPC vs Pt (antip)", HistType::kTH2F, {ptAxis, {160, 0.0, 160.0, "nCluster"}});
    aproton_reg.add("histNClusterITS", "Number of Clusters in ITS vs Pt (antip)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    aproton_reg.add("histNClusterITSib", "Number of Clusters in ib of ITS vs Pt (antip)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    aproton_reg.add("histChi2TPC", "chi^2 TPC vs Pt (antip)", HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    aproton_reg.add("histChi2ITS", "chi^2 ITS vs Pt (antip)", HistType::kTH2F, {ptAxis, {500, 0.0, 50.0, "chi^2"}});
    aproton_reg.add("histTpcNsigmaData_cent", "n-sigma TPC (antip) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{antip}"}, centralityAxis});
    aproton_reg.add("histTofNsigmaData_cent", "n-sigma TOF (antip) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{antip}"}, centralityAxis});
    aproton_reg.add("histTofm2_cent", "mass^2 TOF (antip) centrality", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{antip}"}, centralityAxis});
    aproton_reg.add("histTpcNsigmaData_eta", "n-sigma TPC (antip) vs eta", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{p}"}, etaAxis});
    aproton_reg.add("histTofNsigmaData_eta", "n-sigma TOF (antip) vs eta", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{p}"}, etaAxis});
    aproton_reg.add("histTofm2_eta", "mass^2 TOF (antip) vs eta", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{p}"}, etaAxis});

    // histograms for Deuterons
    deuteron_reg.add("histKeepEventData", "skimming histogram (d)", HistType::kTH1F, {{2, -0.5, +1.5, "true: keep event, false: reject event"}});
    deuteron_reg.add("histTpcSignalData", "Specific energy loss (d)", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}});
    deuteron_reg.add("histTofSignalData", "TOF signal (d)", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    deuteron_reg.add("histDcaVsPtData", "dcaXY vs Pt (d)", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    deuteron_reg.add("histDcaZVsPtData", "dcaZ vs Pt (d)", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    deuteron_reg.add("histTOFm2", "TOF m^2 vs Pt (d)", HistType::kTH2F, {ptAxis, {400, 0.0, 10.0, "m^2"}});
    deuteron_reg.add("histTpcNsigmaData", "n-sigma TPC (d)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{d}"}});
    deuteron_reg.add("histTofNsigmaData", "n-sigma TOF (d)", HistType::kTH2F, {ptAxis_reduced, {160, -20., +20., "n#sigma_{d}"}});
    deuteron_reg.add("histNClusterTPC", "Number of Clusters in TPC vs Pt (d)", HistType::kTH2F, {ptAxis, {160, 0.0, 160.0, "nCluster"}});
    deuteron_reg.add("histNClusterITS", "Number of Clusters in ITS vs Pt (d)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    deuteron_reg.add("histNClusterITSib", "Number of Clusters in ib of ITS vs Pt (d)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    deuteron_reg.add("histChi2TPC", "chi^2 TPC vs Pt (d)", HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    deuteron_reg.add("histChi2ITS", "chi^2 ITS vs Pt (d)", HistType::kTH2F, {ptAxis, {500, 0.0, 50.0, "chi^2"}});
    deuteron_reg.add("histTpcNsigmaData_cent", "n-sigma TPC (d) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{d}"}, centralityAxis});
    deuteron_reg.add("histTofNsigmaData_cent", "n-sigma TOF (d) centrality", HistType::kTH3F, {ptAxis_reduced, {160, -20., +20., "n#sigma_{d}"}, centralityAxis});
    deuteron_reg.add("histTofm2_cent", "mass^2 TOF (d) centrality", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{d}"}, centralityAxis});
    deuteron_reg.add("histTpcNsigmaData_eta", "n-sigma TPC (d) vs eta", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{p}"}, etaAxis});
    deuteron_reg.add("histTofNsigmaData_eta", "n-sigma TOF (d) vs eta", HistType::kTH3F, {ptAxis_reduced, {160, -20., +20., "n#sigma_{p}"}, etaAxis});
    deuteron_reg.add("histTofm2_eta", "mass^2 TOF (d) vs eta", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{p}"}, etaAxis});

    // histograms for antiDeuterons
    adeuteron_reg.add("histKeepEventData", "skimming histogram (antid)", HistType::kTH1F, {{2, -0.5, +1.5, "true: keep event, false: reject event"}});
    adeuteron_reg.add("histTpcSignalData", "Specific energy loss (antid)", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}});
    adeuteron_reg.add("histTofSignalData", "TOF signal (antid)", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    adeuteron_reg.add("histDcaVsPtData", "dcaXY vs Pt (antid)", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    adeuteron_reg.add("histDcaZVsPtData", "dcaZ vs Pt (antid)", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    adeuteron_reg.add("histTOFm2", "TOF m^2 vs Pt (antid)", HistType::kTH2F, {ptAxis, {400, 0.0, 10.0, "m^2"}});
    adeuteron_reg.add("histTpcNsigmaData", "n-sigma TPC (antid)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{antid}"}});
    adeuteron_reg.add("histTofNsigmaData", "n-sigma TOF (antid)", HistType::kTH2F, {ptAxis_reduced, {160, -20., +20., "n#sigma_{antid}"}});
    adeuteron_reg.add("histNClusterTPC", "Number of Clusters in TPC vs Pt (antid)", HistType::kTH2F, {ptAxis, {160, 0.0, 160.0, "nCluster"}});
    adeuteron_reg.add("histNClusterITS", "Number of Clusters in ITS vs Pt (antid)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    adeuteron_reg.add("histNClusterITSib", "Number of Clusters in ib of ITS vs Pt (antid)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    adeuteron_reg.add("histChi2TPC", "chi^2 TPC vs Pt (antid)", HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    adeuteron_reg.add("histChi2ITS", "chi^2 ITS vs Pt (antid)", HistType::kTH2F, {ptAxis, {500, 0.0, 50.0, "chi^2"}});
    adeuteron_reg.add("histTpcNsigmaData_cent", "n-sigma TPC (antid) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{antid}"}, centralityAxis});
    adeuteron_reg.add("histTofNsigmaData_cent", "n-sigma TOF (antid) centrality", HistType::kTH3F, {ptAxis_reduced, {160, -20., +20., "n#sigma_{antid}"}, centralityAxis});
    adeuteron_reg.add("histTofm2_cent", "mass^2 TOF (antid) centrality", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{antid}"}, centralityAxis});
    adeuteron_reg.add("histTpcNsigmaData_eta", "n-sigma TPC (antid) vs eta", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{p}"}, etaAxis});
    adeuteron_reg.add("histTofNsigmaData_eta", "n-sigma TOF (antid) vs eta", HistType::kTH3F, {ptAxis_reduced, {160, -20., +20., "n#sigma_{p}"}, etaAxis});
    adeuteron_reg.add("histTofm2_eta", "mass^2 TOF (antid) vs eta", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{p}"}, etaAxis});

    // histograms for Triton
    triton_reg.add("histKeepEventData", "skimming histogram (t)", HistType::kTH1F, {{2, -0.5, +1.5, "true: keep event, false: reject event"}});
    triton_reg.add("histTpcSignalData", "Specific energy loss (t)", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}});
    triton_reg.add("histTofSignalData", "TOF signal (t)", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    triton_reg.add("histDcaVsPtData", "dcaXY vs Pt (t)", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    triton_reg.add("histDcaZVsPtData", "dcaZ vs Pt (t)", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    triton_reg.add("histTOFm2", "TOF m^2 vs Pt (t)", HistType::kTH2F, {ptAxis, {400, 0.0, 10.0, "m^2"}});
    triton_reg.add("histTpcNsigmaData", "n-sigma TPC (t)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{t}"}});
    triton_reg.add("histTofNsigmaData", "n-sigma TOF (t)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{t}"}});
    triton_reg.add("histNClusterTPC", "Number of Clusters in TPC vs Pt (t)", HistType::kTH2F, {ptAxis, {160, 0.0, 160.0, "nCluster"}});
    triton_reg.add("histNClusterITS", "Number of Clusters in ITS vs Pt (t)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    triton_reg.add("histNClusterITSib", "Number of Clusters in ib of ITS vs Pt (t)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    triton_reg.add("histChi2TPC", "chi^2 TPC vs Pt (t)", HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    triton_reg.add("histChi2ITS", "chi^2 ITS vs Pt (t)", HistType::kTH2F, {ptAxis, {500, 0.0, 50.0, "chi^2"}});
    triton_reg.add("histTpcNsigmaData_cent", "n-sigma TPC (t) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{t}"}, centralityAxis});
    triton_reg.add("histTofNsigmaData_cent", "n-sigma TOF (t) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{t}"}, centralityAxis});
    triton_reg.add("histTofm2_cent", "mass^2 TOF (t) centrality", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{t}"}, centralityAxis});
    triton_reg.add("histTpcNsigmaData_eta", "n-sigma TPC (t) vs eta", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{p}"}, etaAxis});
    triton_reg.add("histTofNsigmaData_eta", "n-sigma TOF (t) vs eta", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{p}"}, etaAxis});
    triton_reg.add("histTofm2_eta", "mass^2 TOF (t) vs eta", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{p}"}, etaAxis});

    // histograms for antiTriton
    atriton_reg.add("histKeepEventData", "skimming histogram (antit)", HistType::kTH1F, {{2, -0.5, +1.5, "true: keep event, false: reject event"}});
    atriton_reg.add("histTpcSignalData", "Specific energy loss (antit)", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}});
    atriton_reg.add("histTofSignalData", "TOF signal (antit)", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    atriton_reg.add("histDcaVsPtData", "dcaXY vs Pt (antit)", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    atriton_reg.add("histDcaZVsPtData", "dcaZ vs Pt (antit)", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    atriton_reg.add("histTOFm2", "TOF m^2 vs Pt (antit)", HistType::kTH2F, {ptAxis, {400, 0.0, 10.0, "m^2"}});
    atriton_reg.add("histTpcNsigmaData", "n-sigma TPC (antit)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{antit}"}});
    atriton_reg.add("histTofNsigmaData", "n-sigma TOF (antit)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{antit}"}});
    atriton_reg.add("histNClusterTPC", "Number of Clusters in TPC vs Pt (antit)", HistType::kTH2F, {ptAxis, {160, 0.0, 160.0, "nCluster"}});
    atriton_reg.add("histNClusterITS", "Number of Clusters in ITS vs Pt (antit)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    atriton_reg.add("histNClusterITSib", "Number of Clusters in ib of ITS vs Pt (antit)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    atriton_reg.add("histChi2TPC", "chi^2 TPC vs Pt (antit)", HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    atriton_reg.add("histChi2ITS", "chi^2 ITS vs Pt (antit)", HistType::kTH2F, {ptAxis, {500, 0.0, 50.0, "chi^2"}});
    atriton_reg.add("histTpcNsigmaData_cent", "n-sigma TPC (antit) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{antit}"}, centralityAxis});
    atriton_reg.add("histTofNsigmaData_cent", "n-sigma TOF (antit) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{antit}"}, centralityAxis});
    atriton_reg.add("histTofm2_cent", "mass^2 TOF (antit) centrality", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{antit}"}, centralityAxis});
    atriton_reg.add("histTpcNsigmaData_eta", "n-sigma TPC (antit) vs eta", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{p}"}, etaAxis});
    atriton_reg.add("histTofNsigmaData_eta", "n-sigma TOF (antit) vs eta", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{p}"}, etaAxis});
    atriton_reg.add("histTofm2_eta", "mass^2 TOF (antit) vs eta", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{p}"}, etaAxis});

    // histograms for Helium-3
    Helium3_reg.add("histKeepEventData", "skimming histogram (He-3)", HistType::kTH1F, {{2, -0.5, +1.5, "true: keep event, false: reject event"}});
    Helium3_reg.add("histTpcSignalData", "Specific energy loss (He-3)", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}});
    Helium3_reg.add("histTofSignalData", "TOF signal (He-3)", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    Helium3_reg.add("histDcaVsPtData", "dcaXY vs Pt (He-3)", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    Helium3_reg.add("histDcaZVsPtData", "dcaZ vs Pt (He-3)", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    Helium3_reg.add("histTOFm2", "TOF m^2 vs Pt (He-3)", HistType::kTH2F, {ptAxis, {400, 0.0, 10.0, "m^2"}});
    Helium3_reg.add("histTpcNsigmaData", "n-sigma TPC (He-3)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{He-3}"}});
    Helium3_reg.add("histTofNsigmaData", "n-sigma TOF (He-3)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{He-3}"}});
    Helium3_reg.add("histNClusterTPC", "Number of Clusters in TPC vs Pt (He-3)", HistType::kTH2F, {ptAxis, {160, 0.0, 160.0, "nCluster"}});
    Helium3_reg.add("histNClusterITS", "Number of Clusters in ITS vs Pt (He-3)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    Helium3_reg.add("histNClusterITSib", "Number of Clusters in ib of ITS vs Pt (He-3)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    Helium3_reg.add("histChi2TPC", "chi^2 TPC vs Pt (He-3)", HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    Helium3_reg.add("histChi2ITS", "chi^2 ITS vs Pt (He-3)", HistType::kTH2F, {ptAxis, {500, 0.0, 50.0, "chi^2"}});
    Helium3_reg.add("histTpcNsigmaData_cent", "n-sigma TPC (He-3) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{He-3}"}, centralityAxis});
    Helium3_reg.add("histTofNsigmaData_cent", "n-sigma TOF (He-3) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{He-3}"}, centralityAxis});
    Helium3_reg.add("histTofm2_cent", "mass^2 TOF (He-3) centrality", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{He-3}"}, centralityAxis});
    Helium3_reg.add("histTpcNsigmaData_eta", "n-sigma TPC (He-3) vs eta", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{p}"}, etaAxis});
    Helium3_reg.add("histTofNsigmaData_eta", "n-sigma TOF (He-3) vs eta", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{p}"}, etaAxis});
    Helium3_reg.add("histTofm2_eta", "mass^2 TOF (He-3) vs eta", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{p}"}, etaAxis});

    // histograms for antiHelium-3
    aHelium3_reg.add("histKeepEventData", "skimming histogram (antiHe-3)", HistType::kTH1F, {{2, -0.5, +1.5, "true: keep event, false: reject event"}});
    aHelium3_reg.add("histTpcSignalData", "Specific energy loss (antiHe-3)", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}});
    aHelium3_reg.add("histTofSignalData", "TOF signal (antiHe-3)", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    aHelium3_reg.add("histDcaVsPtData", "dcaXY vs Pt (antiHe-3)", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    aHelium3_reg.add("histDcaZVsPtData", "dcaZ vs Pt (antiHe-3)", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    aHelium3_reg.add("histTOFm2", "TOF m^2 vs Pt (antiHe-3)", HistType::kTH2F, {ptAxis, {400, 0.0, 10.0, "m^2"}});
    aHelium3_reg.add("histTpcNsigmaData", "n-sigma TPC (antiHe-3)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{antiHe-3}"}});
    aHelium3_reg.add("histTofNsigmaData", "n-sigma TOF (antiHe-3)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{antiHe-3}"}});
    aHelium3_reg.add("histNClusterTPC", "Number of Clusters in TPC vs Pt (antiHe-3)", HistType::kTH2F, {ptAxis, {160, 0.0, 160.0, "nCluster"}});
    aHelium3_reg.add("histNClusterITS", "Number of Clusters in ITS vs Pt (antiHe-3)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    aHelium3_reg.add("histNClusterITSib", "Number of Clusters in ib of ITS vs Pt (antiHe-3)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    aHelium3_reg.add("histChi2TPC", "chi^2 TPC vs Pt (antiHe-3)", HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    aHelium3_reg.add("histChi2ITS", "chi^2 ITS vs Pt (antiHe-3)", HistType::kTH2F, {ptAxis, {500, 0.0, 50.0, "chi^2"}});
    aHelium3_reg.add("histTpcNsigmaData_cent", "n-sigma TPC (antiHe-3) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{antiHe-3}"}, centralityAxis});
    aHelium3_reg.add("histTofNsigmaData_cent", "n-sigma TOF (antiHe-3) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{antiHe-3}"}, centralityAxis});
    aHelium3_reg.add("histTofm2_cent", "mass^2 TOF (antiHe-3) centrality", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{antiHe-3}"}, centralityAxis});
    aHelium3_reg.add("histTpcNsigmaData_eta", "n-sigma TPC (antiHe-3) vs eta", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{p}"}, etaAxis});
    aHelium3_reg.add("histTofNsigmaData_eta", "n-sigma TOF (antiHe-3) vs eta", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{p}"}, etaAxis});
    aHelium3_reg.add("histTofm2_eta", "mass^2 TOF (antiHe-3) vs eta", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{p}"}, etaAxis});

    // histograms for Helium-4 (alpha)
    Helium4_reg.add("histKeepEventData", "skimming histogram (He-4)", HistType::kTH1F, {{2, -0.5, +1.5, "true: keep event, false: reject event"}});
    Helium4_reg.add("histTpcSignalData", "Specific energy loss (He-4)", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}});
    Helium4_reg.add("histTofSignalData", "TOF signal (He-4)", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    Helium4_reg.add("histDcaVsPtData", "dcaXY vs Pt (He-4)", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    Helium4_reg.add("histDcaZVsPtData", "dcaZ vs Pt (He-4)", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    Helium4_reg.add("histTOFm2", "TOF m^2 vs Pt (He-4)", HistType::kTH2F, {ptAxis, {400, 0.0, 10.0, "m^2"}});
    Helium4_reg.add("histTpcNsigmaData", "n-sigma TPC (He-4)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{He-4}"}});
    Helium4_reg.add("histTofNsigmaData", "n-sigma TOF (He-4)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{He-4}"}});
    Helium4_reg.add("histNClusterTPC", "Number of Clusters in TPC vs Pt (He-4)", HistType::kTH2F, {ptAxis, {160, 0.0, 160.0, "nCluster"}});
    Helium4_reg.add("histNClusterITS", "Number of Clusters in ITS vs Pt (He-4)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    Helium4_reg.add("histNClusterITSib", "Number of Clusters in ib of ITS vs Pt (He-4)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    Helium4_reg.add("histChi2TPC", "chi^2 TPC vs Pt (He-4)", HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    Helium4_reg.add("histChi2ITS", "chi^2 ITS vs Pt (He-4)", HistType::kTH2F, {ptAxis, {500, 0.0, 50.0, "chi^2"}});
    Helium4_reg.add("histTpcNsigmaData_cent", "n-sigma TPC (He-4) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{He-4}"}, centralityAxis});
    Helium4_reg.add("histTofNsigmaData_cent", "n-sigma TOF (He-4) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{He-4}"}, centralityAxis});
    Helium4_reg.add("histTofm2_cent", "mass^2 TOF (He-4) centrality", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{He-4}"}, centralityAxis});
    Helium4_reg.add("histTpcNsigmaData_eta", "n-sigma TPC (He-4) vs eta", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{p}"}, etaAxis});
    Helium4_reg.add("histTofNsigmaData_eta", "n-sigma TOF (He-4) vs eta", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{p}"}, etaAxis});
    Helium4_reg.add("histTofm2_eta", "mass^2 TOF (He-4) vs eta", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{p}"}, etaAxis});

    // histograms for antiHelium-4 (alpha)
    aHelium4_reg.add("histKeepEventData", "skimming histogram (antiHe-4)", HistType::kTH1F, {{2, -0.5, +1.5, "true: keep event, false: reject event"}});
    aHelium4_reg.add("histTpcSignalData", "Specific energy loss (antiHe-4)", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {5000, 0, 5000, "d#it{E} / d#it{X} (a. u.)"}});
    aHelium4_reg.add("histTofSignalData", "TOF signal (antiHe-4)", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    aHelium4_reg.add("histDcaVsPtData", "dcaXY vs Pt (antiHe-4)", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    aHelium4_reg.add("histDcaZVsPtData", "dcaZ vs Pt (antiHe-4)", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    aHelium4_reg.add("histTOFm2", "TOF m^2 vs Pt (antiHe-4)", HistType::kTH2F, {ptAxis, {400, 0.0, 10.0, "m^2"}});
    aHelium4_reg.add("histTpcNsigmaData", "n-sigma TPC (antiHe-4)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{antiHe-4}"}});
    aHelium4_reg.add("histTofNsigmaData", "n-sigma TOF (antiHe-4)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{antiHe-4}"}});
    aHelium4_reg.add("histNClusterTPC", "Number of Clusters in TPC vs Pt (antiHe-4)", HistType::kTH2F, {ptAxis, {160, 0.0, 160.0, "nCluster"}});
    aHelium4_reg.add("histNClusterITS", "Number of Clusters in ITS vs Pt (antiHe-4)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    aHelium4_reg.add("histNClusterITSib", "Number of Clusters in ib of ITS vs Pt (antiHe-4)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    aHelium4_reg.add("histChi2TPC", "chi^2 TPC vs Pt (antiHe-4)", HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    aHelium4_reg.add("histChi2ITS", "chi^2 ITS vs Pt (antiHe-4)", HistType::kTH2F, {ptAxis, {500, 0.0, 50.0, "chi^2"}});
    aHelium4_reg.add("histTpcNsigmaData_cent", "n-sigma TPC (antiHe-4) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{antiHe-4}"}, centralityAxis});
    aHelium4_reg.add("histTofNsigmaData_cent", "n-sigma TOF (antiHe-4) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{antiHe-4}"}, centralityAxis});
    aHelium4_reg.add("histTofm2_cent", "mass^2 TOF (antiHe-4) centrality", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{antiHe-4}"}, centralityAxis});
    aHelium4_reg.add("histTpcNsigmaData_eta", "n-sigma TPC (antiHe-4) vs eta", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{p}"}, etaAxis});
    aHelium4_reg.add("histTofNsigmaData_eta", "n-sigma TOF (antiHe-4) vs eta", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{p}"}, etaAxis});
    aHelium4_reg.add("histTofm2_eta", "mass^2 TOF (antiHe-4) vs eta", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{p}"}, etaAxis});

    // +++++++++++++++++++++ MC ++++++++++++++++++++++++++

    // MC truth
    MC_truth_reg.add("histPhi", "#phi", HistType::kTH2F, {{100, 0., 2. * TMath::Pi()}, PDGBINNING});
    MC_truth_reg.add("histEta", "#eta", HistType::kTH2F, {{102, -2.01, 2.01}, PDGBINNING});
    MC_truth_reg.add("histPt", "p_{t}", HistType::kTH2F, {ptAxis, PDGBINNING});
    MC_truth_reg.add("histRecVtxMC", "MC collision z position", HistType::kTH1F, {{400, -40., +40., "z position (cm)"}});

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
    MC_recon_reg.add("histTPCFractionSharedCls", "Fraction of shared TPC clusters", HistType::kTH2F, {{100, -2.0, 2.0, "Shared Cls"}, PDGBINNING});
    MC_recon_reg.add("histTpcNsigmaDataPr", "TPC nSigma (p)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{p}"}});
    MC_recon_reg.add("histTpcNsigmaDataaPr", "TPC nSigma (antip)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{p}"}});
    MC_recon_reg.add("histTpcNsigmaDataDe", "TPC nSigma (d)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{d}"}});
    MC_recon_reg.add("histTpcNsigmaDataaDe", "TPC nSigma (antid)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{d}"}});
    MC_recon_reg.add("histTpcNsigmaDataTr", "TPC nSigma (t)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{t}"}});
    MC_recon_reg.add("histTpcNsigmaDataaTr", "TPC nSigma (antit)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{t}"}});
    MC_recon_reg.add("histTpcNsigmaDataHe", "TPC nSigma (He-3)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{He-3}"}});
    MC_recon_reg.add("histTpcNsigmaDataaHe", "TPC nSigma (antiHe-3)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{He-3}"}});
    MC_recon_reg.add("histTpcNsigmaDataAl", "TPC nSigma (He-4)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{He-4}"}});
    MC_recon_reg.add("histTpcNsigmaDataaAl", "TPC nSigma (antiHe-4)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{He-4}"}});
    MC_recon_reg.add("histTofNsigmaDataPr", "TOF nSigma (p)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{p}"}});
    MC_recon_reg.add("histTofNsigmaDataaPr", "TOF nSigma (antip)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{p}"}});
    MC_recon_reg.add("histTofNsigmaDataDe", "TOF nSigma (d)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{d}"}});
    MC_recon_reg.add("histTofNsigmaDataaDe", "TOF nSigma (antid)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{d}"}});
    MC_recon_reg.add("histTofNsigmaDataTr", "TOF nSigma (t)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{t}"}});
    MC_recon_reg.add("histTofNsigmaDataaTr", "TOF nSigma (antit)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{t}"}});
    MC_recon_reg.add("histTofNsigmaDataHe", "TOF nSigma (He-3)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{He-3}"}});
    MC_recon_reg.add("histTofNsigmaDataaHe", "TOF nSigma (antiHe-3)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{He-3}"}});
    MC_recon_reg.add("histTofNsigmaDataAl", "TOF nSigma (He-4)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{He-4}"}});
    MC_recon_reg.add("histTofNsigmaDataaAl", "TOF nSigma (antiHe-4)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{He-4}"}});

    // MC diff (truth - reconstructed)
    MC_recon_diff_reg.add("histPhiDiff", "MC t", HistType::kTH2F, {ptAxis_diff, PDGBINNING});
    MC_recon_diff_reg.add("histEtaDiff", "MC t", HistType::kTH2F, {ptAxis_diff, PDGBINNING});
    MC_recon_diff_reg.add("histPtDiff", "MC t", HistType::kTH2F, {ptAxis_diff, PDGBINNING});
  }

  // Configurables
  Configurable<float> yMin{"yMin", -0.5, "Maximum rapidity"};
  Configurable<float> yMax{"yMax", 0.5, "Minimum rapidity"};
  Configurable<float> pTmin{"pTmin", 0.1f, "min pT"};
  Configurable<float> pTmax{"pTmax", 1e+10f, "max pT"};

  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<float> nsigmacutLow{"nsigmacutLow", -3.0, "Value of the Nsigma cut"};
  Configurable<float> nsigmacutHigh{"nsigmacutHigh", +3.0, "Value of the Nsigma cut"};
  Configurable<float> minCentrality{"minCentrality", 0.0, "min Centrality used"};
  Configurable<float> maxCentrality{"maxCentrality", 80.0, "max Centrality used"};
  Configurable<bool> enable_Centrality_cut_global{"enable_Centrality_cut_global", true, "use Centrality cut"};

  // Replacement for globalTrack filter
  Configurable<float> minReqClusterITS{"minReqClusterITS", 1.0, "min number of clusters required in ITS"};                            // ITS_nCls
  Configurable<float> minReqClusterITSib{"minReqClusterITSib", 1.0, "min number of clusters required in ITS inner barrel"};           // ITS_nCls
  Configurable<float> minTPCnClsFound{"minTPCnClsFound", 0.0f, "min number of crossed rows TPC"};                                     // TPC_nCls_found
  Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 70.0f, "min number of crossed rows TPC"};                              // TPC_nCls_crossed_Rows
  Configurable<float> minRatioCrossedRowsTPC{"minRatioCrossedRowsTPC", 0.8f, "min ratio of crossed rows over findable clusters TPC"}; // TPC_crossed_Rows_over_findable_Cls_min
  Configurable<float> maxRatioCrossedRowsTPC{"maxRatioCrossedRowsTPC", 1.5f, "max ratio of crossed rows over findable clusters TPC"}; // TPC_crossed_Rows_over_findable_Cls_max
  Configurable<float> maxChi2ITS{"maxChi2ITS", 36.0f, "max chi2 per cluster ITS"};
  Configurable<float> maxChi2TPC{"maxChi2TPC", 4.0f, "max chi2 per cluster TPC"};
  Configurable<float> maxDCA_XY{"maxDCA_XY", 0.5f, "max DCA to vertex xy"};
  Configurable<float> maxDCA_Z{"maxDCA_Z", 2.0f, "max DCA to vertex z"};
  Configurable<int> lastRequiredTrdCluster{"lastRequiredTrdCluster", 5, "Last cluster to required in TRD for track selection. -1 does not require any TRD cluster"};

  Configurable<bool> event_selection_sel8{"event_selection_sel8", true, "Enable sel8 event selection"};

  Configurable<bool> enable_PVcontributor_global{"enable_PVcontributor_global", true, "is PV contributor (global)"};
  Configurable<bool> enable_PVcontributor_proton{"enable_PVcontributor_proton", true, "is PV contributor (global)"};
  Configurable<bool> enable_PVcontributor_antiproton{"enable_PVcontributor_antiproton", true, "is PV contributor (global)"};
  Configurable<bool> enable_PVcontributor_deuteron{"enable_PVcontributor_deuteron", true, "is PV contributor (global)"};
  Configurable<bool> enable_PVcontributor_antideuteron{"enable_PVcontributor_antideuteron", true, "is PV contributor (global)"};
  Configurable<bool> enable_PVcontributor_triton{"enable_PVcontributor_triton", true, "is PV contributor (global)"};
  Configurable<bool> enable_PVcontributor_antitriton{"enable_PVcontributor_antitriton", true, "is PV contributor (global)"};
  Configurable<bool> enable_PVcontributor_Helium3{"enable_PVcontributor_Helium3", true, "is PV contributor (global)"};
  Configurable<bool> enable_PVcontributor_antiHelium3{"enable_PVcontributor_antiHelium3", true, "is PV contributor (global)"};
  Configurable<bool> enable_PVcontributor_Helium4{"enable_PVcontributor_Helium4", true, "is PV contributor (global)"};
  Configurable<bool> enable_PVcontributor_antiHelium4{"enable_PVcontributor_antiHelium4", true, "is PV contributor (global)"};

  template <typename CollisionType, typename TracksType>
  void fillHistograms(const CollisionType& event, const TracksType& tracks)
  {

    // collision process loop
    bool keepEvent_p = kFALSE;
    bool keepEvent_d = kFALSE;
    bool keepEvent_t = kFALSE;
    bool keepEvent_He3 = kFALSE;
    bool keepEvent_He4 = kFALSE;

    bool keepEvent_antip = kFALSE;
    bool keepEvent_antid = kFALSE;
    bool keepEvent_antit = kFALSE;
    bool keepEvent_antiHe3 = kFALSE;
    bool keepEvent_antiHe4 = kFALSE;

    if (event_selection_sel8 && event.sel8()) {
      spectra_reg.fill(HIST("histRecVtxZData"), event.posZ());
    }

    if (!event_selection_sel8) {
      spectra_reg.fill(HIST("histRecVtxZData"), event.posZ());
    }

    for (auto track : tracks) { // start loop over tracks

      if (event_selection_sel8 && !event.sel8()) {
        continue;
      }

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

      if (TMath::Abs(track.dcaXY()) > maxDCA_XY || TMath::Abs(track.dcaZ()) > maxDCA_Z) {
        continue;
      }

      if (TPCnumberClsFound < minTPCnClsFound || TPC_nCls_Crossed_Rows < minNCrossedRowsTPC || RatioCrossedRowsOverFindableTPC < minRatioCrossedRowsTPC || RatioCrossedRowsOverFindableTPC > maxRatioCrossedRowsTPC || Chi2perClusterTPC > maxChi2TPC || Chi2perClusterITS > maxChi2ITS || !(track.passedTPCRefit()) || !(track.passedITSRefit()) || (track.itsNClsInnerBarrel()) < minReqClusterITSib || (track.itsNCls()) < minReqClusterITS) {
        continue;
      }
      if (enable_PVcontributor_global && !(track.isPVContributor())) {
        continue;
      }
      if (TMath::Abs(track.dcaXY()) > maxDCA_XY || TMath::Abs(track.dcaZ()) > maxDCA_Z) {
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

      // fill QA histograms
      float nSigmaProton = track.tpcNSigmaPr();
      float nSigmaDeut = track.tpcNSigmaDe();
      float nSigmaTriton = track.tpcNSigmaTr();
      float nSigmaHe3 = track.tpcNSigmaHe();
      float nSigmaHe4 = track.tpcNSigmaAl();

      spectra_reg.fill(HIST("histTpcSignalData"), track.pt() * track.sign(), track.tpcSignal());
      spectra_reg.fill(HIST("histNClusterTPC"), track.pt(), track.tpcNClsCrossedRows());
      spectra_reg.fill(HIST("histNClusterITS"), track.pt(), track.itsNCls());
      spectra_reg.fill(HIST("histNClusterITSib"), track.pt(), track.itsNClsInnerBarrel());
      spectra_reg.fill(HIST("histChi2TPC"), track.pt(), track.tpcChi2NCl());
      spectra_reg.fill(HIST("histChi2ITS"), track.pt(), track.itsChi2NCl());

      if (track.sign() > 0) {
        proton_reg.fill(HIST("histTpcNsigmaData"), track.pt(), nSigmaProton);
        deuteron_reg.fill(HIST("histTpcNsigmaData"), track.pt(), nSigmaDeut);
        triton_reg.fill(HIST("histTpcNsigmaData"), track.pt(), nSigmaTriton);
        Helium3_reg.fill(HIST("histTpcNsigmaData"), track.pt() * 2.0, nSigmaHe3);
        Helium4_reg.fill(HIST("histTpcNsigmaData"), track.pt() * 2.0, nSigmaHe4);

        //  fill TOF m^2 histogram
        if (track.hasTOF()) {

          if (track.hasTRD() && (lastRequiredTrdCluster > 0)) {
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

          spectra_reg.fill(HIST("histTOFm2"), track.pt(), TOFmass2);
        }
      }

      if (track.sign() < 0) {
        aproton_reg.fill(HIST("histTpcNsigmaData"), track.pt(), nSigmaProton);
        adeuteron_reg.fill(HIST("histTpcNsigmaData"), track.pt(), nSigmaDeut);
        atriton_reg.fill(HIST("histTpcNsigmaData"), track.pt(), nSigmaTriton);
        aHelium3_reg.fill(HIST("histTpcNsigmaData"), track.pt() * 2.0, nSigmaHe3);
        aHelium4_reg.fill(HIST("histTpcNsigmaData"), track.pt() * 2.0, nSigmaHe4);

        // fill TOF m^2 histogram
        if (track.hasTOF()) {

          if (track.hasTRD() && (lastRequiredTrdCluster > 0)) {
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

          spectra_reg.fill(HIST("histTOFm2"), track.pt(), TOFmass2);
        }
      }

      //**************   check offline-trigger (skimming) condidition Proton   *******************

      if (nSigmaProton > nsigmacutLow && nSigmaProton < nsigmacutHigh) {

        if (enable_PVcontributor_proton && !(track.isPVContributor())) {
          continue;
        }

        if (track.sign() > 0) {
          keepEvent_p = kTRUE;

          proton_reg.fill(HIST("histDcaVsPtData"), track.pt(), track.dcaXY());
          proton_reg.fill(HIST("histDcaZVsPtData"), track.pt(), track.dcaZ());
          proton_reg.fill(HIST("histTpcSignalData"), track.pt(), track.tpcSignal());
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
              if (lastLayer < lastRequiredTrdCluster) {
                continue;
              }
            }

            Float_t TOFmass2 = ((track.mass()) * (track.mass()));
            Float_t beta = track.beta();

            proton_reg.fill(HIST("histTOFm2"), track.pt(), TOFmass2);
            proton_reg.fill(HIST("histTofSignalData"), track.pt(), beta);
            proton_reg.fill(HIST("histTofNsigmaData"), track.pt(), track.tofNSigmaPr());
          }
        }

        if (enable_PVcontributor_antiproton && !(track.isPVContributor())) {
          continue;
        }

        if (track.sign() < 0) {
          keepEvent_antip = kTRUE;

          aproton_reg.fill(HIST("histDcaVsPtData"), track.pt(), track.dcaXY());
          aproton_reg.fill(HIST("histDcaZVsPtData"), track.pt(), track.dcaZ());
          aproton_reg.fill(HIST("histTpcSignalData"), track.pt(), track.tpcSignal());
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
              if (lastLayer < lastRequiredTrdCluster) {
                continue;
              }
            }

            Float_t TOFmass2 = ((track.mass()) * (track.mass()));
            Float_t beta = track.beta();

            aproton_reg.fill(HIST("histTOFm2"), track.pt(), TOFmass2);
            aproton_reg.fill(HIST("histTofSignalData"), track.pt(), beta);
            aproton_reg.fill(HIST("histTofNsigmaData"), track.pt(), track.tofNSigmaPr());
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
            if (lastLayer < lastRequiredTrdCluster) {
              continue;
            }
          }

          spectra_reg.fill(HIST("histTofSignalData"), track.pt() * track.sign(), track.beta());
        }
      }

      //**************   check offline-trigger (skimming) condidition Deuteron   *******************

      if (nSigmaDeut > nsigmacutLow && nSigmaDeut < nsigmacutHigh) {

        if (enable_PVcontributor_deuteron && !(track.isPVContributor())) {
          continue;
        }

        if (track.sign() > 0) {
          keepEvent_d = kTRUE;

          deuteron_reg.fill(HIST("histDcaVsPtData"), track.pt(), track.dcaXY());
          deuteron_reg.fill(HIST("histDcaZVsPtData"), track.pt(), track.dcaZ());
          deuteron_reg.fill(HIST("histTpcSignalData"), track.pt(), track.tpcSignal());
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
              if (lastLayer < lastRequiredTrdCluster) {
                continue;
              }
            }

            Float_t TOFmass2 = ((track.mass()) * (track.mass()));
            Float_t beta = track.beta();

            deuteron_reg.fill(HIST("histTOFm2"), track.pt(), TOFmass2);
            deuteron_reg.fill(HIST("histTofSignalData"), track.pt(), beta);
            deuteron_reg.fill(HIST("histTofNsigmaData"), track.pt(), track.tofNSigmaDe());
          }
        }

        if (enable_PVcontributor_antideuteron && !(track.isPVContributor())) {
          continue;
        }

        if (track.sign() < 0) {
          keepEvent_antid = kTRUE;

          adeuteron_reg.fill(HIST("histDcaVsPtData"), track.pt(), track.dcaXY());
          adeuteron_reg.fill(HIST("histDcaZVsPtData"), track.pt(), track.dcaZ());
          adeuteron_reg.fill(HIST("histTpcSignalData"), track.pt(), track.tpcSignal());
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
              if (lastLayer < lastRequiredTrdCluster) {
                continue;
              }
            }

            Float_t TOFmass2 = ((track.mass()) * (track.mass()));
            Float_t beta = track.beta();

            adeuteron_reg.fill(HIST("histTOFm2"), track.pt(), TOFmass2);
            adeuteron_reg.fill(HIST("histTofSignalData"), track.pt(), beta);
            adeuteron_reg.fill(HIST("histTofNsigmaData"), track.pt(), track.tofNSigmaDe());
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
            if (lastLayer < lastRequiredTrdCluster) {
              continue;
            }
          }

          spectra_reg.fill(HIST("histTofSignalData"), track.pt() * track.sign(), track.beta());
        }
      }

      //**************   check offline-trigger (skimming) condidition Triton   *******************

      if (nSigmaTriton > nsigmacutLow && nSigmaTriton < nsigmacutHigh) {

        if (enable_PVcontributor_triton && !(track.isPVContributor())) {
          continue;
        }

        if (track.sign() > 0) {
          keepEvent_t = kTRUE;

          triton_reg.fill(HIST("histDcaVsPtData"), track.pt(), track.dcaXY());
          triton_reg.fill(HIST("histDcaZVsPtData"), track.pt(), track.dcaZ());
          triton_reg.fill(HIST("histTpcSignalData"), track.pt(), track.tpcSignal());
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
              if (lastLayer < lastRequiredTrdCluster) {
                continue;
              }
            }

            Float_t TOFmass2 = ((track.mass()) * (track.mass()));
            Float_t beta = track.beta();

            triton_reg.fill(HIST("histTOFm2"), track.pt(), TOFmass2);
            triton_reg.fill(HIST("histTofSignalData"), track.pt(), beta);
            triton_reg.fill(HIST("histTofNsigmaData"), track.pt(), track.tofNSigmaTr());
          }
        }

        if (enable_PVcontributor_antitriton && !(track.isPVContributor())) {
          continue;
        }

        if (track.sign() < 0) {
          keepEvent_antit = kTRUE;

          atriton_reg.fill(HIST("histDcaVsPtData"), track.pt(), track.dcaXY());
          atriton_reg.fill(HIST("histDcaZVsPtData"), track.pt(), track.dcaZ());
          atriton_reg.fill(HIST("histTpcSignalData"), track.pt(), track.tpcSignal());
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
              if (lastLayer < lastRequiredTrdCluster) {
                continue;
              }
            }

            Float_t TOFmass2 = ((track.mass()) * (track.mass()));
            Float_t beta = track.beta();

            atriton_reg.fill(HIST("histTOFm2"), track.pt(), TOFmass2);
            atriton_reg.fill(HIST("histTofSignalData"), track.pt(), beta);
            atriton_reg.fill(HIST("histTofNsigmaData"), track.pt(), track.tofNSigmaTr());
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
            if (lastLayer < lastRequiredTrdCluster) {
              continue;
            }
          }

          spectra_reg.fill(HIST("histTofSignalData"), track.pt() * track.sign(), track.beta());
        }
      }

      //**************   check offline-trigger (skimming) condidition Helium-3   *******************

      if (nSigmaHe3 > nsigmacutLow && nSigmaHe3 < nsigmacutHigh) {

        if (enable_PVcontributor_Helium3 && !(track.isPVContributor())) {
          continue;
        }

        if (track.sign() > 0) {
          keepEvent_He3 = kTRUE;

          Helium3_reg.fill(HIST("histDcaVsPtData"), track.pt() * 2.0, track.dcaXY());
          Helium3_reg.fill(HIST("histDcaZVsPtData"), track.pt() * 2.0, track.dcaZ());
          Helium3_reg.fill(HIST("histTpcSignalData"), track.pt() * 2.0, track.tpcSignal());
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
              if (lastLayer < lastRequiredTrdCluster) {
                continue;
              }
            }

            Float_t TOFmass2 = ((track.mass()) * (track.mass()));
            Float_t beta = track.beta();

            Helium3_reg.fill(HIST("histTOFm2"), track.pt() * 2.0, TOFmass2);
            Helium3_reg.fill(HIST("histTofSignalData"), track.pt() * 2.0, beta);
            Helium3_reg.fill(HIST("histTofNsigmaData"), track.pt() * 2.0, track.tofNSigmaHe());
          }
        }

        if (enable_PVcontributor_antiHelium3 && !(track.isPVContributor())) {
          continue;
        }

        if (track.sign() < 0) {
          keepEvent_antiHe3 = kTRUE;
          aHelium3_reg.fill(HIST("histDcaVsPtData"), track.pt() * 2.0, track.dcaXY());
          aHelium3_reg.fill(HIST("histDcaZVsPtData"), track.pt() * 2.0, track.dcaZ());
          aHelium3_reg.fill(HIST("histTpcSignalData"), track.pt() * 2.0, track.tpcSignal());
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
              if (lastLayer < lastRequiredTrdCluster) {
                continue;
              }
            }

            Float_t TOFmass2 = ((track.mass()) * (track.mass()));
            Float_t beta = track.beta();

            aHelium3_reg.fill(HIST("histTOFm2"), track.pt() * 2.0, TOFmass2);
            aHelium3_reg.fill(HIST("histTofSignalData"), track.pt() * 2.0, beta);
            aHelium3_reg.fill(HIST("histTofNsigmaData"), track.pt() * 2.0, track.tofNSigmaHe());
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
            if (lastLayer < lastRequiredTrdCluster) {
              continue;
            }
          }

          spectra_reg.fill(HIST("histTofSignalData"), track.pt() * 2.0 * track.sign(), track.beta());
        }
      }

      //**************   check offline-trigger (skimming) condidition Helium-4   *******************

      if (nSigmaHe4 > nsigmacutLow && nSigmaHe4 < nsigmacutHigh) {

        if (enable_PVcontributor_Helium4 && !(track.isPVContributor())) {
          continue;
        }

        if (track.sign() > 0) {
          keepEvent_He4 = kTRUE;

          Helium4_reg.fill(HIST("histDcaVsPtData"), track.pt() * 2.0, track.dcaXY());
          Helium4_reg.fill(HIST("histDcaZVsPtData"), track.pt() * 2.0, track.dcaZ());
          Helium4_reg.fill(HIST("histTpcSignalData"), track.pt() * 2.0, track.tpcSignal());
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
              if (lastLayer < lastRequiredTrdCluster) {
                continue;
              }
            }

            Float_t TOFmass2 = ((track.mass()) * (track.mass()));
            Float_t beta = track.beta();

            Helium4_reg.fill(HIST("histTOFm2"), track.pt() * 2.0, TOFmass2);
            Helium4_reg.fill(HIST("histTofSignalData"), track.pt() * 2.0, beta);
            Helium4_reg.fill(HIST("histTofNsigmaData"), track.pt() * 2.0, track.tofNSigmaAl());
          }
        }

        if (enable_PVcontributor_antiHelium4 && !(track.isPVContributor())) {
          continue;
        }

        if (track.sign() < 0) {
          keepEvent_antiHe4 = kTRUE;
          aHelium4_reg.fill(HIST("histDcaVsPtData"), track.pt() * 2.0, track.dcaXY());
          aHelium4_reg.fill(HIST("histDcaZVsPtData"), track.pt() * 2.0, track.dcaZ());
          aHelium4_reg.fill(HIST("histTpcSignalData"), track.pt() * 2.0, track.tpcSignal());
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
              if (lastLayer < lastRequiredTrdCluster) {
                continue;
              }
            }

            Float_t TOFmass2 = ((track.mass()) * (track.mass()));
            Float_t beta = track.beta();

            aHelium4_reg.fill(HIST("histTOFm2"), track.pt() * 2.0, TOFmass2);
            aHelium4_reg.fill(HIST("histTofSignalData"), track.pt() * 2.0, beta);
            aHelium4_reg.fill(HIST("histTofNsigmaData"), track.pt() * 2.0, track.tofNSigmaAl());
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
            if (lastLayer < lastRequiredTrdCluster) {
              continue;
            }
          }

          spectra_reg.fill(HIST("histTofSignalData"), track.pt() * 2.0 * track.sign(), track.beta());
        }
      }

    } // end loop over tracks

    // fill trigger (skimming) results
    proton_reg.fill(HIST("histKeepEventData"), keepEvent_p);
    aproton_reg.fill(HIST("histKeepEventData"), keepEvent_antip);
    deuteron_reg.fill(HIST("histKeepEventData"), keepEvent_d);
    adeuteron_reg.fill(HIST("histKeepEventData"), keepEvent_antid);
    triton_reg.fill(HIST("histKeepEventData"), keepEvent_t);
    atriton_reg.fill(HIST("histKeepEventData"), keepEvent_antit);
    Helium3_reg.fill(HIST("histKeepEventData"), keepEvent_He3);
    aHelium3_reg.fill(HIST("histKeepEventData"), keepEvent_antiHe3);
    Helium4_reg.fill(HIST("histKeepEventData"), keepEvent_He4);
    aHelium4_reg.fill(HIST("histKeepEventData"), keepEvent_antiHe4);
  }

  //****************************************************************************************************

  template <typename CollisionType, typename TracksType>
  void fillCentHistorgrams(const CollisionType& event, const TracksType& tracks)
  {

    if (event_selection_sel8 && event.sel8()) {
      spectra_reg.fill(HIST("histCentrality"), event.centFT0C());
    }

    if (!event_selection_sel8) {
      spectra_reg.fill(HIST("histCentrality"), event.centFT0C());
    }

    for (auto track : tracks) { // start loop over tracks

      if (event_selection_sel8 && !event.sel8()) {
        continue;
      }

      if (enable_Centrality_cut_global && (event.centFT0C() < minCentrality) && (event.centFT0C() > maxCentrality)) {
        continue;
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

      // track cuts
      if (enable_PVcontributor_global && !(track.isPVContributor())) {
        continue;
      }

      if (TPCnumberClsFound < minTPCnClsFound || TPC_nCls_Crossed_Rows < minNCrossedRowsTPC || RatioCrossedRowsOverFindableTPC < minRatioCrossedRowsTPC || RatioCrossedRowsOverFindableTPC > maxRatioCrossedRowsTPC || Chi2perClusterTPC > maxChi2TPC || Chi2perClusterITS > maxChi2ITS || !(track.passedTPCRefit()) || !(track.passedITSRefit()) || (track.itsNClsInnerBarrel()) < minReqClusterITSib || (track.itsNCls()) < minReqClusterITS) {
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

      // fill 3D centrality histograms
      if (track.sign() > 0) {

        proton_reg.fill(HIST("histTpcNsigmaData_cent"), track.pt(), track.tpcNSigmaPr(), event.centFT0C());
        deuteron_reg.fill(HIST("histTpcNsigmaData_cent"), track.pt(), track.tpcNSigmaDe(), event.centFT0C());
        triton_reg.fill(HIST("histTpcNsigmaData_cent"), track.pt(), track.tpcNSigmaTr(), event.centFT0C());
        Helium3_reg.fill(HIST("histTpcNsigmaData_cent"), track.pt() * 2.0, track.tpcNSigmaHe(), event.centFT0C());
        Helium4_reg.fill(HIST("histTpcNsigmaData_cent"), track.pt() * 2.0, track.tpcNSigmaAl(), event.centFT0C());

        if ((event.centFT0C() > minCentrality) && (event.centFT0C() < maxCentrality)) {
          proton_reg.fill(HIST("histTpcNsigmaData_eta"), track.pt(), track.tpcNSigmaPr(), track.eta());
          deuteron_reg.fill(HIST("histTpcNsigmaData_eta"), track.pt(), track.tpcNSigmaDe(), track.eta());
          triton_reg.fill(HIST("histTpcNsigmaData_eta"), track.pt(), track.tpcNSigmaTr(), track.eta());
          Helium3_reg.fill(HIST("histTpcNsigmaData_eta"), track.pt() * 2.0, track.tpcNSigmaHe(), track.eta());
          Helium4_reg.fill(HIST("histTpcNsigmaData_eta"), track.pt() * 2.0, track.tpcNSigmaAl(), track.eta());
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
            if (lastLayer < lastRequiredTrdCluster) {
              continue;
            }
          }

          proton_reg.fill(HIST("histTofNsigmaData_cent"), track.pt(), track.tofNSigmaPr(), event.centFT0C());
          deuteron_reg.fill(HIST("histTofNsigmaData_cent"), track.pt(), track.tofNSigmaDe(), event.centFT0C());
          triton_reg.fill(HIST("histTofNsigmaData_cent"), track.pt(), track.tofNSigmaTr(), event.centFT0C());
          Helium3_reg.fill(HIST("histTofNsigmaData_cent"), track.pt() * 2.0, track.tofNSigmaHe(), event.centFT0C());
          Helium4_reg.fill(HIST("histTofNsigmaData_cent"), track.pt() * 2.0, track.tofNSigmaAl(), event.centFT0C());

          proton_reg.fill(HIST("histTofm2_cent"), track.pt(), track.mass() * track.mass(), event.centFT0C());
          deuteron_reg.fill(HIST("histTofm2_cent"), track.pt(), track.mass() * track.mass(), event.centFT0C());
          triton_reg.fill(HIST("histTofm2_cent"), track.pt(), track.mass() * track.mass(), event.centFT0C());
          Helium3_reg.fill(HIST("histTofm2_cent"), track.pt() * 2.0, track.mass() * track.mass(), event.centFT0C());
          Helium4_reg.fill(HIST("histTofm2_cent"), track.pt() * 2.0, track.mass() * track.mass(), event.centFT0C());

          if ((event.centFT0C() > minCentrality) && (event.centFT0C() < maxCentrality)) {
            proton_reg.fill(HIST("histTofm2_eta"), track.pt(), track.mass() * track.mass(), track.eta());
            proton_reg.fill(HIST("histTofNsigmaData_eta"), track.pt(), track.tofNSigmaPr(), track.eta());
            deuteron_reg.fill(HIST("histTofm2_eta"), track.pt(), track.mass() * track.mass(), track.eta());
            deuteron_reg.fill(HIST("histTofNsigmaData_eta"), track.pt(), track.tofNSigmaDe(), track.eta());
            triton_reg.fill(HIST("histTofm2_eta"), track.pt(), track.mass() * track.mass(), track.eta());
            triton_reg.fill(HIST("histTofNsigmaData_eta"), track.pt(), track.tofNSigmaTr(), track.eta());
            Helium3_reg.fill(HIST("histTofm2_eta"), track.pt() * 2.0, track.mass() * track.mass(), track.eta());
            Helium3_reg.fill(HIST("histTofNsigmaData_eta"), track.pt() * 2.0, track.tofNSigmaHe(), track.eta());
            Helium4_reg.fill(HIST("histTofm2_eta"), track.pt() * 2.0, track.mass() * track.mass(), track.eta());
            Helium4_reg.fill(HIST("histTofNsigmaData_eta"), track.pt() * 2.0, track.tofNSigmaAl(), track.eta());
          }
        }
      }

      if (track.sign() < 0) {

        aproton_reg.fill(HIST("histTpcNsigmaData_cent"), track.pt(), track.tpcNSigmaPr(), event.centFT0C());
        adeuteron_reg.fill(HIST("histTpcNsigmaData_cent"), track.pt(), track.tpcNSigmaDe(), event.centFT0C());
        atriton_reg.fill(HIST("histTpcNsigmaData_cent"), track.pt(), track.tpcNSigmaTr(), event.centFT0C());
        aHelium3_reg.fill(HIST("histTpcNsigmaData_cent"), track.pt() * 2.0, track.tpcNSigmaHe(), event.centFT0C());
        aHelium4_reg.fill(HIST("histTpcNsigmaData_cent"), track.pt() * 2.0, track.tpcNSigmaAl(), event.centFT0C());

        if ((event.centFT0C() > minCentrality) && (event.centFT0C() < maxCentrality)) {
          aproton_reg.fill(HIST("histTpcNsigmaData_eta"), track.pt(), track.tpcNSigmaPr(), track.eta());
          adeuteron_reg.fill(HIST("histTpcNsigmaData_eta"), track.pt(), track.tpcNSigmaDe(), track.eta());
          atriton_reg.fill(HIST("histTpcNsigmaData_eta"), track.pt(), track.tpcNSigmaTr(), track.eta());
          aHelium3_reg.fill(HIST("histTpcNsigmaData_eta"), track.pt() * 2.0, track.tpcNSigmaHe(), track.eta());
          aHelium4_reg.fill(HIST("histTpcNsigmaData_eta"), track.pt() * 2.0, track.tpcNSigmaAl(), track.eta());
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
            if (lastLayer < lastRequiredTrdCluster) {
              continue;
            }
          }

          aproton_reg.fill(HIST("histTofNsigmaData_cent"), track.pt(), track.tofNSigmaPr(), event.centFT0C());
          adeuteron_reg.fill(HIST("histTofNsigmaData_cent"), track.pt(), track.tofNSigmaDe(), event.centFT0C());
          atriton_reg.fill(HIST("histTofNsigmaData_cent"), track.pt(), track.tofNSigmaTr(), event.centFT0C());
          aHelium3_reg.fill(HIST("histTofNsigmaData_cent"), track.pt() * 2.0, track.tofNSigmaHe(), event.centFT0C());
          aHelium4_reg.fill(HIST("histTofNsigmaData_cent"), track.pt() * 2.0, track.tofNSigmaAl(), event.centFT0C());

          aproton_reg.fill(HIST("histTofm2_cent"), track.pt(), track.mass() * track.mass(), event.centFT0C());
          adeuteron_reg.fill(HIST("histTofm2_cent"), track.pt(), track.mass() * track.mass(), event.centFT0C());
          atriton_reg.fill(HIST("histTofm2_cent"), track.pt(), track.mass() * track.mass(), event.centFT0C());
          aHelium3_reg.fill(HIST("histTofm2_cent"), track.pt() * 2.0, track.mass() * track.mass(), event.centFT0C());
          aHelium4_reg.fill(HIST("histTofm2_cent"), track.pt() * 2.0, track.mass() * track.mass(), event.centFT0C());

          if ((event.centFT0C() > minCentrality) && (event.centFT0C() < maxCentrality)) {
            aproton_reg.fill(HIST("histTofm2_eta"), track.pt(), track.mass() * track.mass(), track.eta());
            aproton_reg.fill(HIST("histTofNsigmaData_eta"), track.pt(), track.tofNSigmaPr(), track.eta());
            adeuteron_reg.fill(HIST("histTofm2_eta"), track.pt(), track.mass() * track.mass(), track.eta());
            adeuteron_reg.fill(HIST("histTofNsigmaData_eta"), track.pt(), track.tofNSigmaDe(), track.eta());
            atriton_reg.fill(HIST("histTofm2_eta"), track.pt(), track.mass() * track.mass(), track.eta());
            atriton_reg.fill(HIST("histTofNsigmaData_eta"), track.pt(), track.tofNSigmaTr(), track.eta());
            aHelium3_reg.fill(HIST("histTofm2_eta"), track.pt() * 2.0, track.mass() * track.mass(), track.eta());
            aHelium3_reg.fill(HIST("histTofNsigmaData_eta"), track.pt() * 2.0, track.tofNSigmaHe(), track.eta());
            aHelium4_reg.fill(HIST("histTofm2_eta"), track.pt() * 2.0, track.mass() * track.mass(), track.eta());
            aHelium4_reg.fill(HIST("histTofNsigmaData_eta"), track.pt() * 2.0, track.tofNSigmaAl(), track.eta());
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

  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCLfFullPr, aod::pidTOFFullPr, aod::pidTPCLfFullDe, aod::pidTOFFullDe, aod::pidTPCLfFullTr, aod::pidTOFFullTr, aod::pidTPCLfFullHe, aod::pidTOFFullHe, aod::pidTPCLfFullAl, aod::pidTOFFullAl, aod::TrackSelection, aod::TrackSelectionExtension, aod::TOFSignal, aod::pidTOFmass, aod::pidTOFbeta>>;

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

  void processMC(soa::Join<aod::Collisions, aod::McCollisionLabels>::iterator const& collisions, soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCLfFullPr, aod::pidTOFFullPr, aod::pidTPCLfFullDe, aod::pidTOFFullDe, aod::pidTPCLfFullTr, aod::pidTOFFullTr, aod::pidTPCLfFullHe, aod::pidTOFFullHe, aod::pidTPCLfFullAl, aod::pidTOFFullAl, aod::McTrackLabels, aod::TrackSelection, aod::TrackSelectionExtension, aod::TOFSignal, aod::pidTOFmass, aod::pidTOFbeta>> const& tracks,
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
      MC_recon_reg.fill(HIST("histTpcSignalData"), track.pt() * track.sign(), track.tpcSignal(), pdgbin);
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

      if (track.sign() > 0) {
        MC_recon_reg.fill(HIST("histTpcNsigmaDataPr"), track.pt(), track.tpcNSigmaPr());
        MC_recon_reg.fill(HIST("histTpcNsigmaDataDe"), track.pt(), track.tpcNSigmaDe());
        MC_recon_reg.fill(HIST("histTpcNsigmaDataTr"), track.pt(), track.tpcNSigmaTr());
        MC_recon_reg.fill(HIST("histTpcNsigmaDataHe"), track.pt(), track.tpcNSigmaHe());
        MC_recon_reg.fill(HIST("histTpcNsigmaDataAl"), track.pt(), track.tpcNSigmaAl());
      }
      if (track.sign() < 0) {
        MC_recon_reg.fill(HIST("histTpcNsigmaDataaPr"), track.pt(), track.tpcNSigmaPr());
        MC_recon_reg.fill(HIST("histTpcNsigmaDataaDe"), track.pt(), track.tpcNSigmaDe());
        MC_recon_reg.fill(HIST("histTpcNsigmaDataaTr"), track.pt(), track.tpcNSigmaTr());
        MC_recon_reg.fill(HIST("histTpcNsigmaDataaHe"), track.pt(), track.tpcNSigmaHe());
        MC_recon_reg.fill(HIST("histTpcNsigmaDataaAl"), track.pt(), track.tpcNSigmaAl());
      }

      if (track.hasTOF()) {
        Float_t TOFmass2 = ((track.mass()) * (track.mass()));

        MC_recon_reg.fill(HIST("histTOFm2"), track.pt(), TOFmass2, pdgbin);
        MC_recon_reg.fill(HIST("histTofSignalData"), track.pt() * track.sign(), track.beta(), pdgbin);

        if (track.sign() > 0) {
          MC_recon_reg.fill(HIST("histTofNsigmaDataPr"), track.pt(), track.tofNSigmaPr());
          MC_recon_reg.fill(HIST("histTofNsigmaDataDe"), track.pt(), track.tofNSigmaDe());
          MC_recon_reg.fill(HIST("histTofNsigmaDataTr"), track.pt(), track.tofNSigmaTr());
          MC_recon_reg.fill(HIST("histTofNsigmaDataHe"), track.pt(), track.tofNSigmaHe());
          MC_recon_reg.fill(HIST("histTofNsigmaDataAl"), track.pt(), track.tofNSigmaAl());
        }
        if (track.sign() < 0) {
          MC_recon_reg.fill(HIST("histTofNsigmaDataaPr"), track.pt(), track.tofNSigmaPr());
          MC_recon_reg.fill(HIST("histTofNsigmaDataaDe"), track.pt(), track.tofNSigmaDe());
          MC_recon_reg.fill(HIST("histTofNsigmaDataaTr"), track.pt(), track.tofNSigmaTr());
          MC_recon_reg.fill(HIST("histTofNsigmaDataaHe"), track.pt(), track.tofNSigmaHe());
          MC_recon_reg.fill(HIST("histTofNsigmaDataaAl"), track.pt(), track.tofNSigmaAl());
        }
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
  PROCESS_SWITCH(NucleiHistTask, processMC, "process MC", false);
};

//****************************************************************************************************

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<NucleiHistTask>(cfgc, TaskName{"nuclei-hist"})};
}
