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

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace NucleiTableHist
{
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Sign, sign, int8_t);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(IsPVContributor, isPVcontributor, bool);
DECLARE_SOA_COLUMN(HasITS, hasITS, bool);
DECLARE_SOA_COLUMN(HasTPC, hasTPC, bool);
DECLARE_SOA_COLUMN(HasTRD, hasTRD, bool);
DECLARE_SOA_COLUMN(HasTOF, hasTOF, bool);
DECLARE_SOA_COLUMN(ITSnCls, itsNCls, uint8_t);
DECLARE_SOA_COLUMN(TPCNClsFound, tpcNClsFound, int16_t);
DECLARE_SOA_COLUMN(TPCNClsCrossedRows, tpcNClsCrossedRows, int16_t);
DECLARE_SOA_COLUMN(TPCCrossedRowsOverFindableCls, tpcCrossedRowsOverFindableCls, float);
DECLARE_SOA_COLUMN(PassedITSRefit, passedITSRefit, bool);
DECLARE_SOA_COLUMN(PassedTPCRefit, passedTPCRefit, bool);
DECLARE_SOA_COLUMN(ITSChi2NCl, itsChi2NCl, float);
DECLARE_SOA_COLUMN(TPCChi2NCl, tpcChi2NCl, float);
DECLARE_SOA_COLUMN(TPCNSigmaPr, tpcNSigmaPr, float);
DECLARE_SOA_COLUMN(TOFNSigmaPr, tofNSigmaPr, float);
DECLARE_SOA_COLUMN(TPCNSigmaDe, tpcNSigmaDe, float);
DECLARE_SOA_COLUMN(TOFNSigmaDe, tofNSigmaDe, float);
DECLARE_SOA_COLUMN(TPCNSigmaHe, tpcNSigmaHe, float);
DECLARE_SOA_COLUMN(TOFNSigmaHe, tofNSigmaHe, float);
DECLARE_SOA_COLUMN(TOFmass2, tofmass2, float);
DECLARE_SOA_COLUMN(CentFT0C, centFT0C, float);
DECLARE_SOA_COLUMN(DCAxy, dcaxy, int8_t);
DECLARE_SOA_COLUMN(DCAz, dcaz, int8_t);
} // namespace NucleiTableHist
DECLARE_SOA_TABLE(NucleiTable, "AOD", "NUCLEITABLE",
                  NucleiTableHist::Pt,
                  NucleiTableHist::Sign,
                  NucleiTableHist::Eta,
                  NucleiTableHist::IsPVContributor,
                  NucleiTableHist::HasITS,
                  NucleiTableHist::HasTPC,
                  NucleiTableHist::HasTRD,
                  NucleiTableHist::HasTOF,
                  NucleiTableHist::ITSnCls,
                  NucleiTableHist::TPCNClsFound,
                  NucleiTableHist::TPCNClsCrossedRows,
                  NucleiTableHist::TPCCrossedRowsOverFindableCls,
                  NucleiTableHist::PassedITSRefit,
                  NucleiTableHist::PassedTPCRefit,
                  NucleiTableHist::ITSChi2NCl,
                  NucleiTableHist::TPCChi2NCl,
                  NucleiTableHist::TPCNSigmaPr,
                  NucleiTableHist::TOFNSigmaPr,
                  NucleiTableHist::TPCNSigmaDe,
                  NucleiTableHist::TOFNSigmaDe,
                  NucleiTableHist::TPCNSigmaHe,
                  NucleiTableHist::TOFNSigmaHe,
                  NucleiTableHist::TOFmass2,
                  NucleiTableHist::CentFT0C,
                  NucleiTableHist::DCAxy,
                  NucleiTableHist::DCAz)
} // namespace o2::aod

struct NucleiHistTask {

  Produces<o2::aod::NucleiTable> nucleiTable;

  HistogramRegistry spectra{"spectra", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry proton_erg{"proton", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry aproton_erg{"aproton", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry deuteron_reg{"deuteron", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry adeuteron_reg{"adeuteron", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry Helium3_reg{"Helium3", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry aHelium3_reg{"aHelium3", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  void init(o2::framework::InitContext&)
  {

    if (doprocessData == true && doprocessDataCent == true) {
      LOG(fatal) << "Can't enable processData and processDataCent in the same time, pick one!";
    }

    std::vector<double> ptBinning = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 3.6, 4., 5., 6., 8., 10., 12., 14.};
    std::vector<double> centBinning = {0., 1., 5., 10., 20., 30., 40., 50., 70., 100.};

    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec centAxis = {centBinning, "V0M (%)"};
    AxisSpec centralityAxis = {100, 0.0, 100.0, "VT0C (%)"};

    // QA histograms
    spectra.add("histRecVtxZData", "collision z position", HistType::kTH1F, {{200, -20., +20., "z position (cm)"}});
    spectra.add("histTpcSignalData", "Specific energy loss", HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {1400, 0, 1400, "d#it{E} / d#it{X} (a. u.)"}});
    spectra.add("histTofSignalData", "TOF signal", HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    spectra.add("histDcaVsPtData_particle", "dcaXY vs Pt (particle)", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    spectra.add("histDcaZVsPtData_particle", "dcaZ vs Pt (particle)", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    spectra.add("histDcaVsPtData_antiparticle", "dcaXY vs Pt (antiparticle)", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    spectra.add("histDcaZVsPtData_antiparticle", "dcaZ vs Pt (antiparticle)", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    spectra.add("histTOFm2", "TOF m^2 vs Pt", HistType::kTH2F, {ptAxis, {400, 0.0, 10.0, "m^2"}});
    spectra.add("histNClusterTPC", "Number of Clusters in TPC vs Pt", HistType::kTH2F, {ptAxis, {160, 0.0, 160.0, "nCluster"}});
    spectra.add("histNClusterITS", "Number of Clusters in ITS vs Pt", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    spectra.add("histChi2TPC", "chi^2 TPC vs Pt", HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    spectra.add("histChi2ITS", "chi^2 ITS vs Pt", HistType::kTH2F, {ptAxis, {500, 0.0, 50.0, "chi^2"}});

    // histograms for Proton
    proton_erg.add("histKeepEventData", "skimming histogram (p)", HistType::kTH1F, {{2, -0.5, +1.5, "true: keep event, false: reject event"}});
    proton_erg.add("histTpcSignalData", "Specific energy loss (p)", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {1400, 0, 1400, "d#it{E} / d#it{X} (a. u.)"}});
    proton_erg.add("histTofSignalData", "TOF signal (p)", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    proton_erg.add("histDcaVsPtData", "dcaXY vs Pt (p)", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    proton_erg.add("histDcaZVsPtData", "dcaZ vs Pt (p)", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    proton_erg.add("histTOFm2", "TOF m^2 vs Pt (p)", HistType::kTH2F, {ptAxis, {400, 0.0, 10.0, "m^2"}});
    proton_erg.add("histTpcNsigmaData", "n-sigma TPC (p)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{p}"}});
    proton_erg.add("histTofNsigmaData", "n-sigma TOF (p)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{p}"}});
    proton_erg.add("histNClusterTPC", "Number of Clusters in TPC vs Pt (p)", HistType::kTH2F, {ptAxis, {160, 0.0, 160.0, "nCluster"}});
    proton_erg.add("histNClusterITS", "Number of Clusters in ITS vs Pt (p)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    proton_erg.add("histChi2TPC", "chi^2 TPC vs Pt (p)", HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    proton_erg.add("histChi2ITS", "chi^2 ITS vs Pt (p)", HistType::kTH2F, {ptAxis, {500, 0.0, 50.0, "chi^2"}});
    proton_erg.add("histTpcNsigmaData_cent", "n-sigma TPC (p) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{p}"}, centralityAxis});
    proton_erg.add("histTofNsigmaData_cent", "n-sigma TOF (p) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{p}"}, centralityAxis});
    proton_erg.add("histTofm2_cent", "mass^2 TOF (p) centrality", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{p}"}, centralityAxis});

    // histograms for antiProton
    aproton_erg.add("histKeepEventData", "skimming histogram (antip)", HistType::kTH1F, {{2, -0.5, +1.5, "true: keep event, false: reject event"}});
    aproton_erg.add("histTpcSignalData", "Specific energy loss (antip)", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {1400, 0, 1400, "d#it{E} / d#it{X} (a. u.)"}});
    aproton_erg.add("histTofSignalData", "TOF signal (antip)", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    aproton_erg.add("histDcaVsPtData", "dcaXY vs Pt (antip)", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    aproton_erg.add("histDcaZVsPtData", "dcaZ vs Pt (antip)", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    aproton_erg.add("histTOFm2", "TOF m^2 vs Pt (antip)", HistType::kTH2F, {ptAxis, {400, 0.0, 10.0, "m^2"}});
    aproton_erg.add("histTpcNsigmaData", "n-sigma TPC (antip)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{antip}"}});
    aproton_erg.add("histTofNsigmaData", "n-sigma TOF (antip)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{antip}"}});
    aproton_erg.add("histNClusterTPC", "Number of Clusters in TPC vs Pt (antip)", HistType::kTH2F, {ptAxis, {160, 0.0, 160.0, "nCluster"}});
    aproton_erg.add("histNClusterITS", "Number of Clusters in ITS vs Pt (antip)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    aproton_erg.add("histChi2TPC", "chi^2 TPC vs Pt (antip)", HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    aproton_erg.add("histChi2ITS", "chi^2 ITS vs Pt (antip)", HistType::kTH2F, {ptAxis, {500, 0.0, 50.0, "chi^2"}});
    aproton_erg.add("histTpcNsigmaData_cent", "n-sigma TPC (antip) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{antip}"}, centralityAxis});
    aproton_erg.add("histTofNsigmaData_cent", "n-sigma TOF (antip) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{antip}"}, centralityAxis});
    aproton_erg.add("histTofm2_cent", "mass^2 TOF (antip) centrality", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{antip}"}, centralityAxis});

    // histograms for Deuterons
    deuteron_reg.add("histKeepEventData", "skimming histogram (d)", HistType::kTH1F, {{2, -0.5, +1.5, "true: keep event, false: reject event"}});
    deuteron_reg.add("histTpcSignalData", "Specific energy loss (d)", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {1400, 0, 1400, "d#it{E} / d#it{X} (a. u.)"}});
    deuteron_reg.add("histTofSignalData", "TOF signal (d)", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    deuteron_reg.add("histDcaVsPtData", "dcaXY vs Pt (d)", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    deuteron_reg.add("histDcaZVsPtData", "dcaZ vs Pt (d)", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    deuteron_reg.add("histTOFm2", "TOF m^2 vs Pt (d)", HistType::kTH2F, {ptAxis, {400, 0.0, 10.0, "m^2"}});
    deuteron_reg.add("histTpcNsigmaData", "n-sigma TPC (d)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{d}"}});
    deuteron_reg.add("histTofNsigmaData", "n-sigma TOF (d)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{d}"}});
    deuteron_reg.add("histNClusterTPC", "Number of Clusters in TPC vs Pt (d)", HistType::kTH2F, {ptAxis, {160, 0.0, 160.0, "nCluster"}});
    deuteron_reg.add("histNClusterITS", "Number of Clusters in ITS vs Pt (d)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    deuteron_reg.add("histChi2TPC", "chi^2 TPC vs Pt (d)", HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    deuteron_reg.add("histChi2ITS", "chi^2 ITS vs Pt (d)", HistType::kTH2F, {ptAxis, {500, 0.0, 50.0, "chi^2"}});
    deuteron_reg.add("histTpcNsigmaData_cent", "n-sigma TPC (d) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{d}"}, centralityAxis});
    deuteron_reg.add("histTofNsigmaData_cent", "n-sigma TOF (d) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{d}"}, centralityAxis});
    deuteron_reg.add("histTofm2_cent", "mass^2 TOF (d) centrality", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{d}"}, centralityAxis});

    // histograms for antiDeuterons
    adeuteron_reg.add("histKeepEventData", "skimming histogram (antid)", HistType::kTH1F, {{2, -0.5, +1.5, "true: keep event, false: reject event"}});
    adeuteron_reg.add("histTpcSignalData", "Specific energy loss (antid)", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {1400, 0, 1400, "d#it{E} / d#it{X} (a. u.)"}});
    adeuteron_reg.add("histTofSignalData", "TOF signal (antid)", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    adeuteron_reg.add("histDcaVsPtData", "dcaXY vs Pt (antid)", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    adeuteron_reg.add("histDcaZVsPtData", "dcaZ vs Pt (antid)", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    adeuteron_reg.add("histTOFm2", "TOF m^2 vs Pt (antid)", HistType::kTH2F, {ptAxis, {400, 0.0, 10.0, "m^2"}});
    adeuteron_reg.add("histTpcNsigmaData", "n-sigma TPC (antid)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{antid}"}});
    adeuteron_reg.add("histTofNsigmaData", "n-sigma TOF (antid)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{antid}"}});
    adeuteron_reg.add("histNClusterTPC", "Number of Clusters in TPC vs Pt (antid)", HistType::kTH2F, {ptAxis, {160, 0.0, 160.0, "nCluster"}});
    adeuteron_reg.add("histNClusterITS", "Number of Clusters in ITS vs Pt (antid)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    adeuteron_reg.add("histChi2TPC", "chi^2 TPC vs Pt (antid)", HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    adeuteron_reg.add("histChi2ITS", "chi^2 ITS vs Pt (antid)", HistType::kTH2F, {ptAxis, {500, 0.0, 50.0, "chi^2"}});
    adeuteron_reg.add("histTpcNsigmaData_cent", "n-sigma TPC (antid) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{antid}"}, centralityAxis});
    adeuteron_reg.add("histTofNsigmaData_cent", "n-sigma TOF (antid) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{antid}"}, centralityAxis});
    adeuteron_reg.add("histTofm2_cent", "mass^2 TOF (antid) centrality", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{antid}"}, centralityAxis});

    // histograms for Helium-3
    Helium3_reg.add("histKeepEventData", "skimming histogram (He-3)", HistType::kTH1F, {{2, -0.5, +1.5, "true: keep event, false: reject event"}});
    Helium3_reg.add("histTpcSignalData", "Specific energy loss (He-3)", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {1400, 0, 1400, "d#it{E} / d#it{X} (a. u.)"}});
    Helium3_reg.add("histTofSignalData", "TOF signal (He-3)", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    Helium3_reg.add("histDcaVsPtData", "dcaXY vs Pt (He-3)", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    Helium3_reg.add("histDcaZVsPtData", "dcaZ vs Pt (He-3)", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    Helium3_reg.add("histTOFm2", "TOF m^2 vs Pt (He-3)", HistType::kTH2F, {ptAxis, {400, 0.0, 10.0, "m^2"}});
    Helium3_reg.add("histTpcNsigmaData", "n-sigma TPC (He-3)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{He-3}"}});
    Helium3_reg.add("histTofNsigmaData", "n-sigma TOF (He-3)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{He-3}"}});
    Helium3_reg.add("histNClusterTPC", "Number of Clusters in TPC vs Pt (He-3)", HistType::kTH2F, {ptAxis, {160, 0.0, 160.0, "nCluster"}});
    Helium3_reg.add("histNClusterITS", "Number of Clusters in ITS vs Pt (He-3)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    Helium3_reg.add("histChi2TPC", "chi^2 TPC vs Pt (He-3)", HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    Helium3_reg.add("histChi2ITS", "chi^2 ITS vs Pt (He-3)", HistType::kTH2F, {ptAxis, {500, 0.0, 50.0, "chi^2"}});
    Helium3_reg.add("histTpcNsigmaData_cent", "n-sigma TPC (He-3) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{He-3}"}, centralityAxis});
    Helium3_reg.add("histTofNsigmaData_cent", "n-sigma TOF (He-3) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{He-3}"}, centralityAxis});
    Helium3_reg.add("histTofm2_cent", "mass^2 TOF (He-3) centrality", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{He-3}"}, centralityAxis});

    // histograms for antiHelium-3
    aHelium3_reg.add("histKeepEventData", "skimming histogram (antiHe-3)", HistType::kTH1F, {{2, -0.5, +1.5, "true: keep event, false: reject event"}});
    aHelium3_reg.add("histTpcSignalData", "Specific energy loss (antiHe-3)", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {1400, 0, 1400, "d#it{E} / d#it{X} (a. u.)"}});
    aHelium3_reg.add("histTofSignalData", "TOF signal (antiHe-3)", HistType::kTH2F, {{600, 0., 6., "#it{p} (GeV/#it{c})"}, {550, 0.0, 1.1, "#beta (TOF)"}});
    aHelium3_reg.add("histDcaVsPtData", "dcaXY vs Pt (antiHe-3)", HistType::kTH2F, {ptAxis, {250, -0.5, 0.5, "dca"}});
    aHelium3_reg.add("histDcaZVsPtData", "dcaZ vs Pt (antiHe-3)", HistType::kTH2F, {ptAxis, {1000, -2.0, 2.0, "dca"}});
    aHelium3_reg.add("histTOFm2", "TOF m^2 vs Pt (antiHe-3)", HistType::kTH2F, {ptAxis, {400, 0.0, 10.0, "m^2"}});
    aHelium3_reg.add("histTpcNsigmaData", "n-sigma TPC (antiHe-3)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{antiHe-3}"}});
    aHelium3_reg.add("histTofNsigmaData", "n-sigma TOF (antiHe-3)", HistType::kTH2F, {ptAxis, {160, -20., +20., "n#sigma_{antiHe-3}"}});
    aHelium3_reg.add("histNClusterTPC", "Number of Clusters in TPC vs Pt (antiHe-3)", HistType::kTH2F, {ptAxis, {160, 0.0, 160.0, "nCluster"}});
    aHelium3_reg.add("histNClusterITS", "Number of Clusters in ITS vs Pt (antiHe-3)", HistType::kTH2F, {ptAxis, {10, 0.0, 10.0, "nCluster"}});
    aHelium3_reg.add("histChi2TPC", "chi^2 TPC vs Pt (antiHe-3)", HistType::kTH2F, {ptAxis, {100, 0.0, 5.0, "chi^2"}});
    aHelium3_reg.add("histChi2ITS", "chi^2 ITS vs Pt (antiHe-3)", HistType::kTH2F, {ptAxis, {500, 0.0, 50.0, "chi^2"}});
    aHelium3_reg.add("histTpcNsigmaData_cent", "n-sigma TPC (antiHe-3) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{antiHe-3}"}, centralityAxis});
    aHelium3_reg.add("histTofNsigmaData_cent", "n-sigma TOF (antiHe-3) centrality", HistType::kTH3F, {ptAxis, {160, -20., +20., "n#sigma_{antiHe-3}"}, centralityAxis});
    aHelium3_reg.add("histTofm2_cent", "mass^2 TOF (antiHe-3) centrality", HistType::kTH3F, {ptAxis, {400, 0.0, 10.0, "m^2_{antiHe-3}"}, centralityAxis});
  }

  Configurable<float> yMin{"yMin", -0.5, "Maximum rapidity"};
  Configurable<float> yMax{"yMax", 0.5, "Minimum rapidity"};

  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<float> nsigmacutLow{"nsigmacutLow", -3.0, "Value of the Nsigma cut"};
  Configurable<float> nsigmacutHigh{"nsigmacutHigh", +3.0, "Value of the Nsigma cut"};

  // Replacement for globalTrack filter
  Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 70.0f, "min number of crossed rows TPC"};
  Configurable<float> minRatioCrossedRowsTPC{"minRatioCrossedRowsTPC", 0.8f, "min ratio of crossed rows over findable clusters TPC"};
  Configurable<float> maxChi2TPC{"maxChi2TPC", 4.0f, "max chi2 per cluster TPC"};
  Configurable<float> maxChi2ITS{"maxChi2ITS", 36.0f, "max chi2 per cluster ITS"};
  Configurable<float> maxDCA_Z{"maxDCA_Z", 2.0f, "max DCA to vertex z"};
  Configurable<float> maxDCA_XY{"maxDCA_XY", 0.5f, "max DCA to vertex xy"};
  Configurable<float> minReqClusterITS{"minReqClusterITS", 1.0, "min number of clusters required in ITS"};
  Configurable<float> pTmin{"pTmin", 0.1f, "min pT"};
  Configurable<float> pTmax{"pTmax", 1e+10f, "max pT"};

  template <typename CollisionType, typename TracksType>
  void fillHistograms(const CollisionType& event, const TracksType& tracks)
  {

    // collision process loop
    bool keepEvent_p = kFALSE;
    bool keepEvent_d = kFALSE;
    bool keepEvent_He3 = kFALSE;

    bool keepEvent_antip = kFALSE;
    bool keepEvent_antid = kFALSE;
    bool keepEvent_antiHe3 = kFALSE;

    spectra.fill(HIST("histRecVtxZData"), event.posZ());

    for (auto track : tracks) { // start loop over tracks

      float TPCnumberCls = track.tpcNClsFound();
      float RatioCrossedRowsOverFindableTPC = track.tpcCrossedRowsOverFindableCls();
      float Chi2perClusterTPC = track.tpcChi2NCl();
      float Chi2perClusterITS = track.itsChi2NCl();

      if (TPCnumberCls < minNCrossedRowsTPC || RatioCrossedRowsOverFindableTPC < minRatioCrossedRowsTPC || Chi2perClusterTPC > maxChi2TPC || Chi2perClusterITS > maxChi2ITS || !(track.passedTPCRefit()) || !(track.passedITSRefit()) || (track.itsNCls()) < minReqClusterITS) {
        continue;
      }

      if (track.sign() > 0) {
        spectra.fill(HIST("histDcaVsPtData_particle"), track.pt(), track.dcaXY());
        spectra.fill(HIST("histDcaZVsPtData_particle"), track.pt(), track.dcaZ());
      }

      if (track.sign() < 0) {
        spectra.fill(HIST("histDcaVsPtData_antiparticle"), track.pt(), track.dcaXY());
        spectra.fill(HIST("histDcaZVsPtData_antiparticle"), track.pt(), track.dcaZ());
      }

      if (TMath::Abs(track.dcaXY()) > maxDCA_XY || TMath::Abs(track.dcaZ()) > maxDCA_Z) {
        continue;
      }

      // cut on rapidity
      TLorentzVector lorentzVector_proton{};
      TLorentzVector lorentzVector_deuteron{};
      TLorentzVector lorentzVector_He3{};

      lorentzVector_proton.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassProton);
      lorentzVector_deuteron.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassDeuteron);
      lorentzVector_He3.SetPtEtaPhiM(track.pt() * 2.0, track.eta(), track.phi(), constants::physics::MassHelium3);

      if (lorentzVector_proton.Rapidity() < yMin || lorentzVector_proton.Rapidity() > yMax ||
          lorentzVector_deuteron.Rapidity() < yMin || lorentzVector_deuteron.Rapidity() > yMax ||
          lorentzVector_He3.Rapidity() < yMin || lorentzVector_He3.Rapidity() > yMax) {
        continue;
      }

      // fill QA histograms
      float nSigmaProton = track.tpcNSigmaPr();
      float nSigmaDeut = track.tpcNSigmaDe();
      float nSigmaHe3 = track.tpcNSigmaHe();

      spectra.fill(HIST("histTpcSignalData"), track.tpcInnerParam() * track.sign(), track.tpcSignal());
      spectra.fill(HIST("histNClusterTPC"), track.pt(), track.tpcNClsCrossedRows());
      spectra.fill(HIST("histNClusterITS"), track.pt(), track.itsNCls());
      spectra.fill(HIST("histChi2TPC"), track.pt(), track.tpcChi2NCl());
      spectra.fill(HIST("histChi2ITS"), track.pt(), track.itsChi2NCl());

      if (track.sign() > 0) {
        proton_erg.fill(HIST("histTpcNsigmaData"), track.pt(), nSigmaProton);
        deuteron_reg.fill(HIST("histTpcNsigmaData"), track.pt(), nSigmaDeut);
        Helium3_reg.fill(HIST("histTpcNsigmaData"), track.pt() * 2.0, nSigmaHe3);

        //  fill TOF m^2 histogram
        if (track.hasTOF()) {

          Float_t TOFmass2 = ((track.mass()) * (track.mass()));

          spectra.fill(HIST("histTOFm2"), track.pt(), TOFmass2);
        }
      }

      if (track.sign() < 0) {
        aproton_erg.fill(HIST("histTpcNsigmaData"), track.pt(), nSigmaProton);
        adeuteron_reg.fill(HIST("histTpcNsigmaData"), track.pt(), nSigmaDeut);
        aHelium3_reg.fill(HIST("histTpcNsigmaData"), track.pt() * 2.0, nSigmaHe3);

        // fill TOF m^2 histogram
        if (track.hasTOF()) {

          Float_t TOFmass2 = ((track.mass()) * (track.mass()));

          spectra.fill(HIST("histTOFm2"), track.pt(), TOFmass2);
        }
      }

      //**************   check offline-trigger (skimming) condidition Proton   *******************

      if (nSigmaProton > nsigmacutLow && nSigmaProton < nsigmacutHigh) {

        if (track.sign() > 0) {
          keepEvent_p = kTRUE;

          proton_erg.fill(HIST("histDcaVsPtData"), track.pt(), track.dcaXY());
          proton_erg.fill(HIST("histDcaZVsPtData"), track.pt(), track.dcaZ());
          proton_erg.fill(HIST("histTpcSignalData"), track.tpcInnerParam(), track.tpcSignal());
          proton_erg.fill(HIST("histNClusterTPC"), track.pt(), track.tpcNClsFound());
          proton_erg.fill(HIST("histNClusterITS"), track.pt(), track.itsNCls());
          proton_erg.fill(HIST("histChi2TPC"), track.pt(), track.tpcChi2NCl());
          proton_erg.fill(HIST("histChi2ITS"), track.pt(), track.itsChi2NCl());

          if (track.hasTOF()) {

            Float_t TOFmass2 = ((track.mass()) * (track.mass()));
            Float_t beta = track.beta();

            proton_erg.fill(HIST("histTOFm2"), track.pt(), TOFmass2);
            proton_erg.fill(HIST("histTofSignalData"), track.tpcInnerParam(), beta);
            proton_erg.fill(HIST("histTofNsigmaData"), track.pt(), track.tofNSigmaPr());
          }
        }

        if (track.sign() < 0) {
          keepEvent_antip = kTRUE;

          aproton_erg.fill(HIST("histDcaVsPtData"), track.pt(), track.dcaXY());
          aproton_erg.fill(HIST("histDcaZVsPtData"), track.pt(), track.dcaZ());
          aproton_erg.fill(HIST("histTpcSignalData"), track.tpcInnerParam(), track.tpcSignal());
          aproton_erg.fill(HIST("histNClusterTPC"), track.pt(), track.tpcNClsFound());
          aproton_erg.fill(HIST("histNClusterITS"), track.pt(), track.itsNCls());
          aproton_erg.fill(HIST("histChi2TPC"), track.pt(), track.tpcChi2NCl());
          aproton_erg.fill(HIST("histChi2ITS"), track.pt(), track.itsChi2NCl());

          if (track.hasTOF()) {

            Float_t TOFmass2 = ((track.mass()) * (track.mass()));
            Float_t beta = track.beta();

            aproton_erg.fill(HIST("histTOFm2"), track.pt(), TOFmass2);
            aproton_erg.fill(HIST("histTofSignalData"), track.tpcInnerParam(), beta);
            aproton_erg.fill(HIST("histTofNsigmaData"), track.pt(), track.tofNSigmaPr());
          }
        }

        if (track.hasTOF()) {
          spectra.fill(HIST("histTofSignalData"), track.tpcInnerParam() * track.sign(), track.beta());
        }

      }

      //**************   check offline-trigger (skimming) condidition Deuteron   *******************

      if (nSigmaDeut > nsigmacutLow && nSigmaDeut < nsigmacutHigh) {

        if (track.sign() > 0) {
          keepEvent_d = kTRUE;

          deuteron_reg.fill(HIST("histDcaVsPtData"), track.pt(), track.dcaXY());
          deuteron_reg.fill(HIST("histDcaZVsPtData"), track.pt(), track.dcaZ());
          deuteron_reg.fill(HIST("histTpcSignalData"), track.tpcInnerParam(), track.tpcSignal());
          deuteron_reg.fill(HIST("histNClusterTPC"), track.pt(), track.tpcNClsFound());
          deuteron_reg.fill(HIST("histNClusterITS"), track.pt(), track.itsNCls());
          deuteron_reg.fill(HIST("histChi2TPC"), track.pt(), track.tpcChi2NCl());
          deuteron_reg.fill(HIST("histChi2ITS"), track.pt(), track.itsChi2NCl());

          if (track.hasTOF()) {

            Float_t TOFmass2 = ((track.mass()) * (track.mass()));
            Float_t beta = track.beta();

            deuteron_reg.fill(HIST("histTOFm2"), track.pt(), TOFmass2);
            deuteron_reg.fill(HIST("histTofSignalData"), track.tpcInnerParam(), beta);
            deuteron_reg.fill(HIST("histTofNsigmaData"), track.pt(), track.tofNSigmaDe());
          }
        }

        if (track.sign() < 0) {
          keepEvent_antid = kTRUE;

          adeuteron_reg.fill(HIST("histDcaVsPtData"), track.pt(), track.dcaXY());
          adeuteron_reg.fill(HIST("histDcaZVsPtData"), track.pt(), track.dcaZ());
          adeuteron_reg.fill(HIST("histTpcSignalData"), track.tpcInnerParam(), track.tpcSignal());
          adeuteron_reg.fill(HIST("histNClusterTPC"), track.pt(), track.tpcNClsFound());
          adeuteron_reg.fill(HIST("histNClusterITS"), track.pt(), track.itsNCls());
          adeuteron_reg.fill(HIST("histChi2TPC"), track.pt(), track.tpcChi2NCl());
          adeuteron_reg.fill(HIST("histChi2ITS"), track.pt(), track.itsChi2NCl());

          if (track.hasTOF()) {

            Float_t TOFmass2 = ((track.mass()) * (track.mass()));
            Float_t beta = track.beta();

            adeuteron_reg.fill(HIST("histTOFm2"), track.pt(), TOFmass2);
            adeuteron_reg.fill(HIST("histTofSignalData"), track.tpcInnerParam(), beta);
            adeuteron_reg.fill(HIST("histTofNsigmaData"), track.pt(), track.tofNSigmaDe());
          }
        }

        if (track.hasTOF()) {
          spectra.fill(HIST("histTofSignalData"), track.tpcInnerParam() * track.sign(), track.beta());
        }

      }

      //**************   check offline-trigger (skimming) condidition Helium-3   *******************

      if (nSigmaHe3 > nsigmacutLow && nSigmaHe3 < nsigmacutHigh) {

        if (track.sign() > 0) {
          keepEvent_He3 = kTRUE;

          Helium3_reg.fill(HIST("histDcaVsPtData"), track.pt() * 2.0, track.dcaXY());
          Helium3_reg.fill(HIST("histDcaZVsPtData"), track.pt() * 2.0, track.dcaZ());
          Helium3_reg.fill(HIST("histTpcSignalData"), track.tpcInnerParam(), track.tpcSignal());
          Helium3_reg.fill(HIST("histNClusterTPC"), track.pt() * 2.0, track.tpcNClsFound());
          Helium3_reg.fill(HIST("histNClusterITS"), track.pt() * 2.0, track.itsNCls());
          Helium3_reg.fill(HIST("histChi2TPC"), track.pt() * 2.0, track.tpcChi2NCl());
          Helium3_reg.fill(HIST("histChi2ITS"), track.pt() * 2.0, track.itsChi2NCl());

          if (track.hasTOF()) {

            Float_t TOFmass2 = ((track.mass()) * (track.mass()));
            Float_t beta = track.beta();

            Helium3_reg.fill(HIST("histTOFm2"), track.pt() * 2.0, TOFmass2);
            Helium3_reg.fill(HIST("histTofSignalData"), track.tpcInnerParam(), beta);
            Helium3_reg.fill(HIST("histTofNsigmaData"), track.pt() * 2.0, track.tofNSigmaHe());
          }
        }

        if (track.sign() < 0) {
          keepEvent_antiHe3 = kTRUE;
          aHelium3_reg.fill(HIST("histDcaVsPtData"), track.pt() * 2.0, track.dcaXY());
          aHelium3_reg.fill(HIST("histDcaZVsPtData"), track.pt() * 2.0, track.dcaZ());
          aHelium3_reg.fill(HIST("histTpcSignalData"), track.tpcInnerParam(), track.tpcSignal());
          aHelium3_reg.fill(HIST("histNClusterTPC"), track.pt() * 2.0, track.tpcNClsFound());
          aHelium3_reg.fill(HIST("histNClusterITS"), track.pt() * 2.0, track.itsNCls());
          aHelium3_reg.fill(HIST("histChi2TPC"), track.pt() * 2.0, track.tpcChi2NCl());
          aHelium3_reg.fill(HIST("histChi2ITS"), track.pt() * 2.0, track.itsChi2NCl());

          if (track.hasTOF()) {

            Float_t TOFmass2 = ((track.mass()) * (track.mass()));
            Float_t beta = track.beta();

            aHelium3_reg.fill(HIST("histTOFm2"), track.pt() * 2.0, TOFmass2);
            aHelium3_reg.fill(HIST("histTofSignalData"), track.tpcInnerParam(), beta);
            aHelium3_reg.fill(HIST("histTofNsigmaData"), track.pt() * 2.0, track.tofNSigmaHe());
          }
        }

        if (track.hasTOF()) {
          spectra.fill(HIST("histTofSignalData"), track.tpcInnerParam() * 2.0 * track.sign(), track.beta());
        }

      }

    } // end loop over tracks

    // fill trigger (skimming) results
    proton_erg.fill(HIST("histKeepEventData"), keepEvent_p);
    aproton_erg.fill(HIST("histKeepEventData"), keepEvent_antip);
    deuteron_reg.fill(HIST("histKeepEventData"), keepEvent_d);
    adeuteron_reg.fill(HIST("histKeepEventData"), keepEvent_antid);
    Helium3_reg.fill(HIST("histKeepEventData"), keepEvent_He3);
    aHelium3_reg.fill(HIST("histKeepEventData"), keepEvent_antiHe3);
  }

  //****************************************************************************************************

  template <typename CollisionType, typename TracksType>
  void fillCentHistorgrams(const CollisionType& event, const TracksType& tracks)
  {

    for (auto track : tracks) { // start loop over tracks

      float TPCnumberCls = track.tpcNClsFound();
      float RatioCrossedRowsOverFindableTPC = track.tpcCrossedRowsOverFindableCls();
      float Chi2perClusterTPC = track.tpcChi2NCl();
      float Chi2perClusterITS = track.itsChi2NCl();

      // track cuts
      if (TPCnumberCls < minNCrossedRowsTPC || RatioCrossedRowsOverFindableTPC < minRatioCrossedRowsTPC || Chi2perClusterTPC > maxChi2TPC || Chi2perClusterITS > maxChi2ITS || !(track.passedTPCRefit()) || !(track.passedITSRefit()) || (track.itsNCls()) < minReqClusterITS || (track.pt()) < pTmin || (track.pt()) > pTmax || TMath::Abs(track.dcaXY()) > maxDCA_XY || TMath::Abs(track.dcaZ()) > maxDCA_Z) {
        continue;
      }

      // cut on rapidity
      TLorentzVector lorentzVector_proton{};
      TLorentzVector lorentzVector_deuteron{};
      TLorentzVector lorentzVector_He3{};

      lorentzVector_proton.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassProton);
      lorentzVector_deuteron.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassDeuteron);
      lorentzVector_He3.SetPtEtaPhiM(track.pt() * 2.0, track.eta(), track.phi(), constants::physics::MassHelium3);

      if (lorentzVector_proton.Rapidity() < yMin || lorentzVector_proton.Rapidity() > yMax ||
          lorentzVector_deuteron.Rapidity() < yMin || lorentzVector_deuteron.Rapidity() > yMax ||
          lorentzVector_He3.Rapidity() < yMin || lorentzVector_He3.Rapidity() > yMax) {
        continue;
      }

      // fill 3D centrality histograms
      if (track.sign() > 0) {

        proton_erg.fill(HIST("histTpcNsigmaData_cent"), track.pt(), track.tpcNSigmaPr(), event.centFT0C());
        proton_erg.fill(HIST("histTofNsigmaData_cent"), track.pt(), track.tofNSigmaPr(), event.centFT0C());
        deuteron_reg.fill(HIST("histTpcNsigmaData_cent"), track.pt(), track.tpcNSigmaDe(), event.centFT0C());
        deuteron_reg.fill(HIST("histTofNsigmaData_cent"), track.pt(), track.tofNSigmaDe(), event.centFT0C());
        Helium3_reg.fill(HIST("histTpcNsigmaData_cent"), track.pt() * 2.0, track.tpcNSigmaHe(), event.centFT0C());
        Helium3_reg.fill(HIST("histTofNsigmaData_cent"), track.pt() * 2.0, track.tofNSigmaHe(), event.centFT0C());

        if (track.hasTOF()) {
          proton_erg.fill(HIST("histTofm2_cent"), track.pt(), track.mass() * track.mass(), event.centFT0C());
          deuteron_reg.fill(HIST("histTofm2_cent"), track.pt(), track.mass() * track.mass(), event.centFT0C());
          Helium3_reg.fill(HIST("histTofm2_cent"), track.pt() * 2.0, track.mass() * track.mass(), event.centFT0C());
        }
      }

      if (track.sign() < 0) {

        aproton_erg.fill(HIST("histTpcNsigmaData_cent"), track.pt(), track.tpcNSigmaPr(), event.centFT0C());
        aproton_erg.fill(HIST("histTofNsigmaData_cent"), track.pt(), track.tofNSigmaPr(), event.centFT0C());
        adeuteron_reg.fill(HIST("histTpcNsigmaData_cent"), track.pt(), track.tpcNSigmaDe(), event.centFT0C());
        adeuteron_reg.fill(HIST("histTofNsigmaData_cent"), track.pt(), track.tofNSigmaDe(), event.centFT0C());
        aHelium3_reg.fill(HIST("histTpcNsigmaData_cent"), track.pt() * 2.0, track.tpcNSigmaHe(), event.centFT0C());
        aHelium3_reg.fill(HIST("histTofNsigmaData_cent"), track.pt() * 2.0, track.tofNSigmaHe(), event.centFT0C());

        if (track.hasTOF()) {
          aproton_erg.fill(HIST("histTofm2_cent"), track.pt(), track.mass() * track.mass(), event.centFT0C());
          adeuteron_reg.fill(HIST("histTofm2_cent"), track.pt(), track.mass() * track.mass(), event.centFT0C());
          aHelium3_reg.fill(HIST("histTofm2_cent"), track.pt() * 2.0, track.mass() * track.mass(), event.centFT0C());
        }
      }
    }
  }

  //****************************************************************************************************

  template <typename CollisionType, typename TracksType>
  void fillTable(const CollisionType& event, const TracksType& tracks)
  {

    for (auto track : tracks) {

      if ((track.pt()) < pTmin || (track.pt()) > pTmax) {
        continue;
      }

      // fill table
      nucleiTable(
        track.pt(),
        track.sign(),
        track.eta(),
        track.isPVContributor(),
        track.hasITS(),
        track.hasTPC(),
        track.hasTRD(),
        track.hasTOF(),
        track.itsNCls(),
        track.tpcNClsFound(),
        track.tpcNClsCrossedRows(),
        track.tpcCrossedRowsOverFindableCls(),
        track.passedITSRefit(),
        track.passedTPCRefit(),
        track.itsChi2NCl(),
        track.tpcChi2NCl(),
        track.tpcNSigmaPr(),
        track.tofNSigmaPr(),
        track.tpcNSigmaDe(),
        track.tofNSigmaDe(),
        track.tpcNSigmaHe(),
        track.tofNSigmaHe(),
        track.mass() * track.mass(),
        0.0,
        track.dcaXY(),
        track.dcaZ());
    }
  }

  //****************************************************************************************************

  template <typename CollisionType, typename TracksType>
  void fillTableCent(const CollisionType& event, const TracksType& tracks)
  {

    for (auto track : tracks) {

      if ((track.pt()) < pTmin || (track.pt()) > pTmax) {
        continue;
      }

      // fill table
      nucleiTable(
        track.pt(),
        track.sign(),
        track.eta(),
        track.isPVContributor(),
        track.hasITS(),
        track.hasTPC(),
        track.hasTRD(),
        track.hasTOF(),
        track.itsNCls(),
        track.tpcNClsFound(),
        track.tpcNClsCrossedRows(),
        track.tpcCrossedRowsOverFindableCls(),
        track.passedITSRefit(),
        track.passedTPCRefit(),
        track.itsChi2NCl(),
        track.tpcChi2NCl(),
        track.tpcNSigmaPr(),
        track.tofNSigmaPr(),
        track.tpcNSigmaDe(),
        track.tofNSigmaDe(),
        track.tpcNSigmaHe(),
        track.tofNSigmaHe(),
        track.mass() * track.mass(),
        event.centFT0C(),
        track.dcaXY(),
        track.dcaZ());
    }
  }

  //****************************************************************************************************

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>;

  using EventCandidatesCent = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>>; // aod::CentFV0As, aod::CentFT0Cs

  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCLfFullPr, aod::pidTOFFullPr, aod::pidTPCLfFullDe, aod::pidTOFFullDe, aod::pidTPCLfFullHe, aod::pidTOFFullHe, aod::TrackSelection, aod::TrackSelectionExtension, aod::TOFSignal, aod::pidTOFmass, aod::pidTOFbeta>>; // aod::ReducedTracks

  void processData(EventCandidates::iterator const& event, TrackCandidates const& tracks)
  {
    fillHistograms(event, tracks);
    fillTable(event, tracks);
  }
  PROCESS_SWITCH(NucleiHistTask, processData, "process data", true);

  void processDataCent(EventCandidatesCent::iterator const& event, TrackCandidates const& tracks)
  {
    fillHistograms(event, tracks);
    fillCentHistorgrams(event, tracks);
    fillTableCent(event, tracks);
  }
  PROCESS_SWITCH(NucleiHistTask, processDataCent, "process data with centralities", false);
};
//****************************************************************************************************

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<NucleiHistTask>(cfgc, TaskName{"nuclei-hist"})};
}
