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
/// \file LFNucleiBATask.cxx
///
/// \brief  Analysis task for the measurement of the coalescence parameter B2/B3 in pp collisions for (anti)deutheron/(anti)helium-3
///
/// \author Giovanni Malfattore <giovanni.malfattore@cern.ch> and Rutuparna Rath <rutuparna.rath@cern.ch>
///

#include "ReconstructionDataFormats/Track.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"

#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFNucleiTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct LFNucleiBATask {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry spectraGen{"spectraGen", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};

  Configurable<float> nsigmaTPCcut{"nsigmaTPCcut", 5.f, "Value of the Nsigma TPC cut"};
  Configurable<float> nsigmaTOFcut{"nsigmaTOFcut", 5.f, "Value of the Nsigma TOF cut"};
  Configurable<float> etaCut{"etaCut", 0.8f, "Value of the eta selection for spectra (default 0.8)"};
  Configurable<float> yCut{"yCut", 0.5f, "Value of the rapidity selection for spectra (default 0.5)"};
  Configurable<float> cfgCutVertex{"cfgCutVSertex", 10.0f, "Accepted z-vertex range"};
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0}, ""};
  Configurable<bool> enableEvTimeSplitting{"enableEvTimeSplitting", false, "Flag to enable histograms splitting depending on the Event Time used"};

  static constexpr int PDGPion = 211;
  static constexpr int PDGKaon = 321;
  static constexpr int PDGProton = 2212;
  static constexpr int PDGDeuteron = 1000010020;
  static constexpr int PDGHelium = 1000020030;
  static constexpr float fMassProton /* = o2::track::PID::getMass2Z(4));*/ = 0.938272088f;
  static constexpr float fMassDeuteron = 1.87561f;
  static constexpr float fMassHelium = 2.80839f;

  void init(o2::framework::InitContext&)
  {
    const AxisSpec pAxis{binsPt, "#it{p} (GeV/#it{c})"};
    const AxisSpec ptAxis{binsPt, "#it{p}_{T} (GeV/#it{c})"};

    if (doprocessData == true && doprocessMCReco == true) {
      LOG(fatal) << "Can't enable processData and processMCReco in the same time, pick one!";
    }

    histos.add<TH1>("event/h1VtxZ", "V_{z};V_{z} (in cm); counts", HistType::kTH1F, {{3000, -15, 15}});
    histos.add<TH1>("event/h1CentV0M", "V0M; Multiplicity; counts", HistType::kTH1F, {{27000, 0, 27000}});

    histos.add<TH1>("qa/h1TPCncr", "number of crossed rows in TPC; TPCncr; counts", HistType::kTH1F, {{150, 60, 170}});
    histos.add<TH1>("qa/h1rTPC", "ratio of ncr over findable in TPC; rTPC; counts", HistType::kTH1F, {{200, 0.9, 1.8}});
    histos.add<TH1>("qa/h1chi2ITS", "#chi^{2}_{ITS}/n_{ITS}; #chi^{2}_{ITS}/n_{ITS};counts", HistType::kTH1F, {{51, -0.5, 50.5}});
    histos.add<TH1>("qa/h1chi2TPC", "#chi^{2}_{TPC}/n_{TPC}; #chi^{2}_{TPC}/n_{TPC}; counts", HistType::kTH1F, {{11, -0.5, 10.5}});

    // trackQA
    histos.add<TH1>("tracks/h1Eta", "pseudoRapidity; #eta; counts", HistType::kTH1F, {{400, -1.0, 1.0}});
    histos.add<TH1>("tracks/h1VarPhi", "#phi; #phi; counts", HistType::kTH1F, {{63, 0.0, 6.3}});
    histos.add<TH2>("tracks/h2EtaVsPhi", "#eta vs #phi; #eta; #phi", HistType::kTH2F, {{400, -1.0, 1.0}, {63, 0.0, 6.3}});
    histos.add<TH1>("tracks/h1pT", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{1000, 0., 10}});
    histos.add<TH1>("tracks/h1p", "Track momentum; p (GeV/#it{c}); counts", HistType::kTH1F, {{1000, 0., 10.}});

    // tracks
    // DCAxy,z
    histos.add<TH1>("tracks/hDCAxy", "DCAxy; #DCAxy; counts", HistType::kTH1F, {{400, -1.0, 1.0}});
    histos.add<TH1>("tracks/hDCAz", "DCAz; #DCAz; counts", HistType::kTH1F, {{400, -1.0, 1.0}});
    histos.add<TH2>("tracks/hDCAxyVsDCAz", "DCAxy vs DCAz; DCAxy (cm); DCAz (cm)", HistType::kTH2F, {{400, -1.0, 1.0}, {400, -1.0, 1.0}});

    histos.add<TH2>("tracks/hDCAxyVsPt", "DCAxy vs Pt; #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});

    histos.add<TH2>("tracks/proton/dca/hDCAxyVsPtProton", "DCAxy vs Pt (p); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});
    histos.add<TH2>("tracks/deuteron/dca/hDCAxyVsPtDeuteron", "DCAxy vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});
    histos.add<TH2>("tracks/helium/dca/hDCAxyVsPtHelium", "DCAxy vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});

    histos.add<TH2>("tracks/proton/dca/hDCAxyVsPtantiProton", "DCAxy vs Pt (#bar{p}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});
    histos.add<TH2>("tracks/deuteron/dca/hDCAxyVsPtantiDeuteron", "DCAxy vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});
    histos.add<TH2>("tracks/helium/dca/hDCAxyVsPtantiHelium", "DCAxy vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});

    histos.add<TH2>("tracks/hDCAzVsPt", "DCAz vs Pt; #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});

    histos.add<TH2>("tracks/proton/dca/hDCAzVsPtProton", "DCAz vs Pt (p); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});
    histos.add<TH2>("tracks/deuteron/dca/hDCAzVsPtDeuteron", "DCAz vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});
    histos.add<TH2>("tracks/helium/dca/hDCAzVsPtHelium", "DCAz vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});

    histos.add<TH2>("tracks/proton/dca/hDCAzVsPtantiProton", "DCAz vs Pt (#bar{p}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});
    histos.add<TH2>("tracks/deuteron/dca/hDCAzVsPtantiDeuteron", "DCAz vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});
    histos.add<TH2>("tracks/helium/dca/hDCAzVsPtantiHelium", "DCAz vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});

    histos.add<TH1>("tracks/proton/h1ProtonSpectra", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
    histos.add<TH1>("tracks/deuteron/h1DeuteronSpectra", "#it{p}_{T} (d)", HistType::kTH1F, {ptAxis});
    histos.add<TH1>("tracks/helium/h1HeliumSpectra", "#it{p}_{T} (He)", HistType::kTH1F, {ptAxis});

    histos.add<TH1>("tracks/proton/h1antiProtonSpectra", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});
    histos.add<TH1>("tracks/deuteron/h1antiDeuteronSpectra", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
    histos.add<TH1>("tracks/helium/h1antiHeliumSpectra", "#it{p}_{T} (#bar{He})", HistType::kTH1F, {ptAxis});

    if (doprocessMCReco) {
      histos.add<TH1>("tracks/proton/h1ProtonSpectraTrue", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
      histos.add<TH1>("tracks/proton/h1ProtonSpectraTruePrim", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
      histos.add<TH1>("tracks/proton/h1ProtonSpectraTrueSec", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
      histos.add<TH1>("tracks/proton/h1ProtonSpectraTrueTransport", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});

      histos.add<TH1>("tracks/proton/h1antiProtonSpectraTrue", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
      histos.add<TH1>("tracks/proton/h1antiProtonSpectraTruePrim", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
      histos.add<TH1>("tracks/proton/h1antiProtonSpectraTrueSec", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
      histos.add<TH1>("tracks/proton/h1antiProtonSpectraTrueTransport", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});

      histos.add<TH1>("tracks/deuteron/h1DeuteronSpectraTrue", "#it{p}_{T} (d)", HistType::kTH1F, {ptAxis});
      histos.add<TH1>("tracks/deuteron/h1DeuteronSpectraTruePrim", "#it{p}_{T} (d)", HistType::kTH1F, {ptAxis});
      histos.add<TH1>("tracks/deuteron/h1DeuteronSpectraTrueSec", "#it{p}_{T} (d)", HistType::kTH1F, {ptAxis});
      histos.add<TH1>("tracks/deuteron/h1DeuteronSpectraTrueTransport", "#it{p}_{T} (d)", HistType::kTH1F, {ptAxis});

      histos.add<TH1>("tracks/deuteron/h1antiDeuteronSpectraTrue", "#it{p}_{T} (d)", HistType::kTH1F, {ptAxis});
      histos.add<TH1>("tracks/deuteron/h1antiDeuteronSpectraTruePrim", "#it{p}_{T} (d)", HistType::kTH1F, {ptAxis});
      histos.add<TH1>("tracks/deuteron/h1antiDeuteronSpectraTrueSec", "#it{p}_{T} (d)", HistType::kTH1F, {ptAxis});
      histos.add<TH1>("tracks/deuteron/h1antiDeuteronSpectraTrueTransport", "#it{p}_{T} (d)", HistType::kTH1F, {ptAxis});

      histos.add<TH1>("tracks/helium/h1HeliumSpectraTrue", "#it{p}_{T} (He)", HistType::kTH1F, {ptAxis});
      histos.add<TH1>("tracks/helium/h1HeliumSpectraTruePrim", "#it{p}_{T} (He)", HistType::kTH1F, {ptAxis});
      histos.add<TH1>("tracks/helium/h1HeliumSpectraTrueSec", "#it{p}_{T} (He)", HistType::kTH1F, {ptAxis});
      histos.add<TH1>("tracks/helium/h1HeliumSpectraTrueTransport", "#it{p}_{T} (He)", HistType::kTH1F, {ptAxis});

      histos.add<TH1>("tracks/helium/h1antiHeliumSpectraTrue", "#it{p}_{T} (He)", HistType::kTH1F, {ptAxis});
      histos.add<TH1>("tracks/helium/h1antiHeliumSpectraTruePrim", "#it{p}_{T} (He)", HistType::kTH1F, {ptAxis});
      histos.add<TH1>("tracks/helium/h1antiHeliumSpectraTrueSec", "#it{p}_{T} (He)", HistType::kTH1F, {ptAxis});
      histos.add<TH1>("tracks/helium/h1antiHeliumSpectraTrueTransport", "#it{p}_{T} (He)", HistType::kTH1F, {ptAxis});

      histos.add<TH2>("tracks/proton/dca/hDCAxyVsPtProtonTrue", "DCAxy vs Pt (p); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});
      histos.add<TH2>("tracks/deuteron/dca/hDCAxyVsPtDeuteronTrue", "DCAxy vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});
      histos.add<TH2>("tracks/helium/dca/hDCAxyVsPtHeliumTrue", "DCAxy vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});

      histos.add<TH2>("tracks/proton/dca/hDCAxyVsPtProtonTruePrim", "DCAxy vs Pt (p); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});
      histos.add<TH2>("tracks/proton/dca/hDCAxyVsPtProtonTrueSec", "DCAxy vs Pt (p); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});
      histos.add<TH2>("tracks/proton/dca/hDCAxyVsPtProtonTrueTransport", "DCAxy vs Pt (p); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});

      histos.add<TH2>("tracks/deuteron/dca/hDCAxyVsPtDeuteronTruePrim", "DCAxy vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});
      histos.add<TH2>("tracks/deuteron/dca/hDCAxyVsPtDeuteronTrueSec", "DCAxy vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});
      histos.add<TH2>("tracks/deuteron/dca/hDCAxyVsPtDeuteronTrueTransport", "DCAxy vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});

      histos.add<TH2>("tracks/helium/dca/hDCAxyVsPtHeliumTruePrim", "DCAxy vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});
      histos.add<TH2>("tracks/helium/dca/hDCAxyVsPtHeliumTrueSec", "DCAxy vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});
      histos.add<TH2>("tracks/helium/dca/hDCAxyVsPtHeliumTrueTransport", "DCAxy vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});

      histos.add<TH2>("tracks/proton/dca/hDCAxyVsPtantiProtonTrue", "DCAxy vs Pt (#bar{p}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});
      histos.add<TH2>("tracks/deuteron/dca/hDCAxyVsPtantiDeuteronTrue", "DCAxy vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});
      histos.add<TH2>("tracks/helium/dca/hDCAxyVsPtantiHeliumTrue", "DCAxy vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});

      histos.add<TH2>("tracks/proton/dca/hDCAxyVsPtantiProtonTruePrim", "DCAxy vs Pt (#bar{p}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});
      histos.add<TH2>("tracks/proton/dca/hDCAxyVsPtantiProtonTrueSec", "DCAxy vs Pt (#bar{p}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});
      histos.add<TH2>("tracks/proton/dca/hDCAxyVsPtantiProtonTrueTransport", "DCAxy vs Pt (#bar{p}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});

      histos.add<TH2>("tracks/deuteron/dca/hDCAxyVsPtantiDeuteronTruePrim", "DCAxy vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});
      histos.add<TH2>("tracks/deuteron/dca/hDCAxyVsPtantiDeuteronTrueSec", "DCAxy vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});
      histos.add<TH2>("tracks/deuteron/dca/hDCAxyVsPtantiDeuteronTrueTransport", "DCAxy vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});

      histos.add<TH2>("tracks/helium/dca/hDCAxyVsPtantiHeliumTruePrim", "DCAxy vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});
      histos.add<TH2>("tracks/helium/dca/hDCAxyVsPtantiHeliumTrueSec", "DCAxy vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});
      histos.add<TH2>("tracks/helium/dca/hDCAxyVsPtantiHeliumTrueTransport", "DCAxy vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});

      histos.add<TH2>("tracks/proton/dca/hDCAzVsPtProtonTrue", "DCAz vs Pt (p); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});
      histos.add<TH2>("tracks/deuteron/dca/hDCAzVsPtDeuteronTrue", "DCAz vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});
      histos.add<TH2>("tracks/helium/dca/hDCAzVsPtHeliumTrue", "DCAz vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});

      histos.add<TH2>("tracks/proton/dca/hDCAzVsPtantiProtonTrue", "DCAz vs Pt (#bar{p}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});
      histos.add<TH2>("tracks/deuteron/dca/hDCAzVsPtantiDeuteronTrue", "DCAz vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});
      histos.add<TH2>("tracks/helium/dca/hDCAzVsPtantiHeliumTrue", "DCAz vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});

      histos.add<TH2>("tracks/proton/dca/hDCAzVsPtProtonTruePrim", "DCAz vs Pt (p); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});
      histos.add<TH2>("tracks/proton/dca/hDCAzVsPtProtonTrueSec", "DCAz vs Pt (p); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});
      histos.add<TH2>("tracks/proton/dca/hDCAzVsPtProtonTrueTransport", "DCAz vs Pt (p); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});

      histos.add<TH2>("tracks/deuteron/dca/hDCAzVsPtDeuteronTruePrim", "DCAz vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});
      histos.add<TH2>("tracks/deuteron/dca/hDCAzVsPtDeuteronTrueSec", "DCAz vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});
      histos.add<TH2>("tracks/deuteron/dca/hDCAzVsPtDeuteronTrueTransport", "DCAz vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});

      histos.add<TH2>("tracks/helium/dca/hDCAzVsPtHeliumTruePrim", "DCAz vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});
      histos.add<TH2>("tracks/helium/dca/hDCAzVsPtHeliumTrueSec", "DCAz vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});
      histos.add<TH2>("tracks/helium/dca/hDCAzVsPtHeliumTrueTransport", "DCAz vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});

      histos.add<TH2>("tracks/proton/dca/hDCAzVsPtantiProtonTruePrim", "DCAz vs Pt (#bar{p}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});
      histos.add<TH2>("tracks/proton/dca/hDCAzVsPtantiProtonTrueSec", "DCAz vs Pt (#bar{p}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});
      histos.add<TH2>("tracks/proton/dca/hDCAzVsPtantiProtonTrueTransport", "DCAz vs Pt (#bar{p}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});

      histos.add<TH2>("tracks/deuteron/dca/hDCAzVsPtantiDeuteronTruePrim", "DCAz vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});
      histos.add<TH2>("tracks/deuteron/dca/hDCAzVsPtantiDeuteronTrueSec", "DCAz vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});
      histos.add<TH2>("tracks/deuteron/dca/hDCAzVsPtantiDeuteronTrueTransport", "DCAz vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});

      histos.add<TH2>("tracks/helium/dca/hDCAzVsPtantiHeliumTruePrim", "DCAz vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});
      histos.add<TH2>("tracks/helium/dca/hDCAzVsPtantiHeliumTrueSec", "DCAz vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});
      histos.add<TH2>("tracks/helium/dca/hDCAzVsPtantiHeliumTrueTransport", "DCAz vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {400, -1.0, 1.0}});
    }

    //  Bethe-Bloch TPC distribution and Beta vs pT TOF distribution
    histos.add<TH2>("tracks/h2TPCsignVsTPCmomentum", "TPC <-dE/dX> vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TPC <-dE/dx> (a.u.)", HistType::kTH2F, {{1000, -5.f, 5.f}, {8000, 0.0, 1E3}});
    histos.add<TH2>("tracks/h2TOFbetaVsP", "#TOF beta vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TOF #beta", HistType::kTH2F, {{1000, -5.f, 5.f}, {1200, 0.0, 1.2}});

    //  TPCExpSignal histograms
    histos.add<TH2>("tracks/proton/h2ProtonTPCExpSignalDiffVsPt", "TPC <-dE/dX> - Exp <-dE/dX> (p) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TPC <-dE/dX> ExpDiff (p)", HistType::kTH2F, {{ptAxis}, {16000, -800, 800.}});
    histos.add<TH2>("tracks/proton/h2antiProtonTPCExpSignalDiffVsPt", "TPC <-dE/dX> - Exp <-dE/dX> (#bar{p}) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TPC <-dE/dX> ExpDiff (#bar{p})", HistType::kTH2F, {{ptAxis}, {16000, -800, 800.}});
    histos.add<TH2>("tracks/deuteron/h2DeuteronTPCExpSignalDiffVsPt", "TPC <-dE/dX> - Exp <-dE/dX> (d) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TPC <-dE/dX> ExpDiff (d)", HistType::kTH2F, {{ptAxis}, {16000, -800, 800.}});
    histos.add<TH2>("tracks/deuteron/h2antiDeuteronTPCExpSignalDiffVsPt", "TPC <-dE/dX> - Exp <-dE/dX> (#bar{d}) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TPC <-dE/dX> ExpDiff (#bar{d})", HistType::kTH2F, {{ptAxis}, {16000, -800, 800.}});

    //  NSigmasTPC histograms
    histos.add<TH2>("tracks/pion/h2PionVspTNSigmaTPC", "NSigmaTPC(pi) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptAxis}, {2000, -100, 100}});
    histos.add<TH2>("tracks/kaon/h2KaonVspTNSigmaTPC", "NSigmaTPC(Ka) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptAxis}, {2000, -100, 100}});
    histos.add<TH2>("tracks/proton/h2ProtonVspTNSigmaTPC", "NSigmaTPC(p) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptAxis}, {2000, -100, 100}});
    histos.add<TH2>("tracks/deuteron/h2DeuteronVspTNSigmaTPC", "NSigmaTPC(d) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptAxis}, {2000, -100, 100}});
    histos.add<TH2>("tracks/helium/h2HeliumVspTNSigmaTPC", "NSigmaTPC(He) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptAxis}, {2000, -100, 100}});
    histos.add<TH2>("tracks/proton/h2antiProtonVspTNSigmaTPC", "NSigmaTPC(#bar{p}) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptAxis}, {2000, -100, 100}});
    histos.add<TH2>("tracks/deuteron/h2antiDeuteronVspTNSigmaTPC", "NSigmaTPC(#bar{d}) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptAxis}, {2000, -100, 100}});
    histos.add<TH2>("tracks/helium/h2antiHeliumVspTNSigmaTPC", "NSigmaTPC(#bar{He}) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptAxis}, {2000, -100, 100}});

    //  TOFExpSignal histograms
    histos.add<TH2>("tracks/proton/h2ProtonTOFExpSignalDiffVsPt", "TOF t - t_{exp}(p) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (p)", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
    histos.add<TH2>("tracks/proton/h2antiProtonTOFExpSignalDiffVsPt", "TOF t - t_{exp}(#bar{p}) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (#bar{p})", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
    histos.add<TH2>("tracks/deuteron/h2DeuteronTOFExpSignalDiffVsPt", "TOF t - t_{exp}(d) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (d)", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
    histos.add<TH2>("tracks/deuteron/h2antiDeuteronTOFExpSignalDiffVsPt", "TOF t - t_{exp}(#bar{d}) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (#bar{d})", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});

    histos.add<TH2>("tracks/proton/h2ProtonTOFExpSignalDiffVsPtCut", "TOF t - t_{exp}(p) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (p)", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
    histos.add<TH2>("tracks/proton/h2antiProtonTOFExpSignalDiffVsPtCut", "TOF t - t_{exp}(#bar{p}) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (#bar{p})", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
    histos.add<TH2>("tracks/deuteron/h2DeuteronTOFExpSignalDiffVsPtCut", "TOF t - t_{exp}(d) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (d)", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
    histos.add<TH2>("tracks/deuteron/h2antiDeuteronTOFExpSignalDiffVsPtCut", "TOF t - t_{exp}(#bar{d}) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (#bar{d})", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});

    // NSigmaTOF histograms
    histos.add<TH2>("tracks/pion/h2PionVspTNSigmaTOF", "NSigmaTOF(pi) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {2000, -100, 100}});
    histos.add<TH2>("tracks/kaon/h2KaonVspTNSigmaTOF", "NSigmaTOF(Ka) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {2000, -100, 100}});
    histos.add<TH2>("tracks/proton/h2ProtonVspTNSigmaTOF", "NSigmaTOF(p) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {2000, -100, 100}});
    histos.add<TH2>("tracks/deuteron/h2DeuteronVspTNSigmaTOF", "NSigmaTOF(d) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {2000, -100, 100}});
    histos.add<TH2>("tracks/helium/h2HeliumVspTNSigmaTOF", "NSigmaTOF(He) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {2000, -100, 100}});
    histos.add<TH2>("tracks/proton/h2antiProtonVspTNSigmaTOF", "NSigmaTOF(#bar{p}) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {2000, -100, 100}});
    histos.add<TH2>("tracks/deuteron/h2antiDeuteronVspTNSigmaTOF", "NSigmaTOF(#bar{d}) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {2000, -100, 100}});
    histos.add<TH2>("tracks/helium/h2antiHeliumVspTNSigmaTOF", "NSigmaTOF(#bar{He}) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {2000, -100, 100}});

    // NSigmaTPC vs NSigmaTOF histograms
    histos.add<TH3>("tracks/proton/h3ProtonNSigmaTPCvsNSigmaTOFvsPt", "NSigmaTPC (p) vs NSigmaTOF(p); NSigmaTPC; NSigmaTOF; #it{p}_{T} (GeV/#it{c})", HistType::kTH3F, {{400, -20, 20}, {400, -20, 20}, {ptAxis}});
    histos.add<TH3>("tracks/deuteron/h3DeuteronNSigmaTPCvsNSigmaTOFvsPt", "NSigmaTPC(d) vs NSigmaTOF(d); NSigmaTPC; NSigmaTOF; #it{p}_{T} (GeV/#it{c})", HistType::kTH3F, {{400, -20, 20}, {400, -20, 20}, {ptAxis}});
    histos.add<TH3>("tracks/proton/h3antiProtonNSigmaTPCvsNSigmaTOFvsPt", "NSigmaTPC (#bar{p}) vs NSigmaTOF(#bar{p}); NSigmaTPC; NSigmaTOF; #it{p}_{T} (GeV/#it{c})", HistType::kTH3F, {{400, -20, 20}, {400, -20, 20}, {ptAxis}});
    histos.add<TH3>("tracks/deuteron/h3antiDeuteronNSigmaTPCvsNSigmaTOFvsPt", "NSigmaTPC (#bar{d}) vs NSigmaTOF(#bar{d}); NSigmaTPC; NSigmaTOF; #it{p}_{T} (GeV/#it{c})", HistType::kTH3F, {{400, -20, 20}, {400, -20, 20}, {ptAxis}});

    // TOF mass histograms
    histos.add<TH2>("tracks/h2TOFmassVsPt", "h2TOFmassVsPt; TOFmass; #it{p}_{T} (GeV)", HistType::kTH2F, {{600, 0., 3.}, {500, 0., 5.}});

    // TOF mass squared histograms
    histos.add<TH2>("tracks/proton/h2TOFmass2ProtonVsPt", "#Delta M^{2} (p) vs #it{p}_{T}; #Delta M^{2} (p); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{600, -3., 3.}, {800, 0., 8.}});
    histos.add<TH2>("tracks/deuteron/h2TOFmass2DeuteronVsPt", "#Delta M^{2} (d) vs #it{p}_{T}; #Delta M^{2} (d); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{1000, -5., 5.}, {800, 0., 8.}});
    histos.add<TH2>("tracks/helium/h2TOFmass2HeliumVsPt", "#Delta M^{2} (He) vs #it{p}_{T}; #Delta M^{2} (He); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{1800, -9., 9.}, {800, 0., 8.}});

    histos.add<TH2>("tracks/proton/h2TOFmass2antiProtonVsPt", "#Delta M^{2} (#bar{p}) vs #it{p}_{T}; #Delta M^{2} (#bar{p}); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{600, -3., 3.}, {800, 0., 8.}});
    histos.add<TH2>("tracks/deuteron/h2TOFmass2antiDeuteronVsPt", "#Delta M^{2} (#bar{d}) vs #it{p}_{T}; #Delta M^{2} (#bar{d}); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{1000, -5., 5.}, {800, 0., 8.}});
    histos.add<TH2>("tracks/helium/h2TOFmass2antiHeliumVsPt", "#Delta M^{2} (#bar{He}) vs #it{p}_{T}; #Delta M^{2} (#bar{He}); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{1800, -9., 9.}, {800, 0., 8.}});

    // TOF EvTime Splitting plots
    if (enableEvTimeSplitting) {
      //  Bethe-Bloch TPC distribution - TOF EvTime Splitted
      histos.add<TH2>("tracks/evtime/fill/h2TPCsignVsTPCmomentum", "TPC <-dE/dX> vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TPC <-dE/dx> (a.u.)", HistType::kTH2F, {{1000, -5.f, 5.f}, {8000, 0.0, 1E3}});
      histos.add<TH2>("tracks/evtime/tof/h2TPCsignVsTPCmomentum", "TPC <-dE/dX> vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TPC <-dE/dx> (a.u.)", HistType::kTH2F, {{1000, -5.f, 5.f}, {8000, 0.0, 1E3}});
      histos.add<TH2>("tracks/evtime/ft0/h2TPCsignVsTPCmomentum", "TPC <-dE/dX> vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TPC <-dE/dx> (a.u.)", HistType::kTH2F, {{1000, -5.f, 5.f}, {8000, 0.0, 1E3}});
      histos.add<TH2>("tracks/evtime/ft0tof/h2TPCsignVsTPCmomentum", "TPC <-dE/dX> vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TPC <-dE/dx> (a.u.)", HistType::kTH2F, {{1000, -5.f, 5.f}, {8000, 0.0, 1E3}});

      // Beta vs pT TOF distribution - TOF EvTime Splitted
      histos.add<TH2>("tracks/evtime/fill/h2TOFbetaVsP", "#TOF beta vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TOF #beta", HistType::kTH2F, {{1000, -5.f, 5.f}, {1200, 0.0, 1.2}});
      histos.add<TH2>("tracks/evtime/tof/h2TOFbetaVsP", "#TOF beta vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TOF #beta", HistType::kTH2F, {{1000, -5.f, 5.f}, {1200, 0.0, 1.2}});
      histos.add<TH2>("tracks/evtime/ft0/h2TOFbetaVsP", "#TOF beta vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TOF #beta", HistType::kTH2F, {{1000, -5.f, 5.f}, {1200, 0.0, 1.2}});
      histos.add<TH2>("tracks/evtime/ft0tof/h2TOFbetaVsP", "#TOF beta vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TOF #beta", HistType::kTH2F, {{1000, -5.f, 5.f}, {1200, 0.0, 1.2}});

      //  TOFExpSignal histograms - TOF EvTime Splitted
      histos.add<TH2>("tracks/evtime/fill/proton/h2ProtonTOFExpSignalDiffVsPt", "TOF t - t_{exp}(p) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (p)", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
      histos.add<TH2>("tracks/evtime/fill/proton/h2antiProtonTOFExpSignalDiffVsPt", "TOF t - t_{exp}(#bar{p}) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (#bar{p})", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
      histos.add<TH2>("tracks/evtime/fill/deuteron/h2DeuteronTOFExpSignalDiffVsPt", "TOF t - t_{exp}(d) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (d)", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
      histos.add<TH2>("tracks/evtime/fill/deuteron/h2antiDeuteronTOFExpSignalDiffVsPt", "TOF t - t_{exp}(#bar{d}) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (#bar{d})", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});

      histos.add<TH2>("tracks/evtime/tof/proton/h2ProtonTOFExpSignalDiffVsPt", "TOF t - t_{exp}(p) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (p)", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
      histos.add<TH2>("tracks/evtime/tof/proton/h2antiProtonTOFExpSignalDiffVsPt", "TOF t - t_{exp}(#bar{p}) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (#bar{p})", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
      histos.add<TH2>("tracks/evtime/tof/deuteron/h2DeuteronTOFExpSignalDiffVsPt", "TOF t - t_{exp}(d) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (d)", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
      histos.add<TH2>("tracks/evtime/tof/deuteron/h2antiDeuteronTOFExpSignalDiffVsPt", "TOF t - t_{exp}(#bar{d}) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (#bar{d})", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});

      histos.add<TH2>("tracks/evtime/ft0/proton/h2ProtonTOFExpSignalDiffVsPt", "TOF t - t_{exp}(p) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (p)", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
      histos.add<TH2>("tracks/evtime/ft0/proton/h2antiProtonTOFExpSignalDiffVsPt", "TOF t - t_{exp}(#bar{p}) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (#bar{p})", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
      histos.add<TH2>("tracks/evtime/ft0/deuteron/h2DeuteronTOFExpSignalDiffVsPt", "TOF t - t_{exp}(d) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (d)", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
      histos.add<TH2>("tracks/evtime/ft0/deuteron/h2antiDeuteronTOFExpSignalDiffVsPt", "TOF t - t_{exp}(#bar{d}) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (#bar{d})", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});

      histos.add<TH2>("tracks/evtime/ft0tof/proton/h2ProtonTOFExpSignalDiffVsPt", "TOF t - t_{exp}(p) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (p)", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
      histos.add<TH2>("tracks/evtime/ft0tof/proton/h2antiProtonTOFExpSignalDiffVsPt", "TOF t - t_{exp}(#bar{p}) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (#bar{p})", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
      histos.add<TH2>("tracks/evtime/ft0tof/deuteron/h2DeuteronTOFExpSignalDiffVsPt", "TOF t - t_{exp}(d) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (d)", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
      histos.add<TH2>("tracks/evtime/ft0tof/deuteron/h2antiDeuteronTOFExpSignalDiffVsPt", "TOF t - t_{exp}(#bar{d}) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (#bar{d})", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});

      // NSigmaTOF histograms - TOF EvTime Splitted
      histos.add<TH2>("tracks/evtime/fill/proton/h2ProtonVspTNSigmaTOF", "NSigmaTOF(p) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {2000, -100, 100}});
      histos.add<TH2>("tracks/evtime/fill/deuteron/h2DeuteronVspTNSigmaTOF", "NSigmaTOF(d) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {2000, -100, 100}});
      histos.add<TH2>("tracks/evtime/fill/proton/h2antiProtonVspTNSigmaTOF", "NSigmaTOF(#bar{p}) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {2000, -100, 100}});
      histos.add<TH2>("tracks/evtime/fill/deuteron/h2antiDeuteronVspTNSigmaTOF", "NSigmaTOF(#bar{d}) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {2000, -100, 100}});

      histos.add<TH2>("tracks/evtime/tof/proton/h2ProtonVspTNSigmaTOF", "NSigmaTOF(p) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {2000, -100, 100}});
      histos.add<TH2>("tracks/evtime/tof/deuteron/h2DeuteronVspTNSigmaTOF", "NSigmaTOF(d) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {2000, -100, 100}});
      histos.add<TH2>("tracks/evtime/tof/proton/h2antiProtonVspTNSigmaTOF", "NSigmaTOF(#bar{p}) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {2000, -100, 100}});
      histos.add<TH2>("tracks/evtime/tof/deuteron/h2antiDeuteronVspTNSigmaTOF", "NSigmaTOF(#bar{d}) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {2000, -100, 100}});

      histos.add<TH2>("tracks/evtime/ft0/proton/h2ProtonVspTNSigmaTOF", "NSigmaTOF(p) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {2000, -100, 100}});
      histos.add<TH2>("tracks/evtime/ft0/deuteron/h2DeuteronVspTNSigmaTOF", "NSigmaTOF(d) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {2000, -100, 100}});
      histos.add<TH2>("tracks/evtime/ft0/proton/h2antiProtonVspTNSigmaTOF", "NSigmaTOF(#bar{p}) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {2000, -100, 100}});
      histos.add<TH2>("tracks/evtime/ft0/deuteron/h2antiDeuteronVspTNSigmaTOF", "NSigmaTOF(#bar{d}) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {2000, -100, 100}});

      histos.add<TH2>("tracks/evtime/ft0tof/proton/h2ProtonVspTNSigmaTOF", "NSigmaTOF(p) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {2000, -100, 100}});
      histos.add<TH2>("tracks/evtime/ft0tof/deuteron/h2DeuteronVspTNSigmaTOF", "NSigmaTOF(d) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {2000, -100, 100}});
      histos.add<TH2>("tracks/evtime/ft0tof/proton/h2antiProtonVspTNSigmaTOF", "NSigmaTOF(#bar{p}) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {2000, -100, 100}});
      histos.add<TH2>("tracks/evtime/ft0tof/deuteron/h2antiDeuteronVspTNSigmaTOF", "NSigmaTOF(#bar{d}) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {2000, -100, 100}});

      // NSigmaTPC vs NSigmaTOF histograms - TOF EvTime Splitted
      histos.add<TH3>("tracks/evtime/fill/proton/h3ProtonNSigmaTPCvsNSigmaTOFvsPt", "NSigmaTPC (p) vs NSigmaTOF(p); NSigmaTPC; NSigmaTOF; #it{p}_{T} (GeV/#it{c})", HistType::kTH3F, {{400, -20, 20}, {400, -20, 20}, {ptAxis}});
      histos.add<TH3>("tracks/evtime/fill/deuteron/h3DeuteronNSigmaTPCvsNSigmaTOFvsPt", "NSigmaTPC(d) vs NSigmaTOF(d); NSigmaTPC; NSigmaTOF; #it{p}_{T} (GeV/#it{c})", HistType::kTH3F, {{400, -20, 20}, {400, -20, 20}, {ptAxis}});
      histos.add<TH3>("tracks/evtime/fill/proton/h3antiProtonNSigmaTPCvsNSigmaTOFvsPt", "NSigmaTPC(#bar{p}) vs NSigmaTOF(#bar{p}); NSigmaTPC; NSigmaTOF; #it{p}_{T} (GeV/#it{c})", HistType::kTH3F, {{400, -20, 20}, {400, -20, 20}, {ptAxis}});
      histos.add<TH3>("tracks/evtime/fill/deuteron/h3antiDeuteronNSigmaTPCvsNSigmaTOFvsPt", "NSigmaTPC(#bar{d}) vs NSigmaTOF(#bar{p}); NSigmaTPC; NSigmaTOF; #it{p}_{T} (GeV/#it{c})", HistType::kTH3F, {{400, -20, 20}, {400, -20, 20}, {ptAxis}});

      histos.add<TH3>("tracks/evtime/tof/proton/h3ProtonNSigmaTPCvsNSigmaTOFvsPt", "NSigmaTPC (p) vs NSigmaTOF(p); NSigmaTPC; NSigmaTOF; #it{p}_{T} (GeV/#it{c})", HistType::kTH3F, {{400, -20, 20}, {400, -20, 20}, {ptAxis}});
      histos.add<TH3>("tracks/evtime/tof/deuteron/h3DeuteronNSigmaTPCvsNSigmaTOFvsPt", "NSigmaTPC(d) vs NSigmaTOF(d); NSigmaTPC; NSigmaTOF; #it{p}_{T} (GeV/#it{c})", HistType::kTH3F, {{400, -20, 20}, {400, -20, 20}, {ptAxis}});
      histos.add<TH3>("tracks/evtime/tof/proton/h3antiProtonNSigmaTPCvsNSigmaTOFvsPt", "NSigmaTPC(#bar{p}) vs NSigmaTOF(#bar{p}); NSigmaTPC; NSigmaTOF; #it{p}_{T} (GeV/#it{c})", HistType::kTH3F, {{400, -20, 20}, {400, -20, 20}, {ptAxis}});
      histos.add<TH3>("tracks/evtime/tof/deuteron/h3antiDeuteronNSigmaTPCvsNSigmaTOFvsPt", "NSigmaTPC(#bar{d}) vs NSigmaTOF(#bar{p}); NSigmaTPC; NSigmaTOF; #it{p}_{T} (GeV/#it{c})", HistType::kTH3F, {{400, -20, 20}, {400, -20, 20}, {ptAxis}});

      histos.add<TH3>("tracks/evtime/ft0/proton/h3ProtonNSigmaTPCvsNSigmaTOFvsPt", "NSigmaTPC (p) vs NSigmaTOF(p); NSigmaTPC; NSigmaTOF; #it{p}_{T} (GeV/#it{c})", HistType::kTH3F, {{400, -20, 20}, {400, -20, 20}, {ptAxis}});
      histos.add<TH3>("tracks/evtime/ft0/deuteron/h3DeuteronNSigmaTPCvsNSigmaTOFvsPt", "NSigmaTPC(d) vs NSigmaTOF(d); NSigmaTPC; NSigmaTOF; #it{p}_{T} (GeV/#it{c})", HistType::kTH3F, {{400, -20, 20}, {400, -20, 20}, {ptAxis}});
      histos.add<TH3>("tracks/evtime/ft0/proton/h3antiProtonNSigmaTPCvsNSigmaTOFvsPt", "NSigmaTPC(#bar{p}) vs NSigmaTOF(#bar{p}); NSigmaTPC; NSigmaTOF; #it{p}_{T} (GeV/#it{c})", HistType::kTH3F, {{400, -20, 20}, {400, -20, 20}, {ptAxis}});
      histos.add<TH3>("tracks/evtime/ft0/deuteron/h3antiDeuteronNSigmaTPCvsNSigmaTOFvsPt", "NSigmaTPC(#bar{d}) vs NSigmaTOF(#bar{p}); NSigmaTPC; NSigmaTOF; #it{p}_{T} (GeV/#it{c})", HistType::kTH3F, {{400, -20, 20}, {400, -20, 20}, {ptAxis}});

      histos.add<TH3>("tracks/evtime/ft0tof/proton/h3ProtonNSigmaTPCvsNSigmaTOFvsPt", "NSigmaTPC (p) vs NSigmaTOF(p); NSigmaTPC; NSigmaTOF; #it{p}_{T} (GeV/#it{c})", HistType::kTH3F, {{400, -20, 20}, {400, -20, 20}, {ptAxis}});
      histos.add<TH3>("tracks/evtime/ft0tof/deuteron/h3DeuteronNSigmaTPCvsNSigmaTOFvsPt", "NSigmaTPC(d) vs NSigmaTOF(d); NSigmaTPC; NSigmaTOF; #it{p}_{T} (GeV/#it{c})", HistType::kTH3F, {{400, -20, 20}, {400, -20, 20}, {ptAxis}});
      histos.add<TH3>("tracks/evtime/ft0tof/proton/h3antiProtonNSigmaTPCvsNSigmaTOFvsPt", "NSigmaTPC(#bar{p}) vs NSigmaTOF(#bar{p}); NSigmaTPC; NSigmaTOF; #it{p}_{T} (GeV/#it{c})", HistType::kTH3F, {{400, -20, 20}, {400, -20, 20}, {ptAxis}});
      histos.add<TH3>("tracks/evtime/ft0tof/deuteron/h3antiDeuteronNSigmaTPCvsNSigmaTOFvsPt", "NSigmaTPC(#bar{d}) vs NSigmaTOF(#bar{p}); NSigmaTPC; NSigmaTOF; #it{p}_{T} (GeV/#it{c})", HistType::kTH3F, {{400, -20, 20}, {400, -20, 20}, {ptAxis}});

      // TOF mass histograms - TOF EvTime Splitted
      histos.add<TH2>("tracks/evtime/fill/h2TOFmassVsPt", "h2TOFmassVsPt; TOFmass; #it{p}_{T} (GeV)", HistType::kTH2F, {{600, 0., 3.}, {500, 0., 5.}});
      histos.add<TH2>("tracks/evtime/tof/h2TOFmassVsPt", "h2TOFmassVsPt; TOFmass; #it{p}_{T} (GeV)", HistType::kTH2F, {{600, 0., 3.}, {500, 0., 5.}});
      histos.add<TH2>("tracks/evtime/ft0/h2TOFmassVsPt", "h2TOFmassVsPt; TOFmass; #it{p}_{T} (GeV)", HistType::kTH2F, {{600, 0., 3.}, {500, 0., 5.}});
      histos.add<TH2>("tracks/evtime/ft0tof/h2TOFmassVsPt", "h2TOFmassVsPt; TOFmass; #it{p}_{T} (GeV)", HistType::kTH2F, {{600, 0., 3.}, {500, 0., 5.}});

      // TOF mass squared histograms - TOF EvTime Splitted
      histos.add<TH2>("tracks/evtime/fill/proton/h2TOFmass2ProtonVsPt", "#Delta M^{2} (p) vs #it{p}_{T}; #Delta M^{2} (p); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{600, -3., 3.}, {800, 0., 8.}});
      histos.add<TH2>("tracks/evtime/fill/deuteron/h2TOFmass2DeuteronVsPt", "#Delta M^{2} (d) vs #it{p}_{T}; #Delta M^{2} (d); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{1000, -5., 5.}, {800, 0., 8.}});

      histos.add<TH2>("tracks/evtime/tof/proton/h2TOFmass2ProtonVsPt", "#Delta M^{2} (p) vs #it{p}_{T}; #Delta M^{2} (p); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{600, -3., 3.}, {800, 0., 8.}});
      histos.add<TH2>("tracks/evtime/tof/deuteron/h2TOFmass2DeuteronVsPt", "#Delta M^{2} (d) vs #it{p}_{T}; #Delta M^{2} (d); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{1000, -5., 5.}, {800, 0., 8.}});

      histos.add<TH2>("tracks/evtime/ft0/proton/h2TOFmass2ProtonVsPt", "#Delta M^{2} (p) vs #it{p}_{T}; #Delta M^{2} (p); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{600, -3., 3.}, {800, 0., 8.}});
      histos.add<TH2>("tracks/evtime/ft0/deuteron/h2TOFmass2DeuteronVsPt", "#Delta M^{2} (d) vs #it{p}_{T}; #Delta M^{2} (d); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{1000, -5., 5.}, {800, 0., 8.}});

      histos.add<TH2>("tracks/evtime/ft0tof/proton/h2TOFmass2ProtonVsPt", "#Delta M^{2} (p) vs #it{p}_{T}; #Delta M^{2} (p); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{600, -3., 3.}, {800, 0., 8.}});
      histos.add<TH2>("tracks/evtime/ft0tof/deuteron/h2TOFmass2DeuteronVsPt", "#Delta M^{2} (d) vs #it{p}_{T}; #Delta M^{2} (d); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{1000, -5., 5.}, {800, 0., 8.}});

      histos.add<TH2>("tracks/evtime/fill/proton/h2TOFmass2antiProtonVsPt", "#Delta M^{2} (#bar{p}) vs #it{p}_{T}; #Delta M^{2} (#bar{p}); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{600, -3., 3.}, {800, 0., 8.}});
      histos.add<TH2>("tracks/evtime/fill/deuteron/h2TOFmass2antiDeuteronVsPt", "#Delta M^{2} (#bar{d}) vs #it{p}_{T}; #Delta M^{2} (#bar{d}); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{1000, -5., 5.}, {800, 0., 8.}});

      histos.add<TH2>("tracks/evtime/tof/proton/h2TOFmass2antiProtonVsPt", "#Delta M^{2} (#bar{p}) vs #it{p}_{T}; #Delta M^{2} (#bar{p}); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{600, -3., 3.}, {800, 0., 8.}});
      histos.add<TH2>("tracks/evtime/tof/deuteron/h2TOFmass2antiDeuteronVsPt", "#Delta M^{2} (#bar{d}) vs #it{p}_{T}; #Delta M^{2} (#bar{d}); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{1000, -5., 5.}, {800, 0., 8.}});

      histos.add<TH2>("tracks/evtime/ft0/proton/h2TOFmass2antiProtonVsPt", "#Delta M^{2} (#bar{p}) vs #it{p}_{T}; #Delta M^{2} (#bar{p}); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{600, -3., 3.}, {800, 0., 8.}});
      histos.add<TH2>("tracks/evtime/ft0/deuteron/h2TOFmass2antiDeuteronVsPt", "#Delta M^{2} (#bar{d}) vs #it{p}_{T}; #Delta M^{2} (#bar{d}); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{1000, -5., 5.}, {800, 0., 8.}});

      histos.add<TH2>("tracks/evtime/ft0tof/proton/h2TOFmass2antiProtonVsPt", "#Delta M^{2} (#bar{p}) vs #it{p}_{T}; #Delta M^{2} (#bar{p}); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{600, -3., 3.}, {800, 0., 8.}});
      histos.add<TH2>("tracks/evtime/ft0tof/deuteron/h2TOFmass2antiDeuteronVsPt", "#Delta M^{2} (#bar{d}) vs #it{p}_{T}; #Delta M^{2} (#bar{d}); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{1000, -5., 5.}, {800, 0., 8.}});
    }

    // MC histograms  -   all, primary, sec. from weak decay, sec. from material
    histos.add("spectraGen/histGenVetxZ", "PosZ generated events", HistType::kTH1F, {{2000, -20.f, 20.f, "Vertex Z (cm)"}});

    histos.add("spectraGen/histGenPtPion", "generated particles", HistType::kTH1F, {ptAxis});
    histos.add("spectraGen/histGenPtPionPrim", "generated particles", HistType::kTH1F, {ptAxis});
    histos.add("spectraGen/histGenPtPionSec", "generated particles", HistType::kTH1F, {ptAxis});
    histos.add("spectraGen/histSecTransportPtPion", "generated particles", HistType::kTH1F, {ptAxis});

    histos.add("spectraGen/histGenPtKaon", "generated particles", HistType::kTH1F, {ptAxis});
    histos.add("spectraGen/histGenPtKaonPrim", "generated particles", HistType::kTH1F, {ptAxis});
    histos.add("spectraGen/histGenPtKaonSec", "generated particles", HistType::kTH1F, {ptAxis});
    histos.add("spectraGen/histSecTransportPtKaon", "generated particles", HistType::kTH1F, {ptAxis});

    histos.add("spectraGen/histGenPtProton", "generated particles", HistType::kTH1F, {ptAxis});
    histos.add("spectraGen/histGenPtProtonPrim", "generated particles", HistType::kTH1F, {ptAxis});
    histos.add("spectraGen/histGenPtProtonSec", "generated particles", HistType::kTH1F, {ptAxis});
    histos.add("spectraGen/histSecTransportPtProton", "generated particles", HistType::kTH1F, {ptAxis});

    histos.add("spectraGen/histGenPtantiProton", "generated particles", HistType::kTH1F, {ptAxis});
    histos.add("spectraGen/histGenPtantiProtonPrim", "generated particles", HistType::kTH1F, {ptAxis});
    histos.add("spectraGen/histGenPtantiProtonSec", "generated particles", HistType::kTH1F, {ptAxis});
    histos.add("spectraGen/histSecTransportPtantiProton", "generated particles", HistType::kTH1F, {ptAxis});

    histos.add("spectraGen/histGenPtD", "generated particles", HistType::kTH1F, {ptAxis});
    histos.add("spectraGen/histGenPtDPrim", "generated particles", HistType::kTH1F, {ptAxis});
    histos.add("spectraGen/histGenPtDSec", "generated particles", HistType::kTH1F, {ptAxis});
    histos.add("spectraGen/histSecTransportPtD", "generated particles", HistType::kTH1F, {ptAxis});

    histos.add("spectraGen/histGenPtantiD", "generated particles", HistType::kTH1F, {ptAxis});
    histos.add("spectraGen/histGenPtantiDPrim", "generated particles", HistType::kTH1F, {ptAxis});
    histos.add("spectraGen/histGenPtantiDSec", "generated particles", HistType::kTH1F, {ptAxis});
    histos.add("spectraGen/histSecTransportPtantiD", "generated particles", HistType::kTH1F, {ptAxis});

    histos.add("spectraGen/histGenPtHe", "generated particles", HistType::kTH1F, {ptAxis});
    histos.add("spectraGen/histGenPtHePrim", "generated particles", HistType::kTH1F, {ptAxis});
    histos.add("spectraGen/histGenPtHeSec", "generated particles", HistType::kTH1F, {ptAxis});
    histos.add("spectraGen/histSecTransportPtHe", "generated particles", HistType::kTH1F, {ptAxis});

    histos.add("spectraGen/histGenPtantiHe", "generated particles", HistType::kTH1F, {ptAxis});
    histos.add("spectraGen/histGenPtantiHePrim", "generated particles", HistType::kTH1F, {ptAxis});
    histos.add("spectraGen/histGenPtantiHeSec", "generated particles", HistType::kTH1F, {ptAxis});
    histos.add("spectraGen/histSecTransportPtantiHe", "generated particles", HistType::kTH1F, {ptAxis});
  }

  template <bool IsMC, typename CollisionType, typename TracksType>
  void fillHistograms(const CollisionType& event, const TracksType& tracks)
  {

    float gamma = 0., massTOF = 0.;

    // Event histos fill
    histos.fill(HIST("event/h1VtxZ"), event.posZ());
    histos.fill(HIST("event/h1CentV0M"), event.v0m());

    for (auto& track : tracks) {
      // LOG(info)<<"\n collisionId ============>"<<track.collisionId();

      // QA histos fill
      histos.fill(HIST("qa/h1TPCncr"), track.ncrTPC());
      histos.fill(HIST("qa/h1rTPC"), track.rTPC());
      histos.fill(HIST("qa/h1chi2ITS"), track.chi2TPC());
      histos.fill(HIST("qa/h1chi2TPC"), track.chi2ITS());

      // Tracks DCA histos fill
      histos.fill(HIST("tracks/hDCAxy"), track.dcaxy());
      histos.fill(HIST("tracks/hDCAz"), track.dcaz());
      histos.fill(HIST("tracks/hDCAxyVsDCAz"), track.dcaxy(), track.dcaz());
      histos.fill(HIST("tracks/hDCAxyVsPt"), track.pt(), track.dcaxy());
      histos.fill(HIST("tracks/hDCAzVsPt"), track.pt(), track.dcaz());

      // Tracks histos fill
      histos.fill(HIST("tracks/h1Eta"), track.eta());
      histos.fill(HIST("tracks/h1VarPhi"), track.phi());
      histos.fill(HIST("tracks/h2EtaVsPhi"), track.eta(), track.phi());
      histos.fill(HIST("tracks/h1pT"), track.pt());
      histos.fill(HIST("tracks/h1p"), track.p());

      //  TPC
      histos.fill(HIST("tracks/h2TPCsignVsTPCmomentum"), track.tpcInnerParam() / (1.f * track.sign()), track.tpcSignal());

      histos.fill(HIST("tracks/pion/h2PionVspTNSigmaTPC"), track.pt(), track.nsigTPCPi());
      histos.fill(HIST("tracks/kaon/h2KaonVspTNSigmaTPC"), track.pt(), track.nsigTPCKa());

      if (track.sign() > 0) {
        histos.fill(HIST("tracks/proton/h2ProtonTPCExpSignalDiffVsPt"), track.pt(), track.tpcExpSignalDiffPr());
        histos.fill(HIST("tracks/deuteron/h2DeuteronTPCExpSignalDiffVsPt"), track.pt(), track.tpcExpSignalDiffDe());
        histos.fill(HIST("tracks/proton/h2ProtonVspTNSigmaTPC"), track.pt(), track.nsigTPCPr());
        histos.fill(HIST("tracks/deuteron/h2DeuteronVspTNSigmaTPC"), track.pt(), track.nsigTPCD());
        histos.fill(HIST("tracks/helium/h2HeliumVspTNSigmaTPC"), track.pt(), track.nsigTPC3He());
      } else {
        histos.fill(HIST("tracks/proton/h2antiProtonTPCExpSignalDiffVsPt"), track.pt(), track.tpcExpSignalDiffPr());
        histos.fill(HIST("tracks/deuteron/h2antiDeuteronTPCExpSignalDiffVsPt"), track.pt(), track.tpcExpSignalDiffDe());
        histos.fill(HIST("tracks/proton/h2antiProtonVspTNSigmaTPC"), track.pt(), track.nsigTPCPr());
        histos.fill(HIST("tracks/deuteron/h2antiDeuteronVspTNSigmaTPC"), track.pt(), track.nsigTPCD());
        histos.fill(HIST("tracks/helium/h2antiHeliumVspTNSigmaTPC"), track.pt(), track.nsigTPC3He());
      }
      //  TOF
      histos.fill(HIST("tracks/pion/h2PionVspTNSigmaTOF"), track.pt(), track.nsigTOFPi());
      histos.fill(HIST("tracks/kaon/h2KaonVspTNSigmaTOF"), track.pt(), track.nsigTOFKa());
      if (track.sign() > 0) {
        histos.fill(HIST("tracks/proton/h2ProtonVspTNSigmaTOF"), track.pt(), track.nsigTOFPr());
        histos.fill(HIST("tracks/deuteron/h2DeuteronVspTNSigmaTOF"), track.pt(), track.nsigTOFD());
        histos.fill(HIST("tracks/helium/h2HeliumVspTNSigmaTOF"), track.pt(), track.nsigTOF3He());

        histos.fill(HIST("tracks/proton/h2ProtonTOFExpSignalDiffVsPt"), track.pt(), track.tofExpSignalDiffPr());
        histos.fill(HIST("tracks/deuteron/h2DeuteronTOFExpSignalDiffVsPt"), track.pt(), track.tofExpSignalDiffDe());

        if (enableEvTimeSplitting && track.hasTOF()) {
          if (track.isEvTimeTOF() && track.isEvTimeT0AC()) {
            histos.fill(HIST("tracks/evtime/ft0tof/proton/h2ProtonVspTNSigmaTOF"), track.pt(), track.nsigTOFPr());
            histos.fill(HIST("tracks/evtime/ft0tof/deuteron/h2DeuteronVspTNSigmaTOF"), track.pt(), track.nsigTOFD());
            histos.fill(HIST("tracks/evtime/ft0tof/proton/h2ProtonTOFExpSignalDiffVsPt"), track.pt(), track.tofExpSignalDiffPr());
            histos.fill(HIST("tracks/evtime/ft0tof/deuteron/h2DeuteronTOFExpSignalDiffVsPt"), track.pt(), track.tofExpSignalDiffDe());
            histos.fill(HIST("tracks/evtime/ft0tof/proton/h3ProtonNSigmaTPCvsNSigmaTOFvsPt"), track.nsigTPCPr(), track.nsigTOFPr(), track.pt());
            histos.fill(HIST("tracks/evtime/ft0tof/deuteron/h3DeuteronNSigmaTPCvsNSigmaTOFvsPt"), track.nsigTPCD(), track.nsigTOFD(), track.pt());
          } else if (track.isEvTimeT0AC()) {
            histos.fill(HIST("tracks/evtime/ft0/proton/h2ProtonVspTNSigmaTOF"), track.pt(), track.nsigTOFPr());
            histos.fill(HIST("tracks/evtime/ft0/deuteron/h2DeuteronVspTNSigmaTOF"), track.pt(), track.nsigTOFD());
            histos.fill(HIST("tracks/evtime/ft0/proton/h2ProtonTOFExpSignalDiffVsPt"), track.pt(), track.tofExpSignalDiffPr());
            histos.fill(HIST("tracks/evtime/ft0/deuteron/h2DeuteronTOFExpSignalDiffVsPt"), track.pt(), track.tofExpSignalDiffDe());
            histos.fill(HIST("tracks/evtime/ft0/proton/h3ProtonNSigmaTPCvsNSigmaTOFvsPt"), track.nsigTPCPr(), track.nsigTOFPr(), track.pt());
            histos.fill(HIST("tracks/evtime/ft0/deuteron/h3DeuteronNSigmaTPCvsNSigmaTOFvsPt"), track.nsigTPCD(), track.nsigTOFD(), track.pt());

          } else if (track.isEvTimeTOF()) {
            histos.fill(HIST("tracks/evtime/tof/proton/h2ProtonVspTNSigmaTOF"), track.pt(), track.nsigTOFPr());
            histos.fill(HIST("tracks/evtime/tof/deuteron/h2DeuteronVspTNSigmaTOF"), track.pt(), track.nsigTOFD());
            histos.fill(HIST("tracks/evtime/tof/proton/h2ProtonTOFExpSignalDiffVsPt"), track.pt(), track.tofExpSignalDiffPr());
            histos.fill(HIST("tracks/evtime/tof/deuteron/h2DeuteronTOFExpSignalDiffVsPt"), track.pt(), track.tofExpSignalDiffDe());
            histos.fill(HIST("tracks/evtime/tof/proton/h3ProtonNSigmaTPCvsNSigmaTOFvsPt"), track.nsigTPCPr(), track.nsigTOFPr(), track.pt());
            histos.fill(HIST("tracks/evtime/tof/deuteron/h3DeuteronNSigmaTPCvsNSigmaTOFvsPt"), track.nsigTPCD(), track.nsigTOFD(), track.pt());

          } else {
            histos.fill(HIST("tracks/evtime/fill/proton/h2ProtonVspTNSigmaTOF"), track.pt(), track.nsigTOFPr());
            histos.fill(HIST("tracks/evtime/fill/deuteron/h2DeuteronVspTNSigmaTOF"), track.pt(), track.nsigTOFD());
            histos.fill(HIST("tracks/evtime/fill/proton/h2ProtonTOFExpSignalDiffVsPt"), track.pt(), track.tofExpSignalDiffPr());
            histos.fill(HIST("tracks/evtime/fill/deuteron/h2DeuteronTOFExpSignalDiffVsPt"), track.pt(), track.tofExpSignalDiffDe());
            histos.fill(HIST("tracks/evtime/fill/proton/h3ProtonNSigmaTPCvsNSigmaTOFvsPt"), track.nsigTPCPr(), track.nsigTOFPr(), track.pt());
            histos.fill(HIST("tracks/evtime/fill/deuteron/h3DeuteronNSigmaTPCvsNSigmaTOFvsPt"), track.nsigTPCD(), track.nsigTOFD(), track.pt());
          }
        }

      } else {
        histos.fill(HIST("tracks/proton/h2antiProtonVspTNSigmaTOF"), track.pt(), track.nsigTOFPr());
        histos.fill(HIST("tracks/deuteron/h2antiDeuteronVspTNSigmaTOF"), track.pt(), track.nsigTOFD());
        histos.fill(HIST("tracks/helium/h2antiHeliumVspTNSigmaTOF"), track.pt(), track.nsigTOF3He());

        histos.fill(HIST("tracks/proton/h2antiProtonTOFExpSignalDiffVsPt"), track.pt(), track.tofExpSignalDiffPr());
        histos.fill(HIST("tracks/deuteron/h2antiDeuteronTOFExpSignalDiffVsPt"), track.pt(), track.tofExpSignalDiffDe());

        if (enableEvTimeSplitting && track.hasTOF()) {
          if (track.isEvTimeTOF() && track.isEvTimeT0AC()) {
            histos.fill(HIST("tracks/evtime/ft0tof/proton/h2antiProtonVspTNSigmaTOF"), track.pt(), track.nsigTOFPr());
            histos.fill(HIST("tracks/evtime/ft0tof/deuteron/h2antiDeuteronVspTNSigmaTOF"), track.pt(), track.nsigTOFD());
            histos.fill(HIST("tracks/evtime/ft0tof/proton/h2antiProtonTOFExpSignalDiffVsPt"), track.pt(), track.tofExpSignalDiffPr());
            histos.fill(HIST("tracks/evtime/ft0tof/deuteron/h2antiDeuteronTOFExpSignalDiffVsPt"), track.pt(), track.tofExpSignalDiffDe());
            histos.fill(HIST("tracks/evtime/ft0tof/proton/h3antiProtonNSigmaTPCvsNSigmaTOFvsPt"), track.nsigTPCPr(), track.nsigTOFPr(), track.pt());
            histos.fill(HIST("tracks/evtime/ft0tof/deuteron/h3antiDeuteronNSigmaTPCvsNSigmaTOFvsPt"), track.nsigTPCD(), track.nsigTOFD(), track.pt());
          } else if (track.isEvTimeT0AC()) {
            histos.fill(HIST("tracks/evtime/ft0/proton/h2antiProtonVspTNSigmaTOF"), track.pt(), track.nsigTOFPr());
            histos.fill(HIST("tracks/evtime/ft0/deuteron/h2antiDeuteronVspTNSigmaTOF"), track.pt(), track.nsigTOFD());
            histos.fill(HIST("tracks/evtime/ft0/proton/h2antiProtonTOFExpSignalDiffVsPt"), track.pt(), track.tofExpSignalDiffPr());
            histos.fill(HIST("tracks/evtime/ft0/deuteron/h2antiDeuteronTOFExpSignalDiffVsPt"), track.pt(), track.tofExpSignalDiffDe());
            histos.fill(HIST("tracks/evtime/ft0/proton/h3antiProtonNSigmaTPCvsNSigmaTOFvsPt"), track.nsigTPCPr(), track.nsigTOFPr(), track.pt());
            histos.fill(HIST("tracks/evtime/ft0/deuteron/h3antiDeuteronNSigmaTPCvsNSigmaTOFvsPt"), track.nsigTPCD(), track.nsigTOFD(), track.pt());
          } else if (track.isEvTimeTOF()) {
            histos.fill(HIST("tracks/evtime/tof/proton/h2antiProtonVspTNSigmaTOF"), track.pt(), track.nsigTOFPr());
            histos.fill(HIST("tracks/evtime/tof/deuteron/h2antiDeuteronVspTNSigmaTOF"), track.pt(), track.nsigTOFD());
            histos.fill(HIST("tracks/evtime/tof/proton/h2antiProtonTOFExpSignalDiffVsPt"), track.pt(), track.tofExpSignalDiffPr());
            histos.fill(HIST("tracks/evtime/tof/deuteron/h2antiDeuteronTOFExpSignalDiffVsPt"), track.pt(), track.tofExpSignalDiffDe());
            histos.fill(HIST("tracks/evtime/tof/proton/h3antiProtonNSigmaTPCvsNSigmaTOFvsPt"), track.nsigTPCPr(), track.nsigTOFPr(), track.pt());
            histos.fill(HIST("tracks/evtime/tof/deuteron/h3antiDeuteronNSigmaTPCvsNSigmaTOFvsPt"), track.nsigTPCD(), track.nsigTOFD(), track.pt());
          } else {
            histos.fill(HIST("tracks/evtime/fill/proton/h2antiProtonVspTNSigmaTOF"), track.pt(), track.nsigTOFPr());
            histos.fill(HIST("tracks/evtime/fill/deuteron/h2antiDeuteronVspTNSigmaTOF"), track.pt(), track.nsigTOFD());
            histos.fill(HIST("tracks/evtime/fill/proton/h2antiProtonTOFExpSignalDiffVsPt"), track.pt(), track.tofExpSignalDiffPr());
            histos.fill(HIST("tracks/evtime/fill/deuteron/h2antiDeuteronTOFExpSignalDiffVsPt"), track.pt(), track.tofExpSignalDiffDe());
            histos.fill(HIST("tracks/evtime/fill/proton/h3antiProtonNSigmaTPCvsNSigmaTOFvsPt"), track.nsigTPCPr(), track.nsigTOFPr(), track.pt());
            histos.fill(HIST("tracks/evtime/fill/deuteron/h3antiDeuteronNSigmaTPCvsNSigmaTOFvsPt"), track.nsigTPCD(), track.nsigTOFD(), track.pt());
          }
        }
      }

      // PID
      if (std::abs(track.nsigTPCPr()) < nsigmaTPCcut) {
        if (track.sign() > 0) {
          histos.fill(HIST("tracks/proton/h1ProtonSpectra"), track.pt());
          histos.fill(HIST("tracks/proton/dca/hDCAxyVsPtProton"), track.pt(), track.dcaxy());
          histos.fill(HIST("tracks/proton/dca/hDCAzVsPtProton"), track.pt(), track.dcaz());
        } else {
          histos.fill(HIST("tracks/proton/h1antiProtonSpectra"), track.pt());
          histos.fill(HIST("tracks/proton/dca/hDCAxyVsPtantiProton"), track.pt(), track.dcaxy());
          histos.fill(HIST("tracks/proton/dca/hDCAzVsPtantiProton"), track.pt(), track.dcaz());
        }
      }

      if (std::abs(track.nsigTPCD()) < nsigmaTPCcut) {
        if (track.sign() > 0) {
          histos.fill(HIST("tracks/deuteron/h1DeuteronSpectra"), track.pt());
          histos.fill(HIST("tracks/deuteron/dca/hDCAxyVsPtDeuteron"), track.pt(), track.dcaxy());
          histos.fill(HIST("tracks/deuteron/dca/hDCAzVsPtDeuteron"), track.pt(), track.dcaz());
        } else {
          histos.fill(HIST("tracks/deuteron/h1antiDeuteronSpectra"), track.pt());
          histos.fill(HIST("tracks/deuteron/dca/hDCAxyVsPtantiDeuteron"), track.pt(), track.dcaxy());
          histos.fill(HIST("tracks/deuteron/dca/hDCAzVsPtantiDeuteron"), track.pt(), track.dcaz());
        }
      }

      if (std::abs(track.nsigTPC3He()) < nsigmaTPCcut) {
        if (track.sign() > 0) {
          histos.fill(HIST("tracks/helium/h1HeliumSpectra"), track.pt());
          histos.fill(HIST("tracks/helium/dca/hDCAxyVsPtHelium"), track.pt(), track.dcaxy());
          histos.fill(HIST("tracks/helium/dca/hDCAzVsPtHelium"), track.pt(), track.dcaz());
        } else {
          histos.fill(HIST("tracks/helium/h1antiHeliumSpectra"), track.pt());
          histos.fill(HIST("tracks/helium/dca/hDCAxyVsPtantiHelium"), track.pt(), track.dcaxy());
          histos.fill(HIST("tracks/helium/dca/hDCAzVsPtantiHelium"), track.pt(), track.dcaz());
        }
      }

      if (track.hasTOF()) {
        histos.fill(HIST("tracks/h2TOFbetaVsP"), track.p() / (1.f * track.sign()), track.beta());

        if (track.sign() > 0) {
          histos.fill(HIST("tracks/proton/h3ProtonNSigmaTPCvsNSigmaTOFvsPt"), track.nsigTPCPr(), track.nsigTOFPr(), track.pt());
          histos.fill(HIST("tracks/deuteron/h3DeuteronNSigmaTPCvsNSigmaTOFvsPt"), track.nsigTPCD(), track.nsigTOFD(), track.pt());
        } else {
          histos.fill(HIST("tracks/proton/h3antiProtonNSigmaTPCvsNSigmaTOFvsPt"), track.nsigTPCPr(), track.nsigTOFPr(), track.pt());
          histos.fill(HIST("tracks/deuteron/h3antiDeuteronNSigmaTPCvsNSigmaTOFvsPt"), track.nsigTPCD(), track.nsigTOFD(), track.pt());
        }

        if (enableEvTimeSplitting) {
          if (track.isEvTimeTOF() && track.isEvTimeT0AC()) {
            histos.fill(HIST("tracks/evtime/ft0tof/h2TOFbetaVsP"), track.p() / (1.f * track.sign()), track.beta());
            histos.fill(HIST("tracks/evtime/ft0tof/h2TPCsignVsTPCmomentum"), track.tpcInnerParam() / (1.f * track.sign()), track.tpcSignal());
          } else if (track.isEvTimeT0AC()) {
            histos.fill(HIST("tracks/evtime/ft0/h2TOFbetaVsP"), track.p() / (1.f * track.sign()), track.beta());
            histos.fill(HIST("tracks/evtime/ft0/h2TPCsignVsTPCmomentum"), track.tpcInnerParam() / (1.f * track.sign()), track.tpcSignal());
          } else if (track.isEvTimeTOF()) {
            histos.fill(HIST("tracks/evtime/tof/h2TOFbetaVsP"), track.p() / (1.f * track.sign()), track.beta());
            histos.fill(HIST("tracks/evtime/tof/h2TPCsignVsTPCmomentum"), track.tpcInnerParam() / (1.f * track.sign()), track.tpcSignal());
          } else {
            histos.fill(HIST("tracks/evtime/fill/h2TOFbetaVsP"), track.p() / (1.f * track.sign()), track.beta());
            histos.fill(HIST("tracks/evtime/fill/h2TPCsignVsTPCmomentum"), track.tpcInnerParam() / (1.f * track.sign()), track.tpcSignal());
          }
        }

        if ((track.beta() * track.beta()) < 1.) {
          gamma = 1.f / TMath::Sqrt(1.f - (track.beta() * track.beta()));
          massTOF = track.p() / TMath::Sqrt(gamma * gamma - 1.f);
        } else {
          massTOF = -99.f;
        }

        histos.fill(HIST("tracks/h2TOFmassVsPt"), massTOF, track.pt());

        if (enableEvTimeSplitting) {
          if (track.isEvTimeTOF() && track.isEvTimeT0AC()) {
            histos.fill(HIST("tracks/evtime/ft0tof/h2TOFmassVsPt"), massTOF, track.pt());
          } else if (track.isEvTimeT0AC()) {
            histos.fill(HIST("tracks/evtime/ft0/h2TOFmassVsPt"), massTOF, track.pt());
          } else if (track.isEvTimeTOF()) {
            histos.fill(HIST("tracks/evtime/tof/h2TOFmassVsPt"), massTOF, track.pt());
          } else {
            histos.fill(HIST("tracks/evtime/fill/h2TOFmassVsPt"), massTOF, track.pt());
          }
        }

        if (std::abs(track.nsigTPCPr()) < nsigmaTPCcut) {
          if (track.sign() > 0) {
            histos.fill(HIST("tracks/proton/h2TOFmass2ProtonVsPt"), massTOF * massTOF - fMassProton * fMassProton, track.pt());
            // histos.fill(HIST("tracks/proton/h2ProtonTOFExpSignalDiffVsPtCut"), track.pt(), track.tofExpSignalDiffPr());

            if (enableEvTimeSplitting) {
              if (track.isEvTimeTOF() && track.isEvTimeT0AC()) {
                histos.fill(HIST("tracks/evtime/ft0tof/proton/h2TOFmass2ProtonVsPt"), massTOF * massTOF - fMassProton * fMassProton, track.pt());
                // histos.fill(HIST("tracks/evtime/ft0tof/proton/h2ProtonTOFExpSignalDiffVsPtCut"), track.pt(), track.tofExpSignalDiffPr());
              } else if (track.isEvTimeT0AC()) {
                histos.fill(HIST("tracks/evtime/ft0/proton/h2TOFmass2ProtonVsPt"), massTOF * massTOF - fMassProton * fMassProton, track.pt());
                // histos.fill(HIST("tracks/evtime/ft0/proton/h2ProtonTOFExpSignalDiffVsPtCut"), track.pt(), track.tofExpSignalDiffPr());
              } else if (track.isEvTimeTOF()) {
                histos.fill(HIST("tracks/evtime/tof/proton/h2TOFmass2ProtonVsPt"), massTOF * massTOF - fMassProton * fMassProton, track.pt());
                // histos.fill(HIST("tracks/evtime/tof/proton/h2ProtonTOFExpSignalDiffVsPtCut"), track.pt(), track.tofExpSignalDiffPr());
              } else {
                histos.fill(HIST("tracks/evtime/fill/proton/h2TOFmass2ProtonVsPt"), massTOF * massTOF - fMassProton * fMassProton, track.pt());
                // histos.fill(HIST("tracks/evtime/fill/proton/h2ProtonTOFExpSignalDiffVsPtCut"), track.pt(), track.tofExpSignalDiffPr());
              }
            }
          } else {
            histos.fill(HIST("tracks/proton/h2TOFmass2antiProtonVsPt"), massTOF * massTOF - fMassProton * fMassProton, track.pt());
            histos.fill(HIST("tracks/proton/h2antiProtonTOFExpSignalDiffVsPtCut"), track.pt(), track.tofExpSignalDiffPr());

            if (enableEvTimeSplitting) {
              if (track.isEvTimeTOF() && track.isEvTimeT0AC()) {
                histos.fill(HIST("tracks/evtime/ft0tof/proton/h2TOFmass2antiProtonVsPt"), massTOF * massTOF - fMassProton * fMassProton, track.pt());
                // histos.fill(HIST("tracks/evtime/ft0tof/proton/h2antiProtonTOFExpSignalDiffVsPtCut"), track.pt(), track.tofExpSignalDiffPr());
              } else if (track.isEvTimeT0AC()) {
                histos.fill(HIST("tracks/evtime/ft0/proton/h2TOFmass2antiProtonVsPt"), massTOF * massTOF - fMassProton * fMassProton, track.pt());
                // histos.fill(HIST("tracks/evtime/ft0/proton/h2antiProtonTOFExpSignalDiffVsPtCut"), track.pt(), track.tofExpSignalDiffPr());
              } else if (track.isEvTimeTOF()) {
                histos.fill(HIST("tracks/evtime/tof/proton/h2TOFmass2antiProtonVsPt"), massTOF * massTOF - fMassProton * fMassProton, track.pt());
                // histos.fill(HIST("tracks/evtime/tof/proton/h2antiProtonTOFExpSignalDiffVsPtCut"), track.pt(), track.tofExpSignalDiffPr());
              } else {
                histos.fill(HIST("tracks/evtime/fill/proton/h2TOFmass2antiProtonVsPt"), massTOF * massTOF - fMassProton * fMassProton, track.pt());
                // histos.fill(HIST("tracks/evtime/fill/proton/h2antiProtonTOFExpSignalDiffVsPtCut"), track.pt(), track.tofExpSignalDiffPr());
              }
            }
          }
        }

        if (std::abs(track.nsigTPCD()) < nsigmaTPCcut) {
          if (track.sign() > 0) {
            histos.fill(HIST("tracks/deuteron/h2TOFmass2DeuteronVsPt"), massTOF * massTOF - fMassDeuteron * fMassDeuteron, track.pt());
            histos.fill(HIST("tracks/deuteron/h2DeuteronTOFExpSignalDiffVsPtCut"), track.pt(), track.tofExpSignalDiffDe());

            if (enableEvTimeSplitting) {
              if (track.isEvTimeTOF() && track.isEvTimeT0AC()) {
                histos.fill(HIST("tracks/evtime/ft0tof/deuteron/h2TOFmass2DeuteronVsPt"), massTOF * massTOF - fMassDeuteron * fMassDeuteron, track.pt());
                // histos.fill(HIST("tracks/evtime/ft0tof/deuteron/h2DeuteronTOFExpSignalDiffVsPtCut"), track.pt(), track.tofExpSignalDiffDe());
              } else if (track.isEvTimeT0AC()) {
                histos.fill(HIST("tracks/evtime/ft0/deuteron/h2TOFmass2DeuteronVsPt"), massTOF * massTOF - fMassDeuteron * fMassDeuteron, track.pt());
                // histos.fill(HIST("tracks/evtime/ft0/deuteron/h2DeuteronTOFExpSignalDiffVsPtCut"), track.pt(), track.tofExpSignalDiffDe());
              } else if (track.isEvTimeTOF()) {
                histos.fill(HIST("tracks/evtime/tof/deuteron/h2TOFmass2DeuteronVsPt"), massTOF * massTOF - fMassDeuteron * fMassDeuteron, track.pt());
                // histos.fill(HIST("tracks/evtime/tof/deuteron/h2DeuteronTOFExpSignalDiffVsPtCut"), track.pt(), track.tofExpSignalDiffDe());
              } else {
                histos.fill(HIST("tracks/evtime/fill/deuteron/h2TOFmass2DeuteronVsPt"), massTOF * massTOF - fMassDeuteron * fMassDeuteron, track.pt());
                // histos.fill(HIST("tracks/evtime/fill/deuteron/h2DeuteronTOFExpSignalDiffVsPtCut"), track.pt(), track.tofExpSignalDiffDe());
              }
            }

          } else {
            histos.fill(HIST("tracks/deuteron/h2TOFmass2antiDeuteronVsPt"), massTOF * massTOF - fMassDeuteron * fMassDeuteron, track.pt());
            histos.fill(HIST("tracks/deuteron/h2DeuteronTOFExpSignalDiffVsPtCut"), track.pt(), track.tofExpSignalDiffDe());

            if (enableEvTimeSplitting) {
              if (track.isEvTimeTOF() && track.isEvTimeT0AC()) {
                histos.fill(HIST("tracks/evtime/ft0tof/deuteron/h2TOFmass2antiDeuteronVsPt"), massTOF * massTOF - fMassDeuteron * fMassDeuteron, track.pt());
                // histos.fill(HIST("tracks/evtime/ft0tof/deuteron/h2antiDeuteronTOFExpSignalDiffVsPtCut"), track.pt(), track.tofExpSignalDiffDe());
              } else if (track.isEvTimeT0AC()) {
                histos.fill(HIST("tracks/evtime/ft0/deuteron/h2TOFmass2antiDeuteronVsPt"), massTOF * massTOF - fMassDeuteron * fMassDeuteron, track.pt());
                // histos.fill(HIST("tracks/evtime/ft0/deuteron/h2antiDeuteronTOFExpSignalDiffVsPtCut"), track.pt(), track.tofExpSignalDiffDe());
              } else if (track.isEvTimeTOF()) {
                histos.fill(HIST("tracks/evtime/tof/deuteron/h2TOFmass2antiDeuteronVsPt"), massTOF * massTOF - fMassDeuteron * fMassDeuteron, track.pt());
                // histos.fill(HIST("tracks/evtime/tof/deuteron/h2antiDeuteronTOFExpSignalDiffVsPtCut"), track.pt(), track.tofExpSignalDiffDe());
              } else {
                histos.fill(HIST("tracks/evtime/fill/deuteron/h2TOFmass2antiDeuteronVsPt"), massTOF * massTOF - fMassDeuteron * fMassDeuteron, track.pt());
                // histos.fill(HIST("tracks/evtime/fill/deuteron/h2antiDeuteronTOFExpSignalDiffVsPtCut"), track.pt(), track.tofExpSignalDiffDe());
              }
            }
          }
        }

        if (std::abs(track.nsigTPC3He()) < nsigmaTPCcut) {
          if (track.sign() > 0)
            histos.fill(HIST("tracks/helium/h2TOFmass2HeliumVsPt"), massTOF * massTOF - fMassHelium * fMassHelium, track.pt());
          else
            histos.fill(HIST("tracks/helium/h2TOFmass2antiHeliumVsPt"), massTOF * massTOF - fMassHelium * fMassHelium, track.pt());
        }
      }

      if constexpr (IsMC) {
        bool isPhysPrim = track.isPhysicalPrimary();
        bool isProdByGen = track.producedByGenerator();

        // PID
        if (track.pdgCode() == PDGProton) {
          histos.fill(HIST("tracks/proton/h1ProtonSpectraTrue"), track.pt());
          histos.fill(HIST("tracks/proton/dca/hDCAxyVsPtProtonTrue"), track.pt(), track.dcaxy());
          histos.fill(HIST("tracks/proton/dca/hDCAzVsPtProtonTrue"), track.pt(), track.dcaz());
          if (isPhysPrim) {
            histos.fill(HIST("tracks/proton/h1ProtonSpectraTruePrim"), track.pt());
            histos.fill(HIST("tracks/proton/dca/hDCAxyVsPtProtonTruePrim"), track.pt(), track.dcaxy());
            histos.fill(HIST("tracks/proton/dca/hDCAzVsPtProtonTruePrim"), track.pt(), track.dcaz());
          }
          if (!isPhysPrim && isProdByGen) {
            histos.fill(HIST("tracks/proton/h1ProtonSpectraTrueSec"), track.pt());
            histos.fill(HIST("tracks/proton/dca/hDCAxyVsPtProtonTrueSec"), track.pt(), track.dcaxy());
            histos.fill(HIST("tracks/proton/dca/hDCAzVsPtProtonTrueSec"), track.pt(), track.dcaz());
          }
          if (!isPhysPrim && !isProdByGen) {
            histos.fill(HIST("tracks/proton/h1ProtonSpectraTrueTransport"), track.pt());
            histos.fill(HIST("tracks/proton/dca/hDCAxyVsPtProtonTrueTransport"), track.pt(), track.dcaxy());
            histos.fill(HIST("tracks/proton/dca/hDCAzVsPtProtonTrueTransport"), track.pt(), track.dcaz());
          }
        }

        if (track.pdgCode() == -PDGProton) {
          histos.fill(HIST("tracks/proton/h1antiProtonSpectraTrue"), track.pt());
          histos.fill(HIST("tracks/proton/dca/hDCAxyVsPtantiProtonTrue"), track.pt(), track.dcaxy());
          histos.fill(HIST("tracks/proton/dca/hDCAzVsPtantiProtonTrue"), track.pt(), track.dcaz());
          if (isPhysPrim) {
            histos.fill(HIST("tracks/proton/h1antiProtonSpectraTruePrim"), track.pt());
            histos.fill(HIST("tracks/proton/dca/hDCAxyVsPtantiProtonTruePrim"), track.pt(), track.dcaxy());
            histos.fill(HIST("tracks/proton/dca/hDCAzVsPtantiProtonTruePrim"), track.pt(), track.dcaz());
          }
          if (!isPhysPrim && isProdByGen) {
            histos.fill(HIST("tracks/proton/h1antiProtonSpectraTrueSec"), track.pt());
            histos.fill(HIST("tracks/proton/dca/hDCAxyVsPtantiProtonTrueSec"), track.pt(), track.dcaxy());
            histos.fill(HIST("tracks/proton/dca/hDCAzVsPtantiProtonTrueSec"), track.pt(), track.dcaz());
          }
          if (!isPhysPrim && !isProdByGen) {
            histos.fill(HIST("tracks/proton/h1antiProtonSpectraTrueTransport"), track.pt());
            histos.fill(HIST("tracks/proton/dca/hDCAxyVsPtantiProtonTrueTransport"), track.pt(), track.dcaxy());
            histos.fill(HIST("tracks/proton/dca/hDCAzVsPtantiProtonTrueTransport"), track.pt(), track.dcaz());
          }
        }
        if (track.pdgCode() == PDGDeuteron) {
          histos.fill(HIST("tracks/deuteron/h1DeuteronSpectraTrue"), track.pt());
          histos.fill(HIST("tracks/deuteron/dca/hDCAxyVsPtDeuteronTrue"), track.pt(), track.dcaxy());
          histos.fill(HIST("tracks/deuteron/dca/hDCAzVsPtDeuteronTrue"), track.pt(), track.dcaz());
          if (isPhysPrim) {
            histos.fill(HIST("tracks/deuteron/h1DeuteronSpectraTruePrim"), track.pt());
            histos.fill(HIST("tracks/deuteron/dca/hDCAxyVsPtDeuteronTruePrim"), track.pt(), track.dcaxy());
            histos.fill(HIST("tracks/deuteron/dca/hDCAzVsPtDeuteronTruePrim"), track.pt(), track.dcaz());
          }
          if (!isPhysPrim && isProdByGen) {
            histos.fill(HIST("tracks/deuteron/h1DeuteronSpectraTrueSec"), track.pt());
            histos.fill(HIST("tracks/deuteron/dca/hDCAxyVsPtDeuteronTrueSec"), track.pt(), track.dcaxy());
            histos.fill(HIST("tracks/deuteron/dca/hDCAzVsPtDeuteronTrueSec"), track.pt(), track.dcaz());
          }
          if (!isPhysPrim && !isProdByGen) {
            histos.fill(HIST("tracks/deuteron/h1DeuteronSpectraTrueTransport"), track.pt());
            histos.fill(HIST("tracks/deuteron/dca/hDCAxyVsPtDeuteronTrueTransport"), track.pt(), track.dcaxy());
            histos.fill(HIST("tracks/deuteron/dca/hDCAzVsPtDeuteronTrueTransport"), track.pt(), track.dcaz());
          }
        }

        if (track.pdgCode() == -PDGDeuteron) {
          histos.fill(HIST("tracks/deuteron/h1antiDeuteronSpectraTrue"), track.pt());
          histos.fill(HIST("tracks/deuteron/dca/hDCAxyVsPtantiDeuteronTrue"), track.pt(), track.dcaxy());
          histos.fill(HIST("tracks/deuteron/dca/hDCAzVsPtantiDeuteronTrue"), track.pt(), track.dcaz());
          if (isPhysPrim) {
            histos.fill(HIST("tracks/deuteron/h1antiDeuteronSpectraTruePrim"), track.pt());
            histos.fill(HIST("tracks/deuteron/dca/hDCAxyVsPtantiDeuteronTruePrim"), track.pt(), track.dcaxy());
            histos.fill(HIST("tracks/deuteron/dca/hDCAzVsPtantiDeuteronTruePrim"), track.pt(), track.dcaz());
          }
          if (!isPhysPrim && isProdByGen) {
            histos.fill(HIST("tracks/deuteron/h1antiDeuteronSpectraTrueSec"), track.pt());
            histos.fill(HIST("tracks/deuteron/dca/hDCAxyVsPtantiDeuteronTrueSec"), track.pt(), track.dcaxy());
            histos.fill(HIST("tracks/deuteron/dca/hDCAzVsPtantiDeuteronTrueSec"), track.pt(), track.dcaz());
          }
          if (!isPhysPrim && !isProdByGen) {
            histos.fill(HIST("tracks/deuteron/h1antiDeuteronSpectraTrueTransport"), track.pt());
            histos.fill(HIST("tracks/deuteron/dca/hDCAxyVsPtantiDeuteronTrueTransport"), track.pt(), track.dcaxy());
            histos.fill(HIST("tracks/deuteron/dca/hDCAzVsPtantiDeuteronTrueTransport"), track.pt(), track.dcaz());
          }
        }

        if (track.pdgCode() == PDGHelium) {
          histos.fill(HIST("tracks/helium/h1HeliumSpectraTrue"), track.pt());
          histos.fill(HIST("tracks/helium/dca/hDCAxyVsPtHeliumTrue"), track.pt(), track.dcaxy());
          histos.fill(HIST("tracks/helium/dca/hDCAzVsPtHeliumTrue"), track.pt(), track.dcaz());
          if (isPhysPrim) {
            histos.fill(HIST("tracks/helium/h1HeliumSpectraTruePrim"), track.pt());
            histos.fill(HIST("tracks/helium/dca/hDCAxyVsPtHeliumTruePrim"), track.pt(), track.dcaxy());
            histos.fill(HIST("tracks/helium/dca/hDCAzVsPtHeliumTruePrim"), track.pt(), track.dcaz());
          }
          if (!isPhysPrim && isProdByGen) {
            histos.fill(HIST("tracks/helium/h1HeliumSpectraTrueSec"), track.pt());
            histos.fill(HIST("tracks/helium/dca/hDCAxyVsPtHeliumTrueSec"), track.pt(), track.dcaxy());
            histos.fill(HIST("tracks/helium/dca/hDCAzVsPtHeliumTrueSec"), track.pt(), track.dcaz());
          }
          if (!isPhysPrim && !isProdByGen) {
            histos.fill(HIST("tracks/helium/h1HeliumSpectraTrueTransport"), track.pt());
            histos.fill(HIST("tracks/helium/dca/hDCAxyVsPtHeliumTrueTransport"), track.pt(), track.dcaxy());
            histos.fill(HIST("tracks/helium/dca/hDCAzVsPtHeliumTrueTransport"), track.pt(), track.dcaz());
          }
        }

        if (track.pdgCode() == -PDGHelium) {
          histos.fill(HIST("tracks/helium/h1antiHeliumSpectraTrue"), track.pt());
          histos.fill(HIST("tracks/helium/dca/hDCAxyVsPtantiHeliumTrue"), track.pt(), track.dcaxy());
          histos.fill(HIST("tracks/helium/dca/hDCAzVsPtantiHeliumTrue"), track.pt(), track.dcaz());
          if (isPhysPrim) {
            histos.fill(HIST("tracks/helium/h1antiHeliumSpectraTruePrim"), track.pt());
            histos.fill(HIST("tracks/helium/dca/hDCAxyVsPtantiHeliumTruePrim"), track.pt(), track.dcaxy());
            histos.fill(HIST("tracks/helium/dca/hDCAzVsPtantiHeliumTruePrim"), track.pt(), track.dcaz());
          }
          if (!isPhysPrim && isProdByGen) {
            histos.fill(HIST("tracks/helium/h1antiHeliumSpectraTrueSec"), track.pt());
            histos.fill(HIST("tracks/helium/dca/hDCAxyVsPtantiHeliumTrueSec"), track.pt(), track.dcaxy());
            histos.fill(HIST("tracks/helium/dca/hDCAzVsPtantiHeliumTrueSec"), track.pt(), track.dcaz());
          }
          if (!isPhysPrim && !isProdByGen) {
            histos.fill(HIST("tracks/helium/h1antiHeliumSpectraTrueTransport"), track.pt());
            histos.fill(HIST("tracks/helium/dca/hDCAxyVsPtantiHeliumTrueTransport"), track.pt(), track.dcaxy());
            histos.fill(HIST("tracks/helium/dca/hDCAzVsPtantiHeliumTrueTransport"), track.pt(), track.dcaz());
          }
        }
      }
    }
  }

  void processData(o2::aod::LfCandNucleusFullEvents::iterator const& event,
                   o2::aod::LfCandNucleusFull const& tracks)
  {
    fillHistograms<false>(event, tracks);
  } // CLOSING PROCESS DATA
  PROCESS_SWITCH(LFNucleiBATask, processData, "process data", true);

  void processMCReco(o2::aod::LfCandNucleusFullEvents::iterator const& event,
                     soa::Join<o2::aod::LfCandNucleusFull, o2::aod::LfCandNucleusMC> const& tracks)
  {
    fillHistograms<true>(event, tracks);
  } // CLOSING PROCESS MC RECO
  PROCESS_SWITCH(LFNucleiBATask, processMCReco, "process mc reco", false);
  Int_t nCount = 0;

  // LOOP OVER GENERATED MC PARTICLES
  void processMCGen(aod::McCollision const& mcCollision, aod::McParticles_001& mcParticles)
  {
    nCount++;
    histos.fill(HIST("spectraGen/histGenVetxZ"), mcCollision.posZ());
    for (auto& mcParticleGen : mcParticles) {
      if (abs(mcParticleGen.y()) > std::abs(yCut)) {
        continue;
      }

      bool isPhysPrim = mcParticleGen.isPhysicalPrimary();
      bool isProdByGen = mcParticleGen.producedByGenerator();
      if (std::abs(mcParticleGen.pdgCode()) == PDGPion) {
        histos.fill(HIST("spectraGen/histGenPtPion"), mcParticleGen.pt());
        if (isPhysPrim)
          histos.fill(HIST("spectraGen/histGenPtPionPrim"), mcParticleGen.pt());
        if (!isPhysPrim && isProdByGen)
          histos.fill(HIST("spectraGen/histGenPtPionSec"), mcParticleGen.pt());
        if (!isPhysPrim && !isProdByGen)
          histos.fill(HIST("spectraGen/histSecTransportPtPion"), mcParticleGen.pt());
      }
      if (std::abs(mcParticleGen.pdgCode()) == PDGKaon) {
        histos.fill(HIST("spectraGen/histGenPtKaon"), mcParticleGen.pt());
        if (isPhysPrim)
          histos.fill(HIST("spectraGen/histGenPtKaonPrim"), mcParticleGen.pt());
        if (!isPhysPrim && isProdByGen)
          histos.fill(HIST("spectraGen/histGenPtKaonSec"), mcParticleGen.pt());
        if (!isPhysPrim && !isProdByGen)
          histos.fill(HIST("spectraGen/histSecTransportPtKaon"), mcParticleGen.pt());
      }
      if (mcParticleGen.pdgCode() == PDGProton) {
        histos.fill(HIST("spectraGen/histGenPtProton"), mcParticleGen.pt());
        if (isPhysPrim)
          histos.fill(HIST("spectraGen/histGenPtProtonPrim"), mcParticleGen.pt());
        if (!isPhysPrim && isProdByGen)
          histos.fill(HIST("spectraGen/histGenPtProtonSec"), mcParticleGen.pt());
        if (!isPhysPrim && !isProdByGen)
          histos.fill(HIST("spectraGen/histSecTransportPtProton"), mcParticleGen.pt());
      }
      if (mcParticleGen.pdgCode() == -PDGProton) {
        histos.fill(HIST("spectraGen/histGenPtantiProton"), mcParticleGen.pt());
        if (isPhysPrim)
          histos.fill(HIST("spectraGen/histGenPtantiProtonPrim"), mcParticleGen.pt());
        if (!isPhysPrim && isProdByGen)
          histos.fill(HIST("spectraGen/histGenPtantiProtonSec"), mcParticleGen.pt());
        if (!isPhysPrim && !isProdByGen)
          histos.fill(HIST("spectraGen/histSecTransportPtantiProton"), mcParticleGen.pt());
      }
      if (mcParticleGen.pdgCode() == PDGDeuteron) {
        histos.fill(HIST("spectraGen/histGenPtD"), mcParticleGen.pt());
        if (isPhysPrim)
          histos.fill(HIST("spectraGen/histGenPtDPrim"), mcParticleGen.pt());
        if (!isPhysPrim && isProdByGen)
          histos.fill(HIST("spectraGen/histGenPtDSec"), mcParticleGen.pt());
        if (!isPhysPrim && !isProdByGen)
          histos.fill(HIST("spectraGen/histSecTransportPtD"), mcParticleGen.pt());
      }
      if (mcParticleGen.pdgCode() == -PDGDeuteron) {
        histos.fill(HIST("spectraGen/histGenPtantiD"), mcParticleGen.pt());
        if (isPhysPrim)
          histos.fill(HIST("spectraGen/histGenPtantiDPrim"), mcParticleGen.pt());
        if (!isPhysPrim && isProdByGen)
          histos.fill(HIST("spectraGen/histGenPtantiDSec"), mcParticleGen.pt());
        if (!isPhysPrim && !isProdByGen)
          histos.fill(HIST("spectraGen/histSecTransportPtantiD"), mcParticleGen.pt());
      }
      if (mcParticleGen.pdgCode() == PDGHelium) {
        histos.fill(HIST("spectraGen/histGenPtHe"), mcParticleGen.pt());
        if (isPhysPrim)
          histos.fill(HIST("spectraGen/histGenPtHePrim"), mcParticleGen.pt());
        if (!isPhysPrim && isProdByGen)
          histos.fill(HIST("spectraGen/histGenPtHeSec"), mcParticleGen.pt());
        if (!isPhysPrim && !isProdByGen)
          histos.fill(HIST("spectraGen/histSecTransportPtHe"), mcParticleGen.pt());
      }
      if (mcParticleGen.pdgCode() == -PDGHelium) {
        histos.fill(HIST("spectraGen/histGenPtantiHe"), mcParticleGen.pt());
        if (isPhysPrim)
          histos.fill(HIST("spectraGen/histGenPtantiHePrim"), mcParticleGen.pt());
        if (!isPhysPrim && isProdByGen)
          histos.fill(HIST("spectraGen/histGenPtantiHeSec"), mcParticleGen.pt());
        if (!isPhysPrim && !isProdByGen)
          histos.fill(HIST("spectraGen/histSecTransportPtantiHe"), mcParticleGen.pt());
      }
    }
  } // Close processMCGen
  PROCESS_SWITCH(LFNucleiBATask, processMCGen, "process MC Generated", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<LFNucleiBATask>(cfgc)};
}
