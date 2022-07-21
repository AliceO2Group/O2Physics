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
/// \brief  Analysis task for the measurement of the coalescence parameter
///
/// \author Rutuparna Rath <rutuparna.rath@cern.ch> and Giovanni Malfattore <giovanni.malfattore@cern.ch>
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
    histos.add<TH1>("tracks/h1Eta", "pseudoRapidity; #eta; counts", HistType::kTH1F, {{200, -1.0, 1.0}});
    histos.add<TH1>("tracks/h1VarPhi", "#phi; #phi; counts", HistType::kTH1F, {{63, 0.0, 6.3}});
    histos.add<TH2>("tracks/h2EtaVsPhi", "#eta vs #phi; #eta; #phi", HistType::kTH2F, {{200, -1.0, 1.0}, {63, 0.0, 6.3}});
    histos.add<TH1>("tracks/h1pT", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{1000, 0., 10}});
    histos.add<TH1>("tracks/h1p", "Track momentum; p (GeV/#it{c}); counts", HistType::kTH1F, {{1000, 0., 10.}});

    // tracks
    histos.add<TH1>("tracks/hDCAxy", "DCAxy; #DCAxy; counts", HistType::kTH1F, {{200, -2.0, 2.0}});
    histos.add<TH1>("tracks/hDCAz", "DCAz; #DCAz; counts", HistType::kTH1F, {{200, -2.0, 2.0}});
    histos.add<TH2>("tracks/hDCAxyVsDCAz", "DCAxy vs DCAz; DCAxy (cm); DCAz (cm)", HistType::kTH2F, {{200, -2.0, 2.0}, {200, -2.0, 2.0}});

    histos.add<TH1>("tracks/proton/h1ProtonSpectra", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
    histos.add<TH1>("tracks/deuteron/h1DeuteronSpectra", "#it{p}_{T} (d)", HistType::kTH1F, {ptAxis});
    histos.add<TH1>("tracks/helium/h1HeliumSpectra", "#it{p}_{T} (He)", HistType::kTH1F, {ptAxis});

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
    }

    histos.add<TH2>("tracks/h2TPCsignVsTPCmomentum", "-dE/dX vs p; p (GeV/c); -dE/dx (a.u.)", HistType::kTH2F, {{500, 0.0, 5.0}, {81000, 0.0, 1E3}});
    histos.add<TH2>("tracks/h2TOFbetaVsP", "#beta (TOF) vs p; p (GeV/c); #beta", HistType::kTH2F, {{500, 0.0, 5.0}, {1200, 0.0, 1.2}});

    histos.add<TH2>("tracks/pion/h2PionVspTNSigmaTPC", "NSigmaTPC(pi) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptAxis}, {200, -10, 10.}});
    histos.add<TH2>("tracks/kaon/h2KaonVspTNSigmaTPC", "NSigmaTPC(Ka) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptAxis}, {200, -10, 10.}});
    histos.add<TH2>("tracks/proton/h2ProtonVspTNSigmaTPC", "NSigmaTPC(proton) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptAxis}, {200, -10, 10.}});
    histos.add<TH2>("tracks/deuteron/h2DeuteronVspTNSigmaTPC", "NSigmaTPC(D) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptAxis}, {200, -10, 10.}});
    histos.add<TH2>("tracks/helium/h2HeliumVspTNSigmaTPC", "NSigmaTPC(He) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptAxis}, {200, -10, 10.}});
    histos.add<TH2>("tracks/pion/h2PionVspTNSigmaTOF", "NSigmaTOF(pi) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {200, -10, 10.}});
    histos.add<TH2>("tracks/kaon/h2KaonVspTNSigmaTOF", "NSigmaTOF(Ka) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {200, -10, 10.}});
    histos.add<TH2>("tracks/proton/h2ProtonVspTNSigmaTOF", "NSigmaTOF(proton) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {200, -10, 10.}});
    histos.add<TH2>("tracks/deuteron/h2DeuteronVspTNSigmaTOF", "NSigmaTOF(D) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {200, -10, 10.}});
    histos.add<TH2>("tracks/helium/h2HeliumVspTNSigmaTOF", "NSigmaTOF(He) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {200, -10, 10.}});

    // TOF mass histograms
    histos.add<TH2>("tracks/h2TOFmassVsPt", "h2TOFmassVsPt; TOFmass; #it{p}_{T} (GeV)", HistType::kTH2F, {{600, 0., 3.}, {500, 0., 5.}});

    // TOF mass squared histograms
    histos.add<TH2>("tracks/proton/h2TOFmass2ProtonVsPt", "#Delta M^{2} (p) vs #it{p}_{T}; #Delta M^{2} (p); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{600, -3., 3.}, {800, 0., 8.}});
    histos.add<TH2>("tracks/deuteron/h2TOFmass2DeuteronVsPt", "#Delta M^{2} (d) vs #it{p}_{T}; #Delta M^{2} (d); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{1000, -5., 5.}, {800, 0., 8.}});
    histos.add<TH2>("tracks/helium/h2TOFmass2HeliumVsPt", "#Delta M^{2} (He) vs #it{p}_{T}; #Delta M^{2} (He); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{1800, -9., 9.}, {800, 0., 8.}});

    // MC histograms
    // AxisSpec ptAxis = {2000, 0.f, 20.f, "#it{p}_{T} (GeV/#it{c})"};
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

      // Tracks histos fill
      histos.fill(HIST("tracks/h1Eta"), track.eta());
      histos.fill(HIST("tracks/h1VarPhi"), track.phi());
      histos.fill(HIST("tracks/h2EtaVsPhi"), track.eta(), track.phi());
      histos.fill(HIST("tracks/h1pT"), track.pt());
      histos.fill(HIST("tracks/h1p"), track.p());

      //  TPC
      // if(std::abs(track.nsigTPCD()) < 5. && std::abs(track.nsigTOFD()) < 3.)
      histos.fill(HIST("tracks/h2TPCsignVsTPCmomentum"), track.tpcInnerParam(), track.tpcSignal());

      histos.fill(HIST("tracks/pion/h2PionVspTNSigmaTPC"), track.pt(), track.nsigTPCPi());
      histos.fill(HIST("tracks/kaon/h2KaonVspTNSigmaTPC"), track.pt(), track.nsigTPCKa());
      histos.fill(HIST("tracks/proton/h2ProtonVspTNSigmaTPC"), track.pt(), track.nsigTPCPr());
      histos.fill(HIST("tracks/deuteron/h2DeuteronVspTNSigmaTPC"), track.pt(), track.nsigTPCD());
      histos.fill(HIST("tracks/helium/h2HeliumVspTNSigmaTPC"), track.pt(), track.nsigTPC3He());

      //  TOF
      histos.fill(HIST("tracks/pion/h2PionVspTNSigmaTOF"), track.pt(), track.nsigTOFPi());
      histos.fill(HIST("tracks/kaon/h2KaonVspTNSigmaTOF"), track.pt(), track.nsigTOFKa());
      histos.fill(HIST("tracks/proton/h2ProtonVspTNSigmaTOF"), track.pt(), track.nsigTOFPr());
      histos.fill(HIST("tracks/deuteron/h2DeuteronVspTNSigmaTOF"), track.pt(), track.nsigTOFD());
      histos.fill(HIST("tracks/helium/h2HeliumVspTNSigmaTOF"), track.pt(), track.nsigTOF3He());

      // PID
      if (std::abs(track.nsigTPCPr()) < nsigmaTPCcut)
        histos.fill(HIST("tracks/proton/h1ProtonSpectra"), track.pt());

      if (std::abs(track.nsigTPCD()) < nsigmaTPCcut)
        histos.fill(HIST("tracks/deuteron/h1DeuteronSpectra"), track.pt());

      if (std::abs(track.nsigTPC3He()) < nsigmaTPCcut)
        histos.fill(HIST("tracks/helium/h1HeliumSpectra"), track.pt());

      if (track.hasTOF()) {
        histos.fill(HIST("tracks/h2TOFbetaVsP"), track.p(), track.beta());
        if ((track.beta() * track.beta()) < 1.) {
          gamma = 1.f / TMath::Sqrt(1.f - (track.beta() * track.beta()));
          massTOF = track.p() / TMath::Sqrt(gamma * gamma - 1.f);
        } else {
          massTOF = -99.f;
        }
        histos.fill(HIST("tracks/h2TOFmassVsPt"), massTOF, track.pt());

        if (std::abs(track.nsigTPCPr()) < nsigmaTPCcut)
          histos.fill(HIST("tracks/proton/h2TOFmass2ProtonVsPt"), massTOF * massTOF - fMassProton * fMassProton, track.pt());

        if (std::abs(track.nsigTPCD()) < nsigmaTPCcut)
          histos.fill(HIST("tracks/deuteron/h2TOFmass2DeuteronVsPt"), massTOF * massTOF - fMassDeuteron * fMassDeuteron, track.pt());

        if (std::abs(track.nsigTPC3He()) < nsigmaTPCcut)
          histos.fill(HIST("tracks/helium/h2TOFmass2HeliumVsPt"), massTOF * massTOF - fMassHelium * fMassHelium, track.pt());
      }
      if constexpr (IsMC) {
        bool isPhysPrim = track.isPhysicalPrimary();
        bool isProdByGen = track.producedByGenerator();

        // PID
        if (track.pdgCode() == PDGProton) {
          histos.fill(HIST("tracks/proton/h1ProtonSpectraTrue"), track.pt());
          if (isPhysPrim)
            histos.fill(HIST("tracks/proton/h1ProtonSpectraTruePrim"), track.pt());
          if (!isPhysPrim && isProdByGen)
            histos.fill(HIST("tracks/proton/h1ProtonSpectraTrueSec"), track.pt());
          if (!isPhysPrim && !isProdByGen)
            histos.fill(HIST("tracks/proton/h1ProtonSpectraTrueTransport"), track.pt());
        }

        if (track.pdgCode() == -PDGProton) {
          histos.fill(HIST("tracks/proton/h1antiProtonSpectraTrue"), track.pt());
          if (isPhysPrim)
            histos.fill(HIST("tracks/proton/h1antiProtonSpectraTruePrim"), track.pt());
          if (!isPhysPrim && isProdByGen)
            histos.fill(HIST("tracks/proton/h1antiProtonSpectraTrueSec"), track.pt());
          if (!isPhysPrim && !isProdByGen)
            histos.fill(HIST("tracks/proton/h1antiProtonSpectraTrueTransport"), track.pt());
        }
        if (track.pdgCode() == PDGDeuteron) {
          histos.fill(HIST("tracks/deuteron/h1DeuteronSpectraTrue"), track.pt());
          if (isPhysPrim)
            histos.fill(HIST("tracks/deuteron/h1DeuteronSpectraTruePrim"), track.pt());
          if (!isPhysPrim && isProdByGen)
            histos.fill(HIST("tracks/deuteron/h1DeuteronSpectraTrueSec"), track.pt());
          if (!isPhysPrim && !isProdByGen)
            histos.fill(HIST("tracks/deuteron/h1DeuteronSpectraTrueTransport"), track.pt());
        }

        if (track.pdgCode() == -PDGDeuteron) {
          histos.fill(HIST("tracks/deuteron/h1antiDeuteronSpectraTrue"), track.pt());
          if (isPhysPrim)
            histos.fill(HIST("tracks/deuteron/h1antiDeuteronSpectraTruePrim"), track.pt());
          if (!isPhysPrim && isProdByGen)
            histos.fill(HIST("tracks/deuteron/h1antiDeuteronSpectraTrueSec"), track.pt());
          if (!isPhysPrim && !isProdByGen)
            histos.fill(HIST("tracks/deuteron/h1antiDeuteronSpectraTrueTransport"), track.pt());
        }

        if (track.pdgCode() == PDGHelium) {
          histos.fill(HIST("tracks/helium/h1HeliumSpectraTrue"), track.pt());
          if (isPhysPrim)
            histos.fill(HIST("tracks/helium/h1HeliumSpectraTruePrim"), track.pt());
          if (!isPhysPrim && isProdByGen)
            histos.fill(HIST("tracks/helium/h1HeliumSpectraTrueSec"), track.pt());
          if (!isPhysPrim && !isProdByGen)
            histos.fill(HIST("tracks/helium/h1HeliumSpectraTrueTransport"), track.pt());
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
