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

#include "Common/Core/PID/PIDResponse.h"
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
  Configurable<float> cfgCutVertex{"cfgCutVSertex", 10.0f, "Accepted z-vertex range"};

  void init(o2::framework::InitContext&)
  {
    histos.add<TH1>("event/h1VtxZ", "V_{z};V_{z} (in cm); counts", HistType::kTH1F, {{3000, -15, 15}});
    histos.add<TH1>("event/h1CentV0M", "V0M; Multiplicity; counts", HistType::kTH1F, {{27000, 0, 27000}});

    histos.add<TH1>("qa/h1TPCncr", "number of crossed rows in TPC; TPCncr; counts", HistType::kTH1F, {{150, 60, 170}});
    histos.add<TH1>("qa/h1rTPC", "ratio of ncr over findable in TPC; rTPC; counts", HistType::kTH1F, {{200, 0.9, 1.8}});
    histos.add<TH1>("qa/h1chi2ITS", "#chi^{2}_{ITS}/n_{ITS}; #chi^{2}_{ITS}/n_{ITS};counts", HistType::kTH1F, {{51, -0.5, 50.5}});
    histos.add<TH1>("qa/h1chi2TPC", "#chi^{2}_{TPC}/n_{TPC}; #chi^{2}_{TPC}/n_{TPC}; counts", HistType::kTH1F, {{11, -0.5, 10.5}});

    // trackQA
    histos.add<TH1>("tracks/h1Eta", "pseudoRapidity; #eta; counts", HistType::kTH1F, {{200, -1.0, 1.0}});
    histos.add<TH1>("tracks/h1VarPhi", "#phi; #phi; counts", HistType::kTH1F, {{80, -0.5, 7.5}});
    histos.add<TH2>("tracks/h2EtaVsPhi", "#eta (TOF) vs #phi; #eta; #phi", HistType::kTH2F, {{200, -1.0, 1.0}, {80, -0.5, 7.5}});
    histos.add<TH1>("tracks/h1pT", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{1000, 0., 10}});
    histos.add<TH1>("tracks/h1p", "Track momentum; p (GeV/#it{c}); counts", HistType::kTH1F, {{1000, 0., 10.}});

    // tracks
    histos.add<TH1>("tracks/proton/h1ProtonSpectra", "#it{p}_{T} (p); #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{800, 0, 8}});
    histos.add<TH1>("tracks/deuteron/h1DeuteronSpectra", "#it{p}_{T} (d); #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{800, 0, 8}});
    histos.add<TH1>("tracks/helium/h1HeliumSpectra", "#it{p}_{T} (He); #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{800, 0, 8}});

    histos.add<TH2>("tracks/h2TPCsignVsTPCmomentum", "-dE/dX vs p; p (GeV/c); -dE/dx (a.u.)", HistType::kTH2F, {{500, 0.0, 5.0}, {81000, 0.0, 1E3}});
    histos.add<TH2>("tracks/h2TOFbetaVsP", "#beta (TOF) vs p; p (GeV/c); #beta", HistType::kTH2F, {{500, 0.0, 5.0}, {1200, 0.0, 1.2}});

    histos.add<TH2>("tracks/pion/h2PionVspNSigmaTPC", "NSigmaTPC(pi) vs p; Momentum(p) in GeV; NSigmaTPC", HistType::kTH2F, {{100, 0, 10.}, {200, -10, 10.}});
    histos.add<TH2>("tracks/kaon/h2KaonVspNSigmaTPC", "NSigmaTPC(Ka) vs p; Momentum(p) in GeV; NSigmaTPC", HistType::kTH2F, {{100, 0, 10.}, {200, -10, 10.}});
    histos.add<TH2>("tracks/proton/h2ProtonVspNSigmaTPC", "NSigmaTPC(proton) vs p; Momentum(p) in GeV; NSigmaTPC", HistType::kTH2F, {{100, 0, 10.}, {200, -10, 10.}});
    histos.add<TH2>("tracks/deuteron/h2DeuteronVspNSigmaTPC", "NSigmaTPC(D) vs p; Momentum(p) in GeV; NSigmaTPC", HistType::kTH2F, {{100, 0, 10.}, {200, -10, 10.}});
    histos.add<TH2>("tracks/helium/h2HeliumVspNSigmaTPC", "NSigmaTPC(He) vs p; Momentum(p) in GeV; NSigmaTPC", HistType::kTH2F, {{100, 0, 10.}, {200, -10, 10.}});
    histos.add<TH2>("tracks/pion/h2PionVspNSigmaTOF", "NSigmaTOF(pi) vs p; Momentum(p) in GeV; NSigmaTOF", HistType::kTH2F, {{100, 0, 10.}, {200, -10, 10.}});
    histos.add<TH2>("tracks/kaon/h2KaonVspNSigmaTOF", "NSigmaTOF(Ka) vs p; Momentum(p) in GeV; NSigmaTOF", HistType::kTH2F, {{100, 0, 10.}, {200, -10, 10.}});
    histos.add<TH2>("tracks/proton/h2ProtonVspNSigmaTOF", "NSigmaTOF(proton) vs p; Momentum(p) in GeV; NSigmaTOF", HistType::kTH2F, {{100, 0, 10.}, {200, -10, 10.}});
    histos.add<TH2>("tracks/deuteron/h2DeuteronVspNSigmaTOF", "NSigmaTOF(D) vs p; Momentum(p) in GeV; NSigmaTOF", HistType::kTH2F, {{100, 0, 10.}, {200, -10, 10.}});
    histos.add<TH2>("tracks/helium/h2HeliumVspNSigmaTOF", "NSigmaTOF(He) vs p; Momentum(p) in GeV; NSigmaTOF", HistType::kTH2F, {{100, 0, 10.}, {200, -10, 10.}});

    // TOF mass histograms
    histos.add<TH2>("tracks/h2TOFmassVsPt", "h2TOFmassVsPt; TOFmass; #it{p}_{T} (GeV)", HistType::kTH2F, {{600, 0., 3.}, {500, 0., 5.}});

    // TOF mass squared histograms
    histos.add<TH2>("tracks/proton/h2TOFmass2ProtonVsPt", "#Delta M^{2} (p) vs #it{p}_{T}; #Delta M^{2} (p); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{600, -3., 3.}, {800, 0., 8.}});
    histos.add<TH2>("tracks/deuteron/h2TOFmass2DeuteronVsPt", "#Delta M^{2} (d) vs #it{p}_{T}; #Delta M^{2} (d); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{1000, -5., 5.}, {800, 0., 8.}});
    histos.add<TH2>("tracks/helium/h2TOFmass2HeliumVsPt", "#Delta M^{2} (He) vs #it{p}_{T}; #Delta M^{2} (He); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{1800, -9., 9.}, {800, 0., 8.}});

    // MC histograms
    AxisSpec ptAxis = {2000, 0.f, 20.f, "#it{p}_{T} (GeV/#it{c})"};
    histos.add("spectraGen/histGenVetxZ", "PosZ generated events", HistType::kTH1F, {{2000, -20.f, 20.f, "Vertex Z (cm)"}});
    histos.add("spectraGen/histGenPtPion", "generated particles", HistType::kTH1F, {ptAxis});
    histos.add("spectraGen/histGenPtKaon", "generated particles", HistType::kTH1F, {ptAxis});
    histos.add("spectraGen/histGenPtProton", "generated particles", HistType::kTH1F, {ptAxis});
    histos.add("spectraGen/histGenPtD", "generated particles", HistType::kTH1F, {ptAxis});
    histos.add("spectraGen/histGenPtantiD", "generated particles", HistType::kTH1F, {ptAxis});
    histos.add("spectraGen/histGenPtHe", "generated particles", HistType::kTH1F, {ptAxis});
    histos.add("spectraGen/histGenPtantiHe", "generated particles", HistType::kTH1F, {ptAxis});
  }

  template <bool IsMC, typename CollisionType, typename TracksType>
  void fillHistograms(const CollisionType& event, const TracksType& tracks)
  {
    constexpr float fMassProton = 0.938272088f;
    constexpr float fMassDeuteron = 1.87561f;
    constexpr float fMassHelium = 2.80839f;
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

      // Tracks histos fill
      histos.fill(HIST("tracks/h1Eta"), track.eta());
      histos.fill(HIST("tracks/h1VarPhi"), track.phi());
      histos.fill(HIST("tracks/h2EtaVsPhi"), track.eta(), track.phi());
      histos.fill(HIST("tracks/h1pT"), track.pt());
      histos.fill(HIST("tracks/h1p"), track.p());

      //  TPC
      // if(std::abs(track.nsigTPCD()) < 5. && std::abs(track.nsigTOFD()) < 3.)
      histos.fill(HIST("tracks/h2TPCsignVsTPCmomentum"), track.tpcInnerParam(), track.tpcSignal());

      histos.fill(HIST("tracks/pion/h2PionVspNSigmaTPC"), track.p(), track.nsigTPCPi());
      histos.fill(HIST("tracks/kaon/h2KaonVspNSigmaTPC"), track.p(), track.nsigTPCKa());
      histos.fill(HIST("tracks/proton/h2ProtonVspNSigmaTPC"), track.p(), track.nsigTPCPr());
      histos.fill(HIST("tracks/deuteron/h2DeuteronVspNSigmaTPC"), track.p(), track.nsigTPCD());
      histos.fill(HIST("tracks/helium/h2HeliumVspNSigmaTPC"), track.p(), track.nsigTPC3He());

      //  TOF
      histos.fill(HIST("tracks/pion/h2PionVspNSigmaTOF"), track.p(), track.nsigTOFPi());
      histos.fill(HIST("tracks/kaon/h2KaonVspNSigmaTOF"), track.p(), track.nsigTOFKa());
      histos.fill(HIST("tracks/proton/h2ProtonVspNSigmaTOF"), track.p(), track.nsigTOFPr());
      histos.fill(HIST("tracks/deuteron/h2DeuteronVspNSigmaTOF"), track.p(), track.nsigTOFD());
      histos.fill(HIST("tracks/helium/h2HeliumVspNSigmaTOF"), track.p(), track.nsigTOF3He());

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
        track.pdgCode();
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
    constexpr int PDGPion = 211;
    constexpr int PDGKaon = 321;
    constexpr int PDGProton = 2212;
    constexpr int PDGDeuteron = 1000010020;
    constexpr int PDGHelium = 1000020030;
    nCount++;
    histos.fill(HIST("spectraGen/histGenVetxZ"), mcCollision.posZ());
    for (auto& mcParticleGen : mcParticles) {
      if (abs(mcParticleGen.y()) > std::abs(etaCut)) {
        continue;
      }

      if (std::abs(mcParticleGen.pdgCode()) == PDGPion) {
        histos.fill(HIST("spectraGen/histGenPtPion"), mcParticleGen.pt());
      }

      if (std::abs(mcParticleGen.pdgCode()) == PDGKaon) {
        histos.fill(HIST("spectraGen/histGenPtKaon"), mcParticleGen.pt());
      }

      if (std::abs(mcParticleGen.pdgCode()) == PDGProton) {
        histos.fill(HIST("spectraGen/histGenPtProton"), mcParticleGen.pt());
      }

      if (mcParticleGen.pdgCode() == PDGDeuteron) {
        histos.fill(HIST("spectraGen/histGenPtD"), mcParticleGen.pt());
      }

      if (mcParticleGen.pdgCode() == -PDGDeuteron) {
        histos.fill(HIST("spectraGen/histGenPtantiD"), mcParticleGen.pt());
      }

      if (mcParticleGen.pdgCode() == PDGHelium) {
        histos.fill(HIST("spectraGen/histGenPtHe"), mcParticleGen.pt());
      }

      if (mcParticleGen.pdgCode() == -PDGHelium) {
        histos.fill(HIST("spectraGen/histGenPtantiHe"), mcParticleGen.pt());
      }
    }
  } // Close processMCGen
  PROCESS_SWITCH(LFNucleiBATask, processMCGen, "process MC Generated", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<LFNucleiBATask>(cfgc)};
}