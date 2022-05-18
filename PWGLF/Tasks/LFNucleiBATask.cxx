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

struct LFNucleiDeuteronTask {
  OutputObj<TH1F> h1VtxZ{TH1F("h1VtxZ", "V_{z};V_{z} (in cm); counts", 3000, -15, 15)};
  OutputObj<TH1F> h1CentV0M{TH1F("h1CentV0M", "V0M; Multiplicity; counts", 27000, 0, 27000)};
  OutputObj<TH2F> h2TPCsignVsTPCmomentum{TH2F("h2TPCsignVsTPCmomentum", "-dE/dX vs p; p (GeV/c); -dE/dx (a.u.)", 500, 0.0, 5.0, 1000, 0.0, 1E3)};
  OutputObj<TH2F> h2TOFbetaVsP{TH2F("h2TOFbetaVsP", "#beta (TOF) vs p; p (GeV/c); #beta", 500, 0.0, 5, 1200, 0.0, 1.2)};
  OutputObj<TH2F> h2PionVspNSigmaTPC{TH2F("h2PionVspNSigmaTPC", "NSigmaTPC(pi) vs p; Momentum(p) in GeV; NSigmaTPC", 100, 0, 10., 200, -10, 10.)};
  OutputObj<TH2F> h2KaonVspNSigmaTPC{TH2F("h2KaonVspNSigmaTPC", "NSigmaTPC(Ka) vs p; Momentum(p) in GeV; NSigmaTPC", 100, 0, 10., 200, -10, 10.)};
  OutputObj<TH2F> h2ProtonVspNSigmaTPC{TH2F("h2ProtonVspNSigmaTPC", "NSigmaTPC(proton) vs p; Momentum(p) in GeV; NSigmaTPC", 100, 0, 10., 200, -10, 10.)};
  OutputObj<TH2F> h2DeuteronVspNSigmaTPC{TH2F("h2DeuteronVspNSigmaTPC", "NSigmaTPC(D) vs p; Momentum(p) in GeV; NSigmaTPC", 100, 0, 10., 200, -10, 10.)};
  OutputObj<TH2F> h2HeliumVspNSigmaTPC{TH2F("h2HeliumVspNSigmaTPC", "NSigmaTPC(He) vs p; Momentum(p) in GeV; NSigmaTPC", 100, 0, 10., 200, -10, 10.)};
  OutputObj<TH2F> h2PionVspNSigmaTOF{TH2F("h2PionVspNSigmaTOF", "NSigmaTOF(pi) vs p; Momentum(p) in GeV; NSigmaTOF", 100, 0, 10., 200, -10, 10.)};
  OutputObj<TH2F> h2KaonVspNSigmaTOF{TH2F("h2KaonVspNSigmaTOF", "NSigmaTOF(Ka) vs p; Momentum(p) in GeV; NSigmaTOF", 100, 0, 10., 200, -10, 10.)};
  OutputObj<TH2F> h2ProtonVspNSigmaTOF{TH2F("h2ProtonVspNSigmaTOF", "NSigmaTOF(proton) vs p; Momentum(p) in GeV; NSigmaTOF", 100, 0, 10., 200, -10, 10.)};
  OutputObj<TH2F> h2DeuteronVspNSigmaTOF{TH2F("h2DeuteronVspNSigmaTOF", "NSigmaTOF(D) vs p; Momentum(p) in GeV; NSigmaTOF", 100, 0, 10., 200, -10, 10.)};
  OutputObj<TH2F> h2HeliumVspNSigmaTOF{TH2F("h2HeliumVspNSigmaTOF", "NSigmaTOF(He) vs p; Momentum(p) in GeV; NSigmaTOF", 100, 0, 10., 200, -10, 10.)};
  OutputObj<TH1F> h1ProtonSpectra{TH1F("h1ProtonSpectra", "#it{p}_{T} (p); #it{p}_{T} (GeV/#it{c}); counts", 800, 0, 8)};
  OutputObj<TH1F> h1DeuteronSpectra{TH1F("h1DeuteronSpectra", "#it{p}_{T} (d); #it{p}_{T} (GeV/#it{c}); counts", 800, 0, 8)};
  OutputObj<TH1F> h1HeliumSpectra{TH1F("h1HeliumSpectra", "#it{p}_{T} (He); #it{p}_{T} (GeV/#it{c}); counts", 800, 0, 8)};

  // trackQA
  // OutputObj<TH1F> h1Y{TH1F("h1Y", "rapidity; Y; counts", 200, -10.0, 10)};
  OutputObj<TH1F> h1Eta{TH1F("h1Eta", "pseudoRapidity; #eta; counts", 200, -1.0, 1.0)};
  OutputObj<TH1F> h1VarPhi{TH1F("h1VarPhi", "#phi; #phi; counts", 80, -0.5, 7.5)};
  OutputObj<TH2F> h2EtaVsPhi{TH2F("h2EtaVsPhi", "#eta (TOF) vs #phi; #eta; #phi", 200, -1.0, 1.0, 80, -0.5, 7.5)};
  OutputObj<TH1F> h1pT{TH1F("h1pT", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", 1000, 0., 10)};
  OutputObj<TH1F> h1p{TH1F("h1p", "Track momentum; p (GeV/#it{c}); counts", 1000, 0., 10.)};

  // TPCQA
  OutputObj<TH1F> h1TPCncr{TH1F("h1TPCncr", "number of crossed rows in TPC; TPCncr; counts", 150, 60, 170)};
  OutputObj<TH1F> h1rTPC{TH1F("h1rTPC", "ratio of ncr over findable in TPC; rTPC; counts", 200, 0.9, 1.8)};
  OutputObj<TH1F> h1chi2ITS{TH1F("h1chi2ITS", "#chi^{2}_{ITS}/n_{ITS}; #chi^{2}_{ITS}/n_{ITS};counts", 51, -0.5, 50.5)};
  OutputObj<TH1F> h1chi2TPC{TH1F("h1chi2TPC", "#chi^{2}_{ITS}/n_{TPC}; #chi^{2}_{ITS}/n_{TPC}; counts", 11, -0.5, 10.5)};

  // TOF mass histograms
  OutputObj<TH2F> h2TOFmassVsPt{TH2F("h2TOFmassVsPt", "h2TOFmassVsPt; TOFmass; #it{p}_{T} (GeV)", 600, 0., 3., 500, 0., 5.)};

  // TOF mass squared histograms
  OutputObj<TH2F> h2TOFmass2ProtonVsPt{TH2F("h2TOFmass2ProtonVsPt", "#Delta M^{2} (p) vs #it{p}_{T}; #Delta M^{2} (p); #it{p}_{T} (GeV/#it{c})", 600, -3., 3., 800, 0., 8.)};
  OutputObj<TH2F> h2TOFmass2DeuteronVsPt{TH2F("h2TOFmass2DeuteronVsPt", "#Delta M^{2} (d) vs #it{p}_{T}; #Delta M^{2} (d); #it{p}_{T} (GeV/#it{c})", 1000, -5., 5., 800, 0., 8.)};
  OutputObj<TH2F> h2TOFmass2HeliumVsPt{TH2F("h2TOFmass2HeliumVsPt", "#Delta M^{2} (He) vs #it{p}_{T}; #Delta M^{2} (He); #it{p}_{T} (GeV/#it{c})", 1800, -9., 9., 800, 0., 8.)};
  Configurable<float> nsigmaTPCcut{"nsigmaTPCcut", 5.f, "Value of the Nsigma TPC cut"};
  Configurable<float> nsigmaTOFcut{"nsigmaTOFcut", 5.f, "Value of the Nsigma TOF cut"};

  void init(o2::framework::InitContext&)
  {
  }

  void process(o2::aod::LfCandNucleusFullEvents::iterator const& event, o2::aod::LfCandNucleusFull const& tracks)
  {
    constexpr float fMassProton = 0.938272088f;
    constexpr float fMassDeuteron = 1.87561f;
    constexpr float fMassHelium = 2.80839f;
    float gamma = 0., massTOF = 0.;

    // for(auto& coll:events)
    h1VtxZ->Fill(event.posZ());
    h1CentV0M->Fill(event.v0m());
    for (auto& track : tracks) {
      // LOG(info)<<"\n collisionId ============>"<<track.collisionId();
      h1TPCncr->Fill(track.ncrTPC());
      h1rTPC->Fill(track.rTPC());
      h1chi2TPC->Fill(track.chi2TPC());
      h1chi2ITS->Fill(track.chi2ITS());

      h1Eta->Fill(track.eta());
      h1VarPhi->Fill(track.phi());
      h2EtaVsPhi->Fill(track.eta(), track.phi());
      h1pT->Fill(track.pt());
      h1p->Fill(track.p());
      h2PionVspNSigmaTPC->Fill(track.p(), track.nsigTPCPi());
      h2KaonVspNSigmaTPC->Fill(track.p(), track.nsigTPCKa());
      h2ProtonVspNSigmaTPC->Fill(track.p(), track.nsigTPCPr());
      h2DeuteronVspNSigmaTPC->Fill(track.p(), track.nsigTPCD());
      h2HeliumVspNSigmaTPC->Fill(track.p(), track.nsigTPC3He());
      h2PionVspNSigmaTOF->Fill(track.p(), track.nsigTOFPi());
      h2KaonVspNSigmaTOF->Fill(track.p(), track.nsigTOFKa());
      h2ProtonVspNSigmaTOF->Fill(track.p(), track.nsigTOFPr());
      h2DeuteronVspNSigmaTOF->Fill(track.p(), track.nsigTOFD());
      h2HeliumVspNSigmaTOF->Fill(track.p(), track.nsigTOF3He());
      // if(std::abs(track.nsigTPCD()) < 5. && std::abs(track.nsigTOFD()) < 3.)
      h2TPCsignVsTPCmomentum->Fill(track.tpcInnerParam(), track.tpcSignal());

      if (track.hasTOF()) {
        h2TOFbetaVsP->Fill(track.p(), track.beta());
        if ((track.beta() * track.beta()) < 1.) {
          gamma = 1.f / TMath::Sqrt(1.f - (track.beta() * track.beta()));
          massTOF = track.p() / TMath::Sqrt(gamma * gamma - 1.f);
        } else {
          massTOF = -99.f;
        }
        h2TOFmassVsPt->Fill(massTOF, track.pt());
        h2TOFmass2ProtonVsPt->Fill(massTOF * massTOF - fMassProton * fMassProton, track.pt());
        h2TOFmass2DeuteronVsPt->Fill(massTOF * massTOF - fMassDeuteron * fMassDeuteron, track.pt());
        h2TOFmass2HeliumVsPt->Fill(massTOF * massTOF - fMassHelium * fMassHelium, track.pt());
      }

      if (std::abs(track.nsigTPCPr()) < nsigmaTPCcut)
        h1ProtonSpectra->Fill(track.pt());

      if (std::abs(track.nsigTPCD()) < nsigmaTPCcut)
        h1DeuteronSpectra->Fill(track.pt());

      if (std::abs(track.nsigTPC3He()) < nsigmaTPCcut)
        h1HeliumSpectra->Fill(track.pt());
    }
    // LOG(info)<<"Vertex Z ==="<<coll.posZ();
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<LFNucleiDeuteronTask>(cfgc)};
}
