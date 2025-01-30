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
/// \brief this is a starting point for the third session of the tutorial
/// \author Francesca Ercolessi (francesca.ercolessi@cern.ch)
/// \since

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/v0qaanalysis.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "CommonConstants/PhysicsConstants.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using DauTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCPr, aod::pidTPCKa, aod::pidTOFPi, aod::pidTOFPr, aod::pidTOFKa>;

struct v0postprocessing {

  Configurable<float> radius{"radius", 0.5, "Radius"};
  Configurable<float> dcanegtopv{"dcanegtopv", 0.05, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", 0.05, "DCA Pos To PV"};
  Configurable<double> cospaK0s{"cospaK0s", 0.97, "K0s CosPA"};
  Configurable<double> cospaLambda{"cospaLambda", 0.995, "Lambda CosPA"};
  Configurable<float> dcav0dau{"dcav0dau", 1.0, "DCA V0 Daughters"};
  Configurable<float> rap{"rap", 0.5, "Rapidity"};
  Configurable<float> ctauK0s{"ctauK0s", 20, "C tau K0s(cm)"};
  Configurable<float> ctauLambda{"ctauLambda", 30, "C tau Lambda (cm)"};
  Configurable<float> v0rejK0s{"v0rejK0s", 0.005, "V0 rej K0s"};
  Configurable<float> v0rejLambda{"v0rejLambda", 0.01, "V0 rej K0s"};
  Configurable<float> ntpcsigma{"ntpcsigma", 5, "N sigma TPC"};
  Configurable<float> ntpcsigmaMC{"ntpcsigmaMC", 100, "N sigma TPC for MC"};
  Configurable<float> etadau{"etadau", 0.8, "Eta Daughters"};
  Configurable<float> minITShits{"minITShits", 2, "min ITS hits"};
  Configurable<bool> isMC{"isMC", 1, "isMC"};
  Configurable<bool> evSel{"evSel", 1, "evSel"};
  Configurable<bool> hasTOF2Leg{"hasTOF2Leg", 0, "hasTOF2Leg"};
  Configurable<bool> hasTOF1Leg{"hasTOF1Leg", 1, "hasTOF1Leg"};
  Configurable<float> paramArmenterosCut{"paramArmenterosCut", 0.2, "parameter Armenteros Cut"};
  Configurable<bool> doArmenterosCut{"doArmenterosCut", 1, "do Armenteros Cut"};

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {

    registry.add("hMassK0Short", ";M_{#pi^{+}#pi^{-}} [GeV/c^{2}]", {HistType::kTH1F, {{200, 0.4f, 0.6f}}});
    registry.add("hMassVsPtK0Short", ";p_{T} [GeV/c];M_{#pi^{+}#pi^{-}} [GeV/c^{2}]", {HistType::kTH2F, {{250, 0.0f, 25.0f}, {200, 0.4f, 0.6f}}});
    registry.add("hMassVsPtK0ShortVsCentFT0M", ";p_{T} [GeV/c]; CentFT0M; M_{#pi^{+}#pi^{-}} [GeV/c^{2}]", {HistType::kTH3F, {{250, 0.0f, 25.0f}, {100, 0.f, 100.f}, {200, 0.4f, 0.6f}}});
    registry.add("hMassVsPtK0ShortVsCentFV0A", ";p_{T} [GeV/c]; CentFT0M; M_{#pi^{+}#pi^{-}} [GeV/c^{2}]", {HistType::kTH3F, {{250, 0.0f, 25.0f}, {100, 0.f, 100.f}, {200, 0.4f, 0.6f}}});
    registry.add("hMassLambda", "hMassLambda", {HistType::kTH1F, {{200, 1.016f, 1.216f}}});
    registry.add("hMassVsPtLambda", "hMassVsPtLambda", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {200, 1.016f, 1.216f}}});
    registry.add("hMassVsPtLambdaVsCentFT0M", ";p_{T} [GeV/c]; CentFT0M; M_{#pi^{+}#pi^{-}} [GeV/c^{2}]", {HistType::kTH3F, {{250, 0.0f, 25.0f}, {100, 0.f, 100.f}, {200, 1.016f, 1.216f}}});
    registry.add("hMassAntiLambda", "hMassAntiLambda", {HistType::kTH1F, {{200, 1.016f, 1.216f}}});
    registry.add("hMassVsPtAntiLambda", "hMassVsPtAntiLambda", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {200, 1.016f, 1.216f}}});
    registry.add("hMassVsPtAntiLambdaVsCentFT0M", ";p_{T} [GeV/c]; CentFT0M; M_{#pi^{+}#pi^{-}} [GeV/c^{2}]", {HistType::kTH3F, {{250, 0.0f, 25.0f}, {100, 0.f, 100.f}, {200, 1.016f, 1.216f}}});

    if (isMC) {
      registry.add("hMassK0Short_MC", ";M_{#pi^{+}#pi^{-}} [GeV/c^{2}]", {HistType::kTH1F, {{200, 0.4f, 0.6f}}});
      registry.add("hMassVsPtK0Short_MC", ";p_{T} [GeV/c];M_{#pi^{+}#pi^{-}} [GeV/c^{2}]", {HistType::kTH3F, {{250, 0.0f, 25.0f}, {100, 0.f, 100.f}, {200, 0.4f, 0.6f}}});
      registry.add("hMassLambda_MC", "hMassLambda", {HistType::kTH1F, {{200, 1.016f, 1.216f}}});
      registry.add("hMassVsPtLambda_MC", ";p_{T} [GeV/c];M_{p^{+}#pi^{-}} [GeV/c^{2}]", {HistType::kTH3F, {{250, 0.0f, 25.0f}, {100, 0.f, 100.f}, {200, 1.016f, 1.216f}}});
      registry.add("hMassAntiLambda_MC", "hMassAntiLambda", {HistType::kTH1F, {{200, 1.016f, 1.216f}}});
      registry.add("hMassVsPtAntiLambda_MC", ";p_{T} [GeV/c];M_{p^{-}#pi^{+}} [GeV/c^{2}]", {HistType::kTH3F, {{250, 0.0f, 25.0f}, {100, 0.f, 100.f}, {200, 1.016f, 1.216f}}});
    }

    // QA
    registry.add("hArmenterosPodolanski", "hArmenterosPodolanski", {HistType::kTH2F, {{1000, -1.0f, 1.0f, "#alpha"}, {1000, 0.0f, 0.30f, "#it{Q}_{T}"}}});
    registry.add("hArmenterosPodolanski_Sel", "hArmenterosPodolanski_Sel", {HistType::kTH2F, {{1000, -1.0f, 1.0f, "#alpha"}, {1000, 0.0f, 0.30f, "#it{Q}_{T}"}}});
    registry.add("hK0sV0Radius", "hK0sV0Radius", {HistType::kTH1D, {{200, 0.0f, 40.0f}}});
    registry.add("hK0sCosPA", "hK0sCosPA", {HistType::kTH1F, {{100, 0.9f, 1.0f}}});
    registry.add("hK0sV0DCANegToPV", "hK0sV0DCANegToPV", {HistType::kTH1F, {{200, -1.0f, 1.0f}}});
    registry.add("hK0sV0DCAPosToPV", "hK0sV0DCAPosToPV", {HistType::kTH1F, {{200, -1.0f, 1.0f}}});
    registry.add("hK0sV0DCAV0Daughters", "hK0sV0DCAV0Daughters", {HistType::kTH1F, {{55, 0.0f, 2.20f}}});
    registry.add("hK0sCtau", "hK0sCtau", {HistType::kTH1F, {{100, 0.0f, 50.0f}}});
    registry.add("hK0sEtaDau", "hK0sEtaDau", {HistType::kTH1F, {{100, -1.0f, 1.0f}}});
    registry.add("hK0sRap", "hK0sRap", {HistType::kTH1F, {{100, -1.0f, 1.0f}}});
    registry.add("hK0sTPCNSigmaPosPi", "hK0sTPCNSigmaPosPi", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
    registry.add("hK0sTPCNSigmaNegPi", "hK0sTPCNSigmaNegPi", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});

    registry.add("hLambdaV0Radius", "hLambdaV0Radius", {HistType::kTH1D, {{200, 0.0f, 40.0f}}});
    registry.add("hLambdaCosPA", "hLambdaCosPA", {HistType::kTH1F, {{100, 0.9f, 1.0f}}});
    registry.add("hLambdaV0DCANegToPV", "hLambdaV0DCANegToPV", {HistType::kTH1F, {{200, -1.0f, 1.0f}}});
    registry.add("hLambdaV0DCAPosToPV", "hLambdaV0DCAPosToPV", {HistType::kTH1F, {{200, -1.0f, 1.0f}}});
    registry.add("hLambdaV0DCAV0Daughters", "hLambdaV0DCAV0Daughters", {HistType::kTH1F, {{55, 0.0f, 2.20f}}});
    registry.add("hLambdaCtau", "hLambdaCtau", {HistType::kTH1F, {{100, 0.0f, 50.0f}}});
    registry.add("hLambdaEtaDau", "hLambdaEtaDau", {HistType::kTH1F, {{100, -1.0f, 1.0f}}});
    registry.add("hLambdaRap", "hLambdaRap", {HistType::kTH1F, {{100, -1.0f, 1.0f}}});
    registry.add("hLambdaTPCNSigmaNegPi", "hLambdaTPCNSigmaNegPi", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
    registry.add("hLambdaTPCNSigmaPosPr", "hLambdaTPCNSigmaPosPr", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});

    registry.add("hAntiLambdaV0Radius", "hAntiLambdaV0Radius", {HistType::kTH1D, {{200, 0.0f, 40.0f}}});
    registry.add("hAntiLambdaCosPA", "hAntiLambdaCosPA", {HistType::kTH1F, {{100, 0.9f, 1.0f}}});
    registry.add("hAntiLambdaV0DCANegToPV", "hAntiLambdaV0DCANegToPV", {HistType::kTH1F, {{200, -1.0f, 1.0f}}});
    registry.add("hAntiLambdaV0DCAPosToPV", "hAntiLambdaV0DCAPosToPV", {HistType::kTH1F, {{200, -1.0f, 1.0f}}});
    registry.add("hAntiLambdaV0DCAV0Daughters", "hAntiLambdaV0DCAV0Daughters", {HistType::kTH1F, {{55, 0.0f, 2.20f}}});
    registry.add("hAntiLambdaCtau", "hAntiLambdaCtau", {HistType::kTH1F, {{100, 0.0f, 50.0f}}});
    registry.add("hAntiLambdaEtaDau", "hAntiLambdaEtaDau", {HistType::kTH1F, {{100, -1.0f, 1.0f}}});
    registry.add("hAntiLambdaRap", "hAntiLambdaRap", {HistType::kTH1F, {{100, -1.0f, 1.0f}}});
    registry.add("hAntiLambdaTPCNSigmaPosPi", "hAntiLambdaTPCNSigmaPosPi", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
    registry.add("hAntiLambdaTPCNSigmaNegPr", "hAntiLambdaTPCNSigmaNegPr", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});

    /*registry.add("TPCNSigmaPosPr", "TPCNSigmaPosPr", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
    registry.add("TPCNSigmaNegPr", "TPCNSigmaNegPr", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
    registry.add("TOFNSigmaPosPi", "TOFNSigmaPosPi", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
    registry.add("TOFNSigmaNegPi", "TOFNSigmaNegPi", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
    registry.add("TOFNSigmaPosPr", "TOFNSigmaPosPr", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
    registry.add("TOFNSigmaNegPr", "TOFNSigmaNegPr", {HistType::kTH1F, {{100, -10.0f, 10.0f}}}); */
  }

  void process(aod::MyV0Candidates const& myv0s)
  {
    for (auto& candidate : myv0s) {

      // common selections
      if (candidate.v0radius() < radius)
        continue;
      if (TMath::Abs(candidate.v0poseta()) > etadau)
        continue;
      if (TMath::Abs(candidate.v0negeta()) > etadau)
        continue;
      if (TMath::Abs(candidate.v0positshits()) < minITShits)
        continue;
      if (TMath::Abs(candidate.v0negitshits()) < minITShits)
        continue;
      if (TMath::Abs(candidate.v0dcanegtopv()) < dcanegtopv)
        continue;
      if (TMath::Abs(candidate.v0dcapostopv()) < dcapostopv)
        continue;
      if (candidate.v0dcav0daughters() > dcav0dau)
        continue;
      if (evSel && candidate.evflag() < 1)
        continue;
      if (hasTOF1Leg && !candidate.poshastof() && !candidate.neghastof())
        continue;
      if (hasTOF2Leg && (!candidate.poshastof() || !candidate.neghastof()))
        continue;

      // K0Short analysis
      if (candidate.v0cospa() > cospaK0s &&
          TMath::Abs(candidate.rapk0short()) < rap &&
          candidate.ctauk0short() < ctauK0s &&
          TMath::Abs(candidate.massk0short() - o2::constants::physics::MassK0Short) < 0.075 &&
          TMath::Abs(candidate.masslambda() - o2::constants::physics::MassLambda0) > v0rejK0s &&
          TMath::Abs(candidate.ntpcsigmanegpi()) <= ntpcsigma &&
          TMath::Abs(candidate.ntpcsigmapospi()) <= ntpcsigma) {

        registry.fill(HIST("hArmenterosPodolanski"), candidate.alpha(), candidate.qtarm());

        if (doArmenterosCut && candidate.qtarm() > (paramArmenterosCut * TMath::Abs(candidate.alpha()))) {
          registry.fill(HIST("hArmenterosPodolanski_Sel"), candidate.alpha(), candidate.qtarm());
          registry.fill(HIST("hMassK0Short"), candidate.massk0short());
          registry.fill(HIST("hMassVsPtK0Short"), candidate.v0pt(), candidate.massk0short());
          registry.fill(HIST("hMassVsPtK0ShortVsCentFT0M"), candidate.v0pt(), candidate.multft0m(), candidate.massk0short());
          registry.fill(HIST("hMassVsPtK0ShortVsCentFV0A"), candidate.v0pt(), candidate.multfv0a(), candidate.massk0short());

          // QA
          if (!isMC) {
            registry.fill(HIST("hK0sV0Radius"), candidate.v0radius());
            registry.fill(HIST("hK0sCosPA"), candidate.v0cospa());
            registry.fill(HIST("hK0sV0DCANegToPV"), candidate.v0dcanegtopv());
            registry.fill(HIST("hK0sV0DCAPosToPV"), candidate.v0dcapostopv());
            registry.fill(HIST("hK0sV0DCAV0Daughters"), candidate.v0dcav0daughters());
            registry.fill(HIST("hK0sCtau"), candidate.ctauk0short());
            registry.fill(HIST("hK0sEtaDau"), candidate.v0poseta());
            registry.fill(HIST("hK0sRap"), candidate.rapk0short());
            registry.fill(HIST("hK0sTPCNSigmaPosPi"), candidate.ntpcsigmapospi());
            registry.fill(HIST("hK0sTPCNSigmaNegPi"), candidate.ntpcsigmanegpi());
          }
        }
      }

      // Lambda analysis
      if (candidate.v0cospa() > cospaLambda &&
          TMath::Abs(candidate.raplambda()) < rap &&
          TMath::Abs(candidate.massk0short() - o2::constants::physics::MassK0Short) > v0rejLambda) {

        // Lambda
        if (TMath::Abs(candidate.ntpcsigmanegpi()) <= ntpcsigma && TMath::Abs(candidate.ntpcsigmapospr()) <= ntpcsigma &&
            candidate.ctaulambda() < ctauLambda &&
            TMath::Abs(candidate.masslambda() - o2::constants::physics::MassLambda0) < 0.075) {

          registry.fill(HIST("hMassLambda"), candidate.masslambda());
          registry.fill(HIST("hMassVsPtLambda"), candidate.v0pt(), candidate.masslambda());
          registry.fill(HIST("hMassVsPtLambdaVsCentFT0M"), candidate.v0pt(), candidate.multft0m(), candidate.masslambda());

          // QA
          if (!isMC) {
            registry.fill(HIST("hLambdaV0Radius"), candidate.v0radius());
            registry.fill(HIST("hLambdaCosPA"), candidate.v0cospa());
            registry.fill(HIST("hLambdaV0DCANegToPV"), candidate.v0dcanegtopv());
            registry.fill(HIST("hLambdaV0DCAPosToPV"), candidate.v0dcapostopv());
            registry.fill(HIST("hLambdaV0DCAV0Daughters"), candidate.v0dcav0daughters());
            registry.fill(HIST("hLambdaCtau"), candidate.ctaulambda());
            registry.fill(HIST("hLambdaEtaDau"), candidate.v0poseta());
            registry.fill(HIST("hLambdaRap"), candidate.raplambda());
            registry.fill(HIST("hLambdaTPCNSigmaPosPr"), candidate.ntpcsigmapospr());
            registry.fill(HIST("hLambdaTPCNSigmaNegPi"), candidate.ntpcsigmanegpi());
          }
        }
        // AntiLambda
        if (TMath::Abs(candidate.ntpcsigmanegpr()) <= ntpcsigma && TMath::Abs(candidate.ntpcsigmapospi()) <= ntpcsigma &&
            candidate.ctauantilambda() < ctauLambda &&
            TMath::Abs(candidate.massantilambda() - o2::constants::physics::MassLambda0) < 0.075) {

          registry.fill(HIST("hMassAntiLambda"), candidate.massantilambda());
          registry.fill(HIST("hMassVsPtAntiLambda"), candidate.v0pt(), candidate.massantilambda());
          registry.fill(HIST("hMassVsPtAntiLambdaVsCentFT0M"), candidate.v0pt(), candidate.multft0m(), candidate.massantilambda());

          // QA
          if (!isMC) {
            registry.fill(HIST("hAntiLambdaV0Radius"), candidate.v0radius());
            registry.fill(HIST("hAntiLambdaCosPA"), candidate.v0cospa());
            registry.fill(HIST("hAntiLambdaV0DCANegToPV"), candidate.v0dcanegtopv());
            registry.fill(HIST("hAntiLambdaV0DCAPosToPV"), candidate.v0dcapostopv());
            registry.fill(HIST("hAntiLambdaV0DCAV0Daughters"), candidate.v0dcav0daughters());
            registry.fill(HIST("hAntiLambdaCtau"), candidate.ctauantilambda());
            registry.fill(HIST("hAntiLambdaEtaDau"), candidate.v0poseta());
            registry.fill(HIST("hAntiLambdaRap"), candidate.raplambda());
            registry.fill(HIST("hAntiLambdaTPCNSigmaNegPr"), candidate.ntpcsigmanegpr());
            registry.fill(HIST("hAntiLambdaTPCNSigmaPosPi"), candidate.ntpcsigmapospi());
          }
        }
      }

      if (isMC) {

        if (candidate.isphysprimary() == 0)
          continue;

        // K0Short analysis
        if (candidate.v0cospa() > cospaK0s &&
            TMath::Abs(candidate.rapk0short()) < rap &&
            candidate.ctauk0short() < ctauK0s &&
            TMath::Abs(candidate.massk0short() - o2::constants::physics::MassK0Short) < 0.075 &&
            TMath::Abs(candidate.masslambda() - o2::constants::physics::MassLambda0) > v0rejK0s &&
            TMath::Abs(candidate.ntpcsigmanegpi()) <= ntpcsigmaMC &&
            TMath::Abs(candidate.ntpcsigmapospi()) <= ntpcsigmaMC &&
            (candidate.pdgcode() == 310) && candidate.isdauk0short()) {

          registry.fill(HIST("hArmenterosPodolanski"), candidate.alpha(), candidate.qtarm());

          if (doArmenterosCut && candidate.qtarm() > (paramArmenterosCut * TMath::Abs(candidate.alpha()))) {
            registry.fill(HIST("hArmenterosPodolanski_Sel"), candidate.alpha(), candidate.qtarm());
            registry.fill(HIST("hMassK0Short_MC"), candidate.massk0short());
            registry.fill(HIST("hMassVsPtK0Short_MC"), candidate.v0pt(), candidate.multft0m(), candidate.massk0short());

            registry.fill(HIST("hK0sV0Radius"), candidate.v0radius());
            registry.fill(HIST("hK0sCosPA"), candidate.v0cospa());
            registry.fill(HIST("hK0sV0DCANegToPV"), candidate.v0dcanegtopv());
            registry.fill(HIST("hK0sV0DCAPosToPV"), candidate.v0dcapostopv());
            registry.fill(HIST("hK0sV0DCAV0Daughters"), candidate.v0dcav0daughters());
            registry.fill(HIST("hK0sCtau"), candidate.ctauk0short());
            registry.fill(HIST("hK0sEtaDau"), candidate.v0poseta());
            registry.fill(HIST("hK0sRap"), candidate.rapk0short());
            registry.fill(HIST("hK0sTPCNSigmaPosPi"), candidate.ntpcsigmapospi());
            registry.fill(HIST("hK0sTPCNSigmaNegPi"), candidate.ntpcsigmanegpi());
          }
        } // k0

        // Lambda analysis
        if (candidate.v0cospa() > cospaLambda &&
            TMath::Abs(candidate.raplambda()) < rap &&
            TMath::Abs(candidate.massk0short() - o2::constants::physics::MassK0Short) > v0rejLambda) {

          // Lambda
          if (TMath::Abs(candidate.ntpcsigmanegpi()) <= ntpcsigmaMC && TMath::Abs(candidate.ntpcsigmapospr()) <= ntpcsigmaMC &&
              candidate.ctaulambda() < ctauLambda &&
              TMath::Abs(candidate.masslambda() - o2::constants::physics::MassLambda0) < 0.075 &&
              candidate.pdgcode() == 3122 && candidate.isdaulambda()) {

            registry.fill(HIST("hMassLambda_MC"), candidate.masslambda());
            registry.fill(HIST("hMassVsPtLambda_MC"), candidate.v0pt(), candidate.multft0m(), candidate.masslambda());

            registry.fill(HIST("hLambdaV0Radius"), candidate.v0radius());
            registry.fill(HIST("hLambdaCosPA"), candidate.v0cospa());
            registry.fill(HIST("hLambdaV0DCANegToPV"), candidate.v0dcanegtopv());
            registry.fill(HIST("hLambdaV0DCAPosToPV"), candidate.v0dcapostopv());
            registry.fill(HIST("hLambdaV0DCAV0Daughters"), candidate.v0dcav0daughters());
            registry.fill(HIST("hLambdaCtau"), candidate.ctaulambda());
            registry.fill(HIST("hLambdaEtaDau"), candidate.v0poseta());
            registry.fill(HIST("hLambdaRap"), candidate.raplambda());
            registry.fill(HIST("hLambdaTPCNSigmaPosPr"), candidate.ntpcsigmapospr());
            registry.fill(HIST("hLambdaTPCNSigmaNegPi"), candidate.ntpcsigmanegpi());
          }
          // AntiLambda
          if (TMath::Abs(candidate.ntpcsigmanegpr()) <= ntpcsigmaMC && TMath::Abs(candidate.ntpcsigmapospi()) <= ntpcsigmaMC &&
              candidate.ctauantilambda() < ctauLambda &&
              TMath::Abs(candidate.massantilambda() - o2::constants::physics::MassLambda0) < 0.075 &&
              candidate.pdgcode() == -3122 && candidate.isdauantilambda()) {

            registry.fill(HIST("hMassAntiLambda_MC"), candidate.massantilambda());
            registry.fill(HIST("hMassVsPtAntiLambda_MC"), candidate.v0pt(), candidate.multft0m(), candidate.massantilambda());

            registry.fill(HIST("hAntiLambdaV0Radius"), candidate.v0radius());
            registry.fill(HIST("hAntiLambdaCosPA"), candidate.v0cospa());
            registry.fill(HIST("hAntiLambdaV0DCANegToPV"), candidate.v0dcanegtopv());
            registry.fill(HIST("hAntiLambdaV0DCAPosToPV"), candidate.v0dcapostopv());
            registry.fill(HIST("hAntiLambdaV0DCAV0Daughters"), candidate.v0dcav0daughters());
            registry.fill(HIST("hAntiLambdaCtau"), candidate.ctauantilambda());
            registry.fill(HIST("hAntiLambdaEtaDau"), candidate.v0poseta());
            registry.fill(HIST("hAntiLambdaRap"), candidate.raplambda());
            registry.fill(HIST("hAntiLambdaTPCNSigmaPosPi"), candidate.ntpcsigmapospi());
            registry.fill(HIST("hAntiLambdaTPCNSigmaNegPr"), candidate.ntpcsigmanegpr());
          }
        } // lambda
      } // is MC
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<v0postprocessing>(cfgc, TaskName{"lf-v0postprocessing"})};
}
