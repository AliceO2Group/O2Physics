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
//
// Analysis task to produce smeared pt,eta,phi for electrons/muons in dilepton analysis
//    Please write to: daiki.sekihata@cern.ch

#include <iostream>
#include <vector>
#include <TMath.h>
#include <TH1D.h>
#include <TString.h>
#include <TGrid.h>
#include <TObjArray.h>
#include <TFile.h>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/ASoA.h"
#include "Framework/DataTypes.h"
#include "Framework/HistogramRegistry.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"

using std::cout;
using std::endl;
using std::string;

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

struct ApplySmearing {
  Produces<aod::SmearedTracks> smearedtrack;

  Configurable<std::string> fConfigResFileName{"cfgResFileName", "", "name of resolution file"};
  Configurable<std::string> fConfigResPtHistName{"cfgResPtHistName", "RelPtResArrCocktail", "histogram name for pt in resolution file"};
  Configurable<std::string> fConfigResEtaHistName{"cfgResEtaHistName", "EtaResArr", "histogram name for eta in resolution file"};
  Configurable<std::string> fConfigResPhiPosHistName{"cfgResPhiPosHistName", "PhiPosResArr", "histogram name for phi pos in resolution file"};
  Configurable<std::string> fConfigResPhiNegHistName{"cfgResPhiNegHistName", "PhiEleResArr", "hisogram for phi neg in resolution file"};

  TObjArray* fArrResoPt;
  TObjArray* fArrResoEta;
  TObjArray* fArrResoPhi_Pos;
  TObjArray* fArrResoPhi_Neg;

  void init(InitContext& context)
  {
    if (TString(fConfigResFileName).BeginsWith("alien://")) {
      TGrid::Connect("alien://");
    }

    // get resolutoin histo
    LOGP(info, "Set Resolution histo");
    // Get Resolution map
    TFile* fFile = TFile::Open(TString(fConfigResFileName));
    if (!fFile) {
      LOGP(error, "Could not open Resolution file {}", TString(fConfigResFileName));
      return;
    }
    TObjArray* ArrResoPt = nullptr;
    if (fFile->GetListOfKeys()->Contains(TString(fConfigResPtHistName))) {
      ArrResoPt = reinterpret_cast<TObjArray*>(fFile->Get(TString(fConfigResPtHistName)));
    } else {
      LOGP(error, "Could not open {} from file {}", TString(fConfigResPtHistName), TString(fConfigResFileName));
    }

    TObjArray* ArrResoEta = nullptr;
    if (fFile->GetListOfKeys()->Contains(TString(fConfigResEtaHistName))) {
      ArrResoEta = reinterpret_cast<TObjArray*>(fFile->Get(TString(fConfigResEtaHistName)));
    } else {
      LOGP(error, "Could not open {} from file {}", TString(fConfigResEtaHistName), TString(fConfigResFileName));
    }

    TObjArray* ArrResoPhi_Pos = nullptr;
    if (fFile->GetListOfKeys()->Contains(TString(fConfigResPhiPosHistName))) {
      ArrResoPhi_Pos = reinterpret_cast<TObjArray*>(fFile->Get(TString(fConfigResPhiPosHistName)));
    } else {
      LOGP(error, "Could not open {} from file {}", TString(fConfigResPhiPosHistName), TString(fConfigResFileName));
    }

    TObjArray* ArrResoPhi_Neg = nullptr;
    if (fFile->GetListOfKeys()->Contains(TString(fConfigResPhiNegHistName))) {
      ArrResoPhi_Neg = reinterpret_cast<TObjArray*>(fFile->Get(TString(fConfigResPhiNegHistName)));
    } else {
      LOGP(error, "Could not open {} from file {}", TString(fConfigResPhiNegHistName), TString(fConfigResFileName));
    }

    fArrResoPt = ArrResoPt;
    fArrResoEta = ArrResoEta;
    fArrResoPhi_Pos = ArrResoPhi_Pos;
    fArrResoPhi_Neg = ArrResoPhi_Neg;
    fFile->Close();
  }

  template <typename TTracksMC>
  void applySmearing(TTracksMC const& tracksMC)
  {
    for (auto& mctrack : tracksMC) {
      float ptgen = mctrack.pt();
      float etagen = mctrack.eta();
      float phigen = mctrack.phi();

      if (abs(mctrack.pdgCode()) == 11 || abs(mctrack.pdgCode()) == 13) {
        // apply smearing for electrons or muons.

        // smear pt
        Int_t ptbin = reinterpret_cast<TH2D*>(fArrResoPt->At(0))->GetXaxis()->FindBin(ptgen);
        if (ptbin < 1) {
          ptbin = 1;
        }
        if (ptbin > fArrResoPt->GetLast()) {
          ptbin = fArrResoPt->GetLast();
        }
        float smearing = 0.;
        TH1D* thisHist = reinterpret_cast<TH1D*>(fArrResoPt->At(ptbin));
        if (thisHist->GetEntries() > 0) {
          smearing = thisHist->GetRandom() * ptgen;
        }
        float ptsmeared = ptgen - smearing;

        // smear eta
        ptbin = reinterpret_cast<TH2D*>(fArrResoEta->At(0))->GetXaxis()->FindBin(ptgen);
        if (ptbin < 1) {
          ptbin = 1;
        }
        if (ptbin > fArrResoEta->GetLast()) {
          ptbin = fArrResoEta->GetLast();
        }
        smearing = 0.;
        thisHist = reinterpret_cast<TH1D*>(fArrResoEta->At(ptbin));
        if (thisHist->GetEntries() > 0) {
          smearing = thisHist->GetRandom();
        }
        float etasmeared = etagen - smearing;

        // smear phi
        ptbin = reinterpret_cast<TH2D*>(fArrResoPhi_Pos->At(0))->GetXaxis()->FindBin(ptgen);
        if (ptbin < 1) {
          ptbin = 1;
        }
        if (ptbin > fArrResoPhi_Pos->GetLast()) {
          ptbin = fArrResoPhi_Pos->GetLast();
        }
        smearing = 0.;
        if (mctrack.pdgCode() < 0) { // positron: -11
          thisHist = reinterpret_cast<TH1D*>(fArrResoPhi_Pos->At(ptbin));
        } else if (mctrack.pdgCode() > 0) { // electron: 11
          thisHist = reinterpret_cast<TH1D*>(fArrResoPhi_Neg->At(ptbin));
        }
        if (thisHist->GetEntries() > 0) {
          smearing = thisHist->GetRandom();
        }
        float phismeared = phigen - smearing;

        smearedtrack(ptsmeared, etasmeared, phismeared);
      } else {
        // don't apply smearing
        smearedtrack(ptgen, etagen, phigen);
      }
    }
  }

  void processMCanalysis(ReducedMCTracks const& tracksMC)
  {
    applySmearing(tracksMC);
  }

  void processCocktail(aod::McParticles_001 const& tracksMC)
  {
    applySmearing(tracksMC);
  }

  void processDummy(aod::McParticles_001 const& tracksMC) {}

  PROCESS_SWITCH(ApplySmearing, processMCanalysis, "Run for MC analysis", false);
  PROCESS_SWITCH(ApplySmearing, processCocktail, "Run for cocktail analysis", false);
  PROCESS_SWITCH(ApplySmearing, processDummy, "Dummy process function", true);
};

struct CheckSmearing {
  using MyTracks = soa::Join<ReducedMCTracks, SmearedTracks>;
  HistogramRegistry registry{
    "registry",
    {
      {"hCorrelation_Pt", "pT correlation", {HistType::kTH2F, {{1000, 0.0f, 10.0f}, {1000, 0.0f, 10.0f}}}},
      {"hCorrelation_Eta", "eta correlation", {HistType::kTH2F, {{200, -1.0f, +1.0f}, {200, -1.0f, +1.0f}}}},
      {"hCorrelation_Phi", "phi correlation", {HistType::kTH2F, {{100, 0.0f, TMath::TwoPi()}, {100, 0.0f, TMath::TwoPi()}}}},
    },
  };

  void processCheck(MyTracks const& tracksMC)
  {
    for (auto& mctrack : tracksMC) {
      if (abs(mctrack.pdgCode()) != 11 && abs(mctrack.pdgCode()) != 13) {
        continue;
      }
      registry.fill(HIST("hCorrelation_Pt"), mctrack.pt(), mctrack.ptSmeared());
      registry.fill(HIST("hCorrelation_Eta"), mctrack.eta(), mctrack.etaSmeared());
      registry.fill(HIST("hCorrelation_Phi"), mctrack.phi(), mctrack.phiSmeared());
    } // end of mctrack loop
  }

  void processDummy(MyTracks const& tracksMC) {}

  PROCESS_SWITCH(CheckSmearing, processCheck, "Run for MC analysis", false);
  PROCESS_SWITCH(CheckSmearing, processDummy, "Dummy process function", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ApplySmearing>(cfgc, TaskName{"apply-smearing"}),
    adaptAnalysisTask<CheckSmearing>(cfgc, TaskName{"check-smearing"})};
}
