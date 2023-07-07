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
/// \brief QA task for V0 analysis using derived data
///
/// \author Francesca Ercolessi (francesca.ercolessi@cern.ch)

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/v0qaanalysis.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using DauTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFPr>;
using DauTracksMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::McTrackLabels, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFPr>;

struct v0qaanalysis {

  // Produces
  Produces<aod::MyV0Candidates> myv0s;

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    registry.add("hNEvents", "hNEvents", {HistType::kTH1I, {{1, 0.f, 1.f}}});
    registry.add("hCentFT0M", "hCentFT0M", {HistType::kTH1F, {{1000, 0.f, 100.f}}});
    registry.add("hCentFV0A", "hCentFV0A", {HistType::kTH1F, {{1000, 0.f, 100.f}}});
    if (isMC) {
      registry.add("hNEventsMC_AllColl", "hNEventsMC_AllColl", {HistType::kTH1I, {{1, 0.f, 1.f}}});
      registry.add("hNEventsMC_RecoColl", "hNEventsMC_RecoColl", {HistType::kTH1I, {{1, 0.f, 1.f}}});
      registry.add("Reconstructed_MCRecoColl_K0Short", "Reconstructed_MCRecoColl_K0Short", {HistType::kTH1F, {{250, 0.f, 25.f}}});
      registry.add("Reconstructed_MCRecoColl_Lambda", "Reconstructed_MCRecoColl_Lambda", {HistType::kTH1F, {{250, 0.f, 25.f}}});
      registry.add("Reconstructed_MCRecoColl_AntiLambda", "Reconstructed_MCRecoColl_AntiLambda", {HistType::kTH1F, {{250, 0.f, 25.f}}});
      registry.add("Generated_MCRecoColl_K0Short", "Generated_MCRecoColl_K0Short", {HistType::kTH1F, {{250, 0.f, 25.f}}});
      registry.add("Generated_MCRecoColl_Lambda", "Generated_MCRecoColl_Lambda", {HistType::kTH1F, {{250, 0.f, 25.f}}});
      registry.add("Generated_MCRecoColl_AntiLambda", "Generated_MCRecoColl_AntiLambda", {HistType::kTH1F, {{250, 0.f, 25.f}}});
      registry.add("Generated_MCAllColl_K0Short", "Generated_MCAllColl_K0Short", {HistType::kTH1F, {{250, 0.f, 25.f}}});
      registry.add("Generated_MCAllColl_Lambda", "Generated_MCAllColl_Lambda", {HistType::kTH1F, {{250, 0.f, 25.f}}});
      registry.add("Generated_MCAllColl_AntiLambda", "Generated_MCAllColl_AntiLambda", {HistType::kTH1F, {{250, 0.f, 25.f}}});
    }
  }

  // Event selection criteria
  Configurable<float> cutzvertex{"cutzvertex", 15.0f, "Accepted z-vertex range (cm)"};
  Configurable<bool> sel8{"sel8", 1, "Apply sel8 event selection"};

  // V0 selection criteria
  Configurable<double> v0cospa{"v0cospa", 0.97, "V0 CosPA"};
  Configurable<float> dcav0dau{"dcav0dau", 1.5, "DCA V0 Daughters"};
  Configurable<float> dcanegtopv{"dcanegtopv", 0.06, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", 0.06, "DCA Pos To PV"};
  Configurable<float> v0radius{"v0radius", 0.5, "Radius"};
  Configurable<float> rapidity{"rapidity", 0.5, "Rapidity"};
  Configurable<float> etadau{"etadau", 0.8, "Eta Daughters"};
  Configurable<bool> isMC{"isMC", 0, "Is MC"};

  Filter preFilterV0 = nabs(aod::v0data::dcapostopv) > dcapostopv&&
                                                         nabs(aod::v0data::dcanegtopv) > dcanegtopv&& aod::v0data::dcaV0daughters < dcav0dau;

  void processData(soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms, aod::CentFV0As>::iterator const& collision,
                   soa::Filtered<aod::V0Datas> const& V0s, DauTracks const& tracks)
  {
    // Event selection
    if (sel8 && !collision.sel8()) {
      return;
    }
    if (TMath::Abs(collision.posZ()) > cutzvertex) {
      return;
    }

    registry.fill(HIST("hNEvents"), 0.5);
    registry.fill(HIST("hCentFT0M"), collision.centFT0M());
    registry.fill(HIST("hCentFV0A"), collision.centFV0A());

    for (auto& v0 : V0s) { // loop over V0s

      float ctauLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * RecoDecay::getMassPDG(3122);
      float ctauAntiLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * RecoDecay::getMassPDG(-3122);
      float ctauK0s = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * RecoDecay::getMassPDG(310);

      int posITSNhits = 0, negITSNhits = 0;
      for (unsigned int i = 0; i < 7; i++) {
        if (v0.posTrack_as<DauTracks>().itsClusterMap() & (1 << i)) {
          posITSNhits++;
        }
        if (v0.negTrack_as<DauTracks>().itsClusterMap() & (1 << i)) {
          negITSNhits++;
        }
      }

      int lPDG = 0;
      bool isPhysicalPrimary = isMC;

      if (v0.v0radius() > v0radius &&
          v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > v0cospa &&
          TMath::Abs(v0.posTrack_as<DauTracks>().eta()) < etadau &&
          TMath::Abs(v0.negTrack_as<DauTracks>().eta()) < etadau) {

        // Fill table
        myv0s(v0.globalIndex(), v0.pt(), v0.yLambda(), v0.yK0Short(),
              v0.mLambda(), v0.mAntiLambda(), v0.mK0Short(),
              v0.v0radius(), v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()),
              v0.dcapostopv(), v0.dcanegtopv(), v0.dcaV0daughters(),
              v0.posTrack_as<DauTracks>().eta(), v0.negTrack_as<DauTracks>().eta(),
              v0.posTrack_as<DauTracks>().phi(), v0.negTrack_as<DauTracks>().phi(),
              posITSNhits, negITSNhits, ctauLambda, ctauAntiLambda, ctauK0s,
              v0.negTrack_as<DauTracks>().tpcNSigmaPr(), v0.posTrack_as<DauTracks>().tpcNSigmaPr(),
              v0.negTrack_as<DauTracks>().tpcNSigmaPi(), v0.posTrack_as<DauTracks>().tpcNSigmaPi(),
              v0.negTrack_as<DauTracks>().tofNSigmaPr(), v0.posTrack_as<DauTracks>().tofNSigmaPr(),
              v0.negTrack_as<DauTracks>().tofNSigmaPi(), v0.posTrack_as<DauTracks>().tofNSigmaPi(),
              v0.posTrack_as<DauTracks>().hasTOF(), v0.negTrack_as<DauTracks>().hasTOF(), lPDG, isPhysicalPrimary,
              collision.centFT0M(), collision.centFV0A());
      }
    }
  }
  PROCESS_SWITCH(v0qaanalysis, processData, "Process data", true);

  Preslice<soa::Join<aod::V0Datas, aod::McV0Labels>> perCol = aod::track::collisionId;
  Preslice<aod::McParticles> perMCCol = aod::mcparticle::mcCollisionId;

  SliceCache cache1;
  SliceCache cache2;

  void processMC(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels> const& collisions,
                 aod::McCollisions const& mcCollisions,
                 soa::Join<aod::V0Datas, aod::McV0Labels> const& V0s,
                 aod::McParticles const& mcParticles, DauTracksMC const& tracks)
  {
    for (const auto& collision : collisions) {
      // Event selection
      if (sel8 && !collision.sel8()) {
        continue;
      }
      if (TMath::Abs(collision.posZ()) > cutzvertex) {
        continue;
      }
      if (!collision.has_mcCollision()) {
        continue;
      }

      registry.fill(HIST("hNEventsMC_RecoColl"), 0.5);

      auto v0sThisCollision = V0s.sliceBy(perCol, collision.globalIndex());
      for (auto& v0 : v0sThisCollision) { // loop over V0s

        if (!v0.has_mcParticle()) {
          continue;
        }
        auto v0mcparticle = v0.mcParticle();

        // Highest numerator of efficiency
        if (v0mcparticle.isPhysicalPrimary() && TMath::Abs(v0mcparticle.y()) <= rapidity) {
          if (v0mcparticle.pdgCode() == 310) {
            registry.fill(HIST("Reconstructed_MCRecoColl_K0Short"), v0mcparticle.pt()); // K0s
          }
          if (v0mcparticle.pdgCode() == 3122) {
            registry.fill(HIST("Reconstructed_MCRecoColl_Lambda"), v0mcparticle.pt()); // Lambda
          }
          if (v0mcparticle.pdgCode() == -3122) {
            registry.fill(HIST("Reconstructed_MCRecoColl_AntiLambda"), v0mcparticle.pt()); // AntiLambda
          }
        }

        int lPDG = 0;
        bool isprimary = false;
        if (TMath::Abs(v0mcparticle.pdgCode()) == 310 || TMath::Abs(v0mcparticle.pdgCode()) == 3122) {
          lPDG = v0mcparticle.pdgCode();
          isprimary = v0mcparticle.isPhysicalPrimary();
        }

        int posITSNhits = 0, negITSNhits = 0;
        for (unsigned int i = 0; i < 7; i++) {
          if (v0.posTrack_as<DauTracksMC>().itsClusterMap() & (1 << i)) {
            posITSNhits++;
          }
          if (v0.negTrack_as<DauTracksMC>().itsClusterMap() & (1 << i)) {
            negITSNhits++;
          }
        }

        float ctauLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * RecoDecay::getMassPDG(3122);
        float ctauAntiLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * RecoDecay::getMassPDG(-3122);
        float ctauK0s = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * RecoDecay::getMassPDG(310);

        if (v0.v0radius() > v0radius &&
            v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > v0cospa &&
            TMath::Abs(v0.posTrack_as<DauTracksMC>().eta()) < etadau &&
            TMath::Abs(v0.negTrack_as<DauTracksMC>().eta()) < etadau // &&
        ) {

          float cent = 0.;

          // Fill table
          myv0s(v0.globalIndex(), v0.pt(), v0.yLambda(), v0.yK0Short(),
                v0.mLambda(), v0.mAntiLambda(), v0.mK0Short(),
                v0.v0radius(), v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()),
                v0.dcapostopv(), v0.dcanegtopv(), v0.dcaV0daughters(),
                v0.posTrack_as<DauTracksMC>().eta(), v0.negTrack_as<DauTracksMC>().eta(),
                v0.posTrack_as<DauTracksMC>().phi(), v0.negTrack_as<DauTracksMC>().phi(),
                posITSNhits, negITSNhits, ctauLambda, ctauAntiLambda, ctauK0s,
                v0.negTrack_as<DauTracksMC>().tpcNSigmaPr(), v0.posTrack_as<DauTracksMC>().tpcNSigmaPr(),
                v0.negTrack_as<DauTracksMC>().tpcNSigmaPi(), v0.posTrack_as<DauTracksMC>().tpcNSigmaPi(),
                v0.negTrack_as<DauTracksMC>().tofNSigmaPr(), v0.posTrack_as<DauTracksMC>().tofNSigmaPr(),
                v0.negTrack_as<DauTracksMC>().tofNSigmaPi(), v0.posTrack_as<DauTracksMC>().tofNSigmaPi(),
                v0.posTrack_as<DauTracksMC>().hasTOF(), v0.negTrack_as<DauTracksMC>().hasTOF(), lPDG, isprimary,
                cent, cent);
        }
      }

      // Generated particles
      const auto particlesInCollision = mcParticles.sliceByCached(aod::mcparticle::mcCollisionId, collision.mcCollision().globalIndex(), cache1);

      for (auto& mcParticle : particlesInCollision) {
        if (std::abs(mcParticle.y()) > rapidity) {
          continue;
        }
        if (!mcParticle.isPhysicalPrimary()) {
          continue;
        }
        if (mcParticle.pdgCode() == 310)
          registry.fill(HIST("Generated_MCRecoColl_K0Short"), mcParticle.pt()); // K0s
        if (mcParticle.pdgCode() == 3122)
          registry.fill(HIST("Generated_MCRecoColl_Lambda"), mcParticle.pt()); // Lambda
        if (mcParticle.pdgCode() == -3122)
          registry.fill(HIST("Generated_MCRecoColl_AntiLambda"), mcParticle.pt()); // AntiLambda
      }
    }

    for (const auto& mccollision : mcCollisions) {

      if (TMath::Abs(mccollision.posZ()) > cutzvertex) {
        continue;
      }

      bool isFT0A = false;
      bool isFT0C = false;

      const auto particlesInMCCollision = mcParticles.sliceByCached(aod::mcparticle::mcCollisionId, mccollision.globalIndex(), cache2);

      for (auto& mcParticle : particlesInMCCollision) {

        if (!mcParticle.isPhysicalPrimary()) {
          continue;
        }

        if (std::abs(mcParticle.pdgCode()) == 211) {
          if (mcParticle.eta() <= -2.3 && mcParticle.eta() >= -3.4) {
            isFT0C = true;
          }
          if (mcParticle.eta() <= 5.0 && mcParticle.eta() >= 3.8) {
            isFT0A = true;
          }
        }

        if (std::abs(mcParticle.y()) > rapidity) {
          continue;
        }

        if (mcParticle.pdgCode() == 310)
          registry.fill(HIST("Generated_MCAllColl_K0Short"), mcParticle.pt()); // K0s
        if (mcParticle.pdgCode() == 3122)
          registry.fill(HIST("Generated_MCAllColl_Lambda"), mcParticle.pt()); // Lambda
        if (mcParticle.pdgCode() == -3122)
          registry.fill(HIST("Generated_MCAllColl_AntiLambda"), mcParticle.pt()); // AntiLambda
      }

      if (isFT0A && isFT0C) {
        registry.fill(HIST("hNEventsMC_AllColl"), 0.5);
      }
    }
  }
  PROCESS_SWITCH(v0qaanalysis, processMC, "Process MC", true);
};

struct myV0s {

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    registry.add("hMassLambda", "hMassLambda", {HistType::kTH1F, {{200, 1.016f, 1.216f}}});
    registry.add("hPt", "hPt", {HistType::kTH1F, {{100, 0.0f, 10.0f}}});
    registry.add("hMassVsPtLambda", "hMassVsPtLambda", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {200, 1.016f, 1.216f}}});
    registry.add("hMassAntiLambda", "hMassAntiLambda", {HistType::kTH1F, {{200, 1.016f, 1.216f}}});
    registry.add("hMassVsPtAntiLambda", "hMassVsPtAntiLambda", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {200, 1.016f, 1.216f}}});
    registry.add("hMassK0Short", "hMassK0Short", {HistType::kTH1F, {{200, 0.4f, 0.6f}}});
    registry.add("hMassVsPtK0Short", "hMassVsPtK0Short", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {200, 0.4f, 0.6f}}});
    registry.add("V0Radius", "V0Radius", {HistType::kTH1D, {{100, 0.0f, 20.0f}}});
    registry.add("CosPA", "CosPA", {HistType::kTH1F, {{100, 0.9f, 1.0f}}});
    registry.add("V0DCANegToPV", "V0DCANegToPV", {HistType::kTH1F, {{100, -1.0f, 1.0f}}});
    registry.add("V0DCAPosToPV", "V0DCAPosToPV", {HistType::kTH1F, {{100, -1.0f, 1.0f}}});
    registry.add("V0DCAV0Daughters", "V0DCAV0Daughters", {HistType::kTH1F, {{55, 0.0f, 2.20f}}});
    registry.add("CtauK0s", "CtauK0s", {HistType::kTH1F, {{150, 0.0f, 30.0f}}});
    registry.add("CtauLambda", "CtauLambda", {HistType::kTH1F, {{200, 0.0f, 40.0f}}});
    registry.add("CtauAntiLambda", "CtauAntiLambda", {HistType::kTH1F, {{200, 0.0f, 40.0f}}});
    registry.add("TPCNSigmaPosPi", "TPCNSigmaPosPi", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
    registry.add("TPCNSigmaNegPi", "TPCNSigmaNegPi", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
    registry.add("TPCNSigmaPosPr", "TPCNSigmaPosPr", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
    registry.add("TPCNSigmaNegPr", "TPCNSigmaNegPr", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
    registry.add("PosITSHits", "PosITSHits", {HistType::kTH1F, {{8, -0.5f, 7.5f}}});
    registry.add("NegITSHits", "NegITSHits", {HistType::kTH1F, {{8, -0.5f, 7.5f}}});
  }

  void process(aod::MyV0Candidates const& myv0s)
  {
    for (auto& candidate : myv0s) {

      registry.fill(HIST("hMassLambda"), candidate.masslambda());
      registry.fill(HIST("hPt"), candidate.v0pt());
      registry.fill(HIST("hMassVsPtLambda"), candidate.v0pt(), candidate.masslambda());
      registry.fill(HIST("hMassAntiLambda"), candidate.massantilambda());
      registry.fill(HIST("hMassVsPtAntiLambda"), candidate.v0pt(), candidate.massantilambda());
      registry.fill(HIST("hMassK0Short"), candidate.massk0short());
      registry.fill(HIST("hMassVsPtK0Short"), candidate.v0pt(), candidate.massk0short());
      registry.fill(HIST("V0Radius"), candidate.v0radius());
      registry.fill(HIST("CosPA"), candidate.v0cospa());
      registry.fill(HIST("V0DCANegToPV"), candidate.v0dcanegtopv());
      registry.fill(HIST("V0DCAPosToPV"), candidate.v0dcapostopv());
      registry.fill(HIST("V0DCAV0Daughters"), candidate.v0dcav0daughters());
      registry.fill(HIST("CtauK0s"), candidate.ctauk0short());
      registry.fill(HIST("CtauLambda"), candidate.ctaulambda());
      registry.fill(HIST("CtauAntiLambda"), candidate.ctauantilambda());
      registry.fill(HIST("TPCNSigmaPosPi"), candidate.ntpcsigmapospi());
      registry.fill(HIST("TPCNSigmaNegPi"), candidate.ntpcsigmanegpi());
      registry.fill(HIST("TPCNSigmaPosPr"), candidate.ntpcsigmapospr());
      registry.fill(HIST("TPCNSigmaNegPr"), candidate.ntpcsigmanegpr());
      registry.fill(HIST("PosITSHits"), candidate.v0positshits());
      registry.fill(HIST("NegITSHits"), candidate.v0negitshits());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<v0qaanalysis>(cfgc, TaskName{"lf-v0qaanalysis"}),
    adaptAnalysisTask<myV0s>(cfgc, TaskName{"lf-myv0s"})};
}
