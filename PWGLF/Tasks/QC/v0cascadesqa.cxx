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
/// \brief QA task for V0s and Cascades
///
/// In case of questions please write to:
/// \author Aimeric Landou (aimeric.landou@cern.ch)
/// \author Chiara De Martin (chiara.de.martin@cern.ch)
/// \author Francesca Ercolessi (francesca.ercolessi@cern.ch)

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"

#include <TFile.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <Math/Vector4D.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <cmath>
#include <array>
#include <cstdlib>
#include "Framework/ASoAHelpers.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra>;                       //, aod::McTrackLabels>; aod::pidTPCPi, aod::pidTPCPr
using MyTracksMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::McTrackLabels>; // aod::pidTPCPi, aod::pidTPCPr

struct v0cascadesQA {

  Configurable<bool> isMC{"isMC", false, "does the data have MC info"};
  Configurable<bool> sel8{"sel8", 0, "Apply sel8 event selection"};

  HistogramRegistry histos_eve{
    "histos-eve",
    {
      {"hEventCounter", "hEventCounter", {HistType::kTH1F, {{2, 0.0f, 2.0f}}}},
    },
    OutputObjHandlingPolicy::AnalysisObject,
    false,
    true};

  HistogramRegistry histos_V0{
    "histos-V0",
    {
      {"CosPA", "CosPA", {HistType::kTH1F, {{100, 0.9f, 1.0f}}}},
      {"V0Radius", "V0Radius", {HistType::kTH1D, {{100, 0.0f, 10.0f}}}},
      {"DecayLength", "DecayLength", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
      {"V0DCANegToPV", "V0DCANegToPV", {HistType::kTH1F, {{100, -1.0f, 1.0f}}}}, // this window captures but a small part of the distribution; plus the dcatopv can be negative ---- set to [-1;1] for comparison
      {"V0DCAPosToPV", "V0DCAPosToPV", {HistType::kTH1F, {{100, 0.0f, 1.0f}}}},  // this window captures but a small part of the distribution; plus the dcatopv can be negative
      {"V0DCAV0Daughters", "V0DCAV0Daughters", {HistType::kTH1F, {{55, 0.0f, 2.20f}}}},
      {"CtauK0s", "CtauK0s", {HistType::kTH1F, {{65, 0.0f, 13.0f}}}},
      {"CtauLambda", "CtauLambda", {HistType::kTH1F, {{200, 0.0f, 40.0f}}}},
      {"CtauAntiLambda", "CtauAntiLambda", {HistType::kTH1F, {{200, 0.0f, 40.0f}}}},
      {"DecayLengthK0s", "DecayLengthK0s", {HistType::kTH1F, {{100, 0.0f, 40.0f}}}},
      {"DecayLengthLambda", "DecayLengthLambda", {HistType::kTH1F, {{100, 0.0f, 80.0f}}}},
      {"DecayLengthAntiLambda", "DecayLengthAntiLambda", {HistType::kTH1F, {{100, 0.0f, 80.0f}}}},
      {"V0DCAV0ToPVK0S", "V0DCAV0ToPVK0S", {HistType::kTH1F, {{250, 0.0f, 0.25f}}}},
      {"V0DCAV0ToPVLambda", "V0DCAV0ToPVLambda", {HistType::kTH1F, {{250, 0.0f, 0.25f}}}},
      {"V0DCAV0ToPVAntiLambda", "V0DCAV0ToPVAntiLambda", {HistType::kTH1F, {{250, 0.0f, 0.25f}}}},

      {"ResponsePionFromLambda", "ResponsePionFromLambda", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {400, -20.f, 20.f}}}},
      {"ResponseProtonFromLambda", "ResponseProtonFromLambda", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {400, -20.f, 20.f}}}},

      {"InvMassK0S", "InvMassK0S", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {200, 0.4f, 0.6f}}}},
      {"InvMassLambda", "InvMassLambda", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {200, 1.07f, 1.17f}}}},
      {"InvMassAntiLambda", "InvMassAntiLambda", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {200, 1.07f, 1.17f}}}},
      {"InvMassK0S_Radius", "InvMassK0S_Radius", {HistType::kTH2F, {{100, 0.0f, 50.0f}, {200, 0.4f, 0.6f}}}},
      {"InvMassLambda_Radius", "InvMassLambda_Radius", {HistType::kTH2F, {{100, 0.0f, 50.0f}, {200, 1.07f, 1.17f}}}},
      {"InvMassAntiLambda_Radius", "InvMassAntiLambda_Radius", {HistType::kTH2F, {{100, 0.0f, 50.0f}, {200, 1.07f, 1.17f}}}},
      {"InvMassK0S_EtaDaughters", "InvMassK0S_EtaDaughters", {HistType::kTH3F, {{100, -1.0f, 1.0f}, {100, -1.0f, 1.0f}, {200, 0.4f, 0.6f}}}},
      {"InvMassLambda_EtaDaughters", "InvMassLambda_EtaDaughters", {HistType::kTH3F, {{100, -1.0f, 1.0f}, {100, -1.0f, 1.0f}, {200, 1.07f, 1.17f}}}},
      {"InvMassAntiLambda_EtaDaughters", "InvMassAntiLambda_EtaDaughters", {HistType::kTH3F, {{100, -1.0f, 1.0f}, {100, -1.0f, 1.0f}, {200, 1.07f, 1.17f}}}},
      {"InvMassK0S_PhiDaughters", "InvMassK0S_PhiDaughters", {HistType::kTH3F, {{628, 0.0f, 2. * TMath::Pi()}, {628, 0.0f, 2. * TMath::Pi()}, {200, 0.4f, 0.6f}}}},
      {"InvMassLambda_PhiDaughters", "InvMassLambda_PhiDaughters", {HistType::kTH3F, {{628, 0.0f, 2. * TMath::Pi()}, {628, 0.0f, 2. * TMath::Pi()}, {200, 1.07f, 1.17f}}}},
      {"InvMassAntiLambda_PhiDaughters", "InvMassAntiLambda_PhiDaughters", {HistType::kTH3F, {{628, 0.0f, 2. * TMath::Pi()}, {628, 0.0f, 2. * TMath::Pi()}, {200, 1.07f, 1.17f}}}},
      {"InvMassLambda_Ctau", "InvMassLambda_Ctau", {HistType::kTH2F, {{200, 0.0f, 40.0f}, {200, 1.07f, 1.17f}}}},
      {"InvMassAntiLambda_Ctau", "InvMassAntiLambda_Ctau", {HistType::kTH2F, {{200, 0.0f, 40.0f}, {200, 1.07f, 1.17f}}}},
      {"InvMassK0S_ITSMapDaughters", "InvMassK0S_ITSMapDaughters", {HistType::kTH3F, {{8, -0.5f, 7.5f}, {8, -0.5f, 7.5f}, {200, 0.4f, 0.6f}}}},
      {"InvMassLambda_ITSMapDaughters", "InvMassLambda_ITSMapDaughters", {HistType::kTH3F, {{8, -0.5f, 7.5f}, {8, -0.5f, 7.5f}, {200, 1.07f, 1.17f}}}},
      {"InvMassAntiLambda_ITSMapDaughters", "InvMassAntiLambda_ITSMapDaughters", {HistType::kTH3F, {{8, -0.5f, 7.5f}, {8, -0.5f, 7.5f}, {200, 1.07f, 1.17f}}}},
    },
    OutputObjHandlingPolicy::AnalysisObject,
    false,
    true};

  HistogramRegistry histos_Casc{
    "histos-Casc",
    {
      {"QA_XinusCandidates", "QA_XinusCandidates", {HistType::kTH1F, {{10, 0.f, 10.f}}}},
      {"XiProgSelections", "XiProgSelections", {HistType::kTH2F, {{30, 0.5f, 30.5f}, {2, -2, 2}}}},
      {"OmegaProgSelections", "OmegaProgSelections", {HistType::kTH2F, {{30, 0.5f, 30.5f}, {2, -2, 2}}}},
      {"CascCosPA", "CascCosPA", {HistType::kTH2F, {{350, 0.65f, 1.0f}, {2, -2, 2}}}},
      {"V0CosPA", "V0CosPA", {HistType::kTH2F, {{250, 0.75f, 1.0f}, {2, -2, 2}}}},
      {"V0CosPAToXi", "V0CosPAToXi", {HistType::kTH2F, {{100, 0.9f, 1.0f}, {2, -2, 2}}}},
      {"CascDecayLength", "CascDecayLength", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {2, -2, 2}}}},
      {"CascDecayLengthXi", "CascDecayLengthXi", {HistType::kTH2F, {{200, 0.0f, 20.0f}, {2, -2, 2}}}},
      {"CascDecayLengthOmega", "CascDecayLengthOmega", {HistType::kTH2F, {{200, 0.0f, 20.0f}, {2, -2, 2}}}},
      {"CascRadius", "CascRadius", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {2, -2, 2}}}},
      {"V0Radius", "V0Radius", {HistType::kTH2D, {{100, 0.0f, 10.0f}, {2, -2, 2}}}}, // asked for D instead of F
      {"CascyXi", "CascyXi", {HistType::kTH2F, {{200, -2.0f, 2.0f}, {2, -2, 2}}}},
      {"CascyOmega", "CascyOmega", {HistType::kTH2F, {{200, -2.0f, 2.0f}, {2, -2, 2}}}},
      {"CascCtauXi", "CascCtauXi", {HistType::kTH2F, {{100, 0.0f, 100.0f}, {2, -2, 2}}}},
      {"CascCtauOmega", "CascCtauOmega", {HistType::kTH2F, {{100, 0.0f, 100.0f}, {2, -2, 2}}}},
      {"V0Ctau", "V0Ctau", {HistType::kTH2F, {{100, 0.0f, 100.0f}, {2, -2, 2}}}},
      {"CascPt", "CascPt", {HistType::kTH2F, {{100, 0.0f, 25.0f}, {2, -2, 2}}}},
      {"DcaV0Daughters", "DcaV0Daughters", {HistType::kTH2F, {{110, 0.0f, 2.2f}, {2, -2, 2}}}},
      {"DcaCascDaughters", "DcaCascDaughters", {HistType::kTH2F, {{110, 0.0f, 2.2f}, {2, -2, 2}}}},
      {"DcaV0ToPV", "DcaV0ToPV", {HistType::kTH2F, {{2000, 0.0f, 2.f}, {2, -2, 2}}}},
      {"DcaBachToPV", "DcaBachToPV", {HistType::kTH2F, {{40, 0.0f, 0.2f}, {2, -2, 2}}}},
      {"DcaPosToPV", "DcaPosToPV", {HistType::kTH2F, {{40, 0.0f, 0.2f}, {2, -2, 2}}}},
      {"DcaNegToPV", "DcaNegToPV", {HistType::kTH2F, {{40, 0.0f, 0.2f}, {2, -2, 2}}}},
      {"InvMassLambdaDaughter", "InvMassLambdaDaughter", {HistType::kTH2F, {{100, 1.1f, 1.13f}, {2, -2, 2}}}},
      {"InvMassXiPlus", "InvMassXiPlus", {HistType::kTH2F, {{100, 0.f, 10.f}, {80, 1.28f, 1.36f}}}},
      {"InvMassXiMinus", "InvMassXiMinus", {HistType::kTH2F, {{100, 0.f, 10.f}, {80, 1.28f, 1.36f}}}},
      {"InvMassXiPlus_Radius", "InvMassXiPlus_Radius", {HistType::kTH2F, {{100, 0.f, 50.f}, {80, 1.28f, 1.36f}}}},
      {"InvMassXiMinus_Radius", "InvMassXiMinus_Radius", {HistType::kTH2F, {{100, 0.f, 50.f}, {80, 1.28f, 1.36f}}}},
      {"InvMassOmegaPlus", "InvMassOmegaPlus", {HistType::kTH2F, {{100, 0.f, 10.f}, {80, 1.63f, 1.71f}}}},
      {"InvMassOmegaMinus", "InvMassOmegaMinus", {HistType::kTH2F, {{100, 0.f, 10.f}, {80, 1.63f, 1.71f}}}},
    },
    OutputObjHandlingPolicy::AnalysisObject,
    false,
    true};

  void init(InitContext const&)
  {
    if (isMC) {
      histos_eve.add("GeneratedParticles", "GeneratedParticles", {HistType::kTH3F, {{14, 0.0f, 14.0f}, {100, 0, 10}, {100, 0.f, 50.f}}});

      histos_V0.add("InvMassK0STrue", "InvMassK0STrue", {HistType::kTH3F, {{100, 0.0f, 10.0f}, {100, 0.f, 50.f}, {200, 0.4f, 0.6f}}});
      histos_V0.add("InvMassLambdaTrue", "InvMassLambdaTrue", {HistType::kTH3F, {{100, 0.0f, 10.0f}, {100, 0.f, 50.f}, {200, 1.07f, 1.17f}}});
      histos_V0.add("InvMassAntiLambdaTrue", "InvMassAntiLambdaTrue", {HistType::kTH3F, {{100, 0.0f, 10.0f}, {100, 0.f, 50.f}, {200, 1.07f, 1.17f}}});

      histos_Casc.add("InvMassXiPlusTrue", "InvMassXiPlusTrue", {HistType::kTH3F, {{100, 0.f, 10.f}, {100, 0.f, 50.f}, {80, 1.28f, 1.36f}}});
      histos_Casc.add("InvMassXiMinusTrue", "InvMassXiMinusTrue", {HistType::kTH3F, {{100, 0.f, 10.f}, {100, 0.f, 50.f}, {80, 1.28f, 1.36f}}});
      histos_Casc.add("InvMassOmegaPlusTrue", "InvMassOmegaPlusTrue", {HistType::kTH3F, {{100, 0.f, 10.f}, {100, 0.f, 50.f}, {80, 1.63f, 1.71f}}});
      histos_Casc.add("InvMassOmegaMinusTrue", "InvMassOmegaMinusTrue", {HistType::kTH3F, {{100, 0.f, 10.f}, {100, 0.f, 50.f}, {80, 1.63f, 1.71f}}});
    }
  }

  ///////////////////////////////////////////////////
  ////////// Collisions QA - reconstructed //////////
  ///////////////////////////////////////////////////

  void processReconstructedEvent(soa::Join<aod::Collisions, aod::EvSels>::iterator const& Collision)
  {
    histos_eve.fill(HIST("hEventCounter"), 0.5);

    if (sel8 && !Collision.sel8()) {
      return;
    }

    histos_eve.fill(HIST("hEventCounter"), 1.5);
  }
  PROCESS_SWITCH(v0cascadesQA, processReconstructedEvent, "Process reconstructed level Event", true);

  ///////////////////////////////////////
  ////////// Collision QA - MC //////////
  ///////////////////////////////////////

  void processMcEvent(soa::Join<aod::McCollisions, aod::EvSels>::iterator const& mcCollision, aod::McParticles const& mcParticles)
  {
    if (sel8 && !mcCollision.sel8()) {
      return;
    }

    double posx = mcCollision.posX();
    double posy = mcCollision.posY();

    for (auto& mcparticle : mcParticles) {

      if (!mcparticle.has_daughters()) {
        continue;
      }

      double vx = 0;
      double vy = 0;
      for (auto& mcparticleDaughter0 : mcparticle.daughters_as<aod::McParticles>()) {
        vx = mcparticleDaughter0.vx() - posx;
        vy = mcparticleDaughter0.vy() - posy;
        if (vx != 0 && vy != 0)
          break;
      }
      double R_Decay = TMath::Sqrt(vx * vx + vy * vy);

      if (mcparticle.isPhysicalPrimary() && TMath::Abs(mcparticle.y()) < V0_rapidity) {
        if (mcparticle.pdgCode() == 310)
          histos_eve.fill(HIST("GeneratedParticles"), 0.5, mcparticle.pt(), R_Decay); // K0s
        if (mcparticle.pdgCode() == 3122)
          histos_eve.fill(HIST("GeneratedParticles"), 2.5, mcparticle.pt(), R_Decay); // Lambda
        if (mcparticle.pdgCode() == -3122)
          histos_eve.fill(HIST("GeneratedParticles"), 4.5, mcparticle.pt(), R_Decay); // AntiLambda
      }
      if (mcparticle.isPhysicalPrimary() && TMath::Abs(mcparticle.y()) < Casc_rapidity) {
        if (mcparticle.pdgCode() == 3312)
          histos_eve.fill(HIST("GeneratedParticles"), 6.5, mcparticle.pt(), R_Decay); // Xi-
        if (mcparticle.pdgCode() == -3312)
          histos_eve.fill(HIST("GeneratedParticles"), 8.5, mcparticle.pt(), R_Decay); // Xi+
        if (mcparticle.pdgCode() == 3334)
          histos_eve.fill(HIST("GeneratedParticles"), 10.5, mcparticle.pt(), R_Decay); // Omega-
        if (mcparticle.pdgCode() == -3334)
          histos_eve.fill(HIST("GeneratedParticles"), 12.5, mcparticle.pt(), R_Decay); // Omega+

        // if (!IsParticleFromOutOfBunchPileupCollision){fill the 1.5, 3.5 etc}   AliPhysics analysis
      }
    }
  }
  PROCESS_SWITCH(v0cascadesQA, processMcEvent, "Process MC level Event", true);

  ////////////////////////////////////////////
  ////////// V0 QA - Reconstructed ///////////
  ////////////////////////////////////////////

  Configurable<float> V0_rapidity{"V0_rapidity", 0.5, "rapidity"};
  Configurable<double> V0_cosPA{"V0_cosPA", 0.995, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<float> V0_dcav0dau{"V0_dcav0dau", 1, "DCA V0 Daughters"};
  Configurable<float> V0_dcapostopv{"V0_dcapostopv", 0.1, "DCA Pos To PV"};
  Configurable<float> V0_dcanegtopv{"V0_dcanegtopv", 0.1, "DCA Neg To PV"};
  Configurable<float> V0_radius{"V0_radius", 5, "v0radius"};

  static constexpr float defaultLifetimeCuts[1][2] = {{25., 20.}};
  Configurable<LabeledArray<float>> lifetimecut{"lifetimecut", {defaultLifetimeCuts[0], 2, {"lifetimecutLambda", "lifetimecutK0S"}}, "lifetimecut"};
  using DaughterTracks = soa::Join<aod::Tracks, aod::TracksCov, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTOFPi, aod::pidTPCPi, aod::pidTOFPr, aod::pidTPCPr>;

  void processReconstructedV0(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, MyTracks const& tracks, aod::V0Datas const& fullV0s, DaughterTracks& dtracks)
  {
    if (sel8 && !collision.sel8()) {
      return;
    }

    for (auto& v0 : fullV0s) {

      auto posdau = v0.posTrack_as<DaughterTracks>();
      auto negdau = v0.negTrack_as<DaughterTracks>();

      Int_t posITSNhits = 0, negITSNhits = 0;
      for (unsigned int i = 0; i < 7; i++) {
        if (posdau.itsClusterMap() & (1 << i)) {
          posITSNhits++;
        }
        if (negdau.itsClusterMap() & (1 << i)) {
          negITSNhits++;
        }
      }

      histos_V0.fill(HIST("CosPA"), v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
      histos_V0.fill(HIST("V0Radius"), v0.v0radius());
      histos_V0.fill(HIST("V0DCANegToPV"), v0.dcanegtopv());
      histos_V0.fill(HIST("V0DCAPosToPV"), v0.dcapostopv());
      histos_V0.fill(HIST("V0DCAV0Daughters"), v0.dcaV0daughters());

      float decayLength = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * RecoDecay::sqrtSumOfSquares(v0.px(), v0.py(), v0.pz());
      histos_V0.fill(HIST("DecayLength"), decayLength);

      float CtauLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * RecoDecay::getMassPDG(kLambda0);
      float CtauK0s = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * RecoDecay::getMassPDG(kK0Short);

      if (v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > V0_cosPA) {
        if (v0.v0radius() > V0_radius && v0.dcaV0daughters() < V0_dcav0dau && TMath::Abs(v0.dcapostopv()) > V0_dcapostopv && TMath::Abs(v0.dcanegtopv()) > V0_dcanegtopv) {
          if (TMath::Abs(v0.yK0Short()) < V0_rapidity && CtauK0s < lifetimecut->get("lifetimecutK0S")) { // that s what we have in lambdakzeroanalysis ; discuss with nicolo chiara if we want to have more like in th aliphysics macro or none like in David's O2 QA
            histos_V0.fill(HIST("CtauK0s"), CtauK0s);
            histos_V0.fill(HIST("DecayLengthK0s"), decayLength);
            histos_V0.fill(HIST("InvMassK0S"), v0.pt(), v0.mK0Short());
            histos_V0.fill(HIST("InvMassK0S_Radius"), v0.v0radius(), v0.mK0Short());
            histos_V0.fill(HIST("InvMassK0S_EtaDaughters"), posdau.eta(), negdau.eta(), v0.mK0Short());
            histos_V0.fill(HIST("InvMassK0S_PhiDaughters"), posdau.phi(), negdau.phi(), v0.mK0Short());
            histos_V0.fill(HIST("V0DCAV0ToPVK0S"), v0.dcav0topv(collision.posX(), collision.posY(), collision.posZ()));
            histos_V0.fill(HIST("InvMassK0S_ITSMapDaughters"), posITSNhits, negITSNhits, v0.mK0Short());
          }

          if (TMath::Abs(v0.yLambda()) < V0_rapidity && CtauLambda < lifetimecut->get("lifetimecutLambda")) {
            histos_V0.fill(HIST("DecayLengthLambda"), decayLength);
            histos_V0.fill(HIST("CtauLambda"), CtauLambda);
            histos_V0.fill(HIST("InvMassLambda"), v0.pt(), v0.mLambda());
            histos_V0.fill(HIST("InvMassLambda_Radius"), v0.v0radius(), v0.mLambda());
            histos_V0.fill(HIST("InvMassLambda_Ctau"), CtauLambda, v0.mLambda());
            histos_V0.fill(HIST("InvMassLambda_EtaDaughters"), posdau.eta(), negdau.eta(), v0.mLambda());
            histos_V0.fill(HIST("InvMassLambda_PhiDaughters"), posdau.phi(), negdau.phi(), v0.mLambda());
            histos_V0.fill(HIST("V0DCAV0ToPVLambda"), v0.dcav0topv(collision.posX(), collision.posY(), collision.posZ()));
            if (v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > 0.999 && v0.dcaV0daughters() < 1 && TMath::Abs(v0.mK0Short() - RecoDecay::getMassPDG(310)) > 0.012 && TMath::Abs(v0.mAntiLambda() - RecoDecay::getMassPDG(3122)) > 0.08 && TMath::Abs(v0.mLambda() - RecoDecay::getMassPDG(3122)) < 0.002) {
              histos_V0.fill(HIST("ResponsePionFromLambda"), v0.pt(), negdau.tpcNSigmaPi());
              histos_V0.fill(HIST("ResponseProtonFromLambda"), v0.pt(), posdau.tpcNSigmaPr());
            }
            histos_V0.fill(HIST("InvMassLambda_ITSMapDaughters"), posITSNhits, negITSNhits, v0.mLambda());
          }

          if (TMath::Abs(v0.yLambda()) < V0_rapidity && CtauLambda < lifetimecut->get("lifetimecutLambda")) {
            histos_V0.fill(HIST("DecayLengthAntiLambda"), decayLength);
            histos_V0.fill(HIST("CtauAntiLambda"), CtauLambda);
            histos_V0.fill(HIST("InvMassAntiLambda"), v0.pt(), v0.mAntiLambda());
            histos_V0.fill(HIST("InvMassAntiLambda_Radius"), v0.v0radius(), v0.mAntiLambda());
            histos_V0.fill(HIST("InvMassAntiLambda_Ctau"), CtauLambda, v0.mAntiLambda());
            histos_V0.fill(HIST("InvMassAntiLambda_EtaDaughters"), posdau.eta(), negdau.eta(), v0.mAntiLambda());
            histos_V0.fill(HIST("InvMassAntiLambda_PhiDaughters"), posdau.phi(), negdau.phi(), v0.mAntiLambda());
            histos_V0.fill(HIST("V0DCAV0ToPVAntiLambda"), v0.dcav0topv(collision.posX(), collision.posY(), collision.posZ()));
            histos_V0.fill(HIST("InvMassAntiLambda_ITSMapDaughters"), posITSNhits, negITSNhits, v0.mAntiLambda());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(v0cascadesQA, processReconstructedV0, "Process reconstructed level V0s", true);

  ////////////////////////////////
  ////////// V0 QA - MC //////////
  ////////////////////////////////

  void processMcV0(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, MyTracksMC const& tracks, aod::V0Datas const& fullV0s, aod::McParticles const& mcParticles)
  {
    if (sel8 && !collision.sel8()) {
      return;
    }

    for (auto& v0 : fullV0s) {

      float CtauLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * RecoDecay::getMassPDG(kLambda0);
      float CtauK0s = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * RecoDecay::getMassPDG(kK0Short);

      if (v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > V0_cosPA &&
          v0.v0radius() > V0_radius &&
          v0.dcaV0daughters() < V0_dcav0dau &&
          TMath::Abs(v0.dcapostopv()) > V0_dcapostopv &&
          TMath::Abs(v0.dcanegtopv()) > V0_dcanegtopv) {

        auto reconegtrack = v0.negTrack_as<MyTracksMC>();
        auto recopostrack = v0.posTrack_as<MyTracksMC>();
        if (!reconegtrack.has_mcParticle() || !recopostrack.has_mcParticle()) {
          continue;
        }

        auto mcnegtrack = reconegtrack.mcParticle_as<aod::McParticles>();
        auto mcpostrack = recopostrack.mcParticle_as<aod::McParticles>();
        if (!mcnegtrack.has_mothers() || !mcpostrack.has_mothers()) {
          continue;
        }

        for (auto& particleMotherOfNeg : mcnegtrack.mothers_as<aod::McParticles>()) {
          for (auto& particleMotherOfPos : mcpostrack.mothers_as<aod::McParticles>()) {

            bool MomIsPrimary = particleMotherOfNeg.isPhysicalPrimary();

            bool isK0sV0 = MomIsPrimary &&
                           particleMotherOfNeg == particleMotherOfPos &&
                           particleMotherOfNeg.pdgCode() == 310 &&
                           mcnegtrack.pdgCode() == -211 &&
                           mcpostrack.pdgCode() == +211;

            bool isLambdaV0 = MomIsPrimary &&
                              particleMotherOfNeg == particleMotherOfPos &&
                              particleMotherOfNeg.pdgCode() == 3122 &&
                              mcnegtrack.pdgCode() == -211 &&
                              mcpostrack.pdgCode() == 2212;

            bool isAntiLambdaV0 = MomIsPrimary &&
                                  particleMotherOfNeg == particleMotherOfPos &&
                                  particleMotherOfNeg.pdgCode() == -3122 &&
                                  mcpostrack.pdgCode() == 211 &&
                                  mcnegtrack.pdgCode() == -2212;

            if (isK0sV0) {
              if (TMath::Abs(v0.yK0Short()) < V0_rapidity && CtauK0s < lifetimecut->get("lifetimecutK0S")) {
                histos_V0.fill(HIST("InvMassK0STrue"), v0.pt(), v0.v0radius(), v0.mK0Short());
              }
            }
            if (isLambdaV0) {
              if (TMath::Abs(v0.yLambda()) < V0_rapidity && CtauLambda < lifetimecut->get("lifetimecutLambda")) {
                histos_V0.fill(HIST("InvMassLambdaTrue"), v0.pt(), v0.v0radius(), v0.mLambda());
              }
            }
            if (isAntiLambdaV0) {
              if (TMath::Abs(v0.yLambda()) < V0_rapidity && CtauLambda < lifetimecut->get("lifetimecutLambda")) {
                histos_V0.fill(HIST("InvMassAntiLambdaTrue"), v0.pt(), v0.v0radius(), v0.mAntiLambda());
              }
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(v0cascadesQA, processMcV0, "Process MC level V0s", false);

  //////////////////////////////////////
  ///// Cascade QA - Reconstructed /////
  //////////////////////////////////////

  Configurable<float> Casc_rapidity{"Casc_rapidity", 0.5, "rapidity"};
  Configurable<double> Casc_v0cospa{"Casc_V0cospa", 0.98, "V0 CosPA"};               // double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<double> Casc_casccospa{"Casc_casccospa", 0.98, "Cascade CosPA"};      // AliAnalysisTaskStrAODqa: 0.9992      //double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<float> Casc_dcav0dau{"Casc_dcav0dau", 1.0, "DCA V0 Daughters"};       // AliAnalysisTaskStrAODqa: 1. different scale
  Configurable<float> Casc_dcacascdau{"Casc_dcacascdau", 0.6, "DCA Casc Daughters"}; // AliAnalysisTaskStrAODqa: 0.3 different scale
  Configurable<float> Casc_dcav0topv{"Casc_dcav0topv", 0.1, "DCA Pos To PV"};        // AliAnalysisTaskStrAODqa: 0.15 different scale
  Configurable<float> Casc_dcabachtopv{"Casc_dcabachtopv", .1, "DCA Bach To PV"};    // AliAnalysisTaskStrAODqa: 0.17 different scale
  Configurable<float> Casc_dcapostopv{"Casc_dcapostopv", 0.1, "DCA V0 To PV"};       // AliAnalysisTaskStrAODqa:    if( fCasc_charge>0 &&(fCasc_DcaPosToPV < 0.3 || fCasc_DcaNegToPV < 0.11)) return kFALSE;  different scale
  Configurable<float> Casc_dcanegtopv{"Casc_dcanegtopv", 0.1, "DCA Neg To PV"};      // AliAnalysisTaskStrAODqa:    if( fCasc_charge<0 &&(fCasc_DcaPosToPV < 0.11 || fCasc_DcaNegToPV < 0.3)) return kFALSE;  different scale
  Configurable<float> Casc_v0radius{"Casc_v0radius", 0.9, "v0 radius"};              // AliAnalysisTaskStrAODqa: 5.
  Configurable<float> Casc_cascradius{"Casc_cascradius", 1.0, "cascade radius"};     // AliAnalysisTaskStrAODqa: 1.

  // additional AliAnalysisTaskStrAODqa.cxx cuts not present here
  //  if( fCasc_LeastCRaws<70 ) return kFALSE; //assume LeastCRows

  // if( fCasc_V0CosPAToXi<0.99 ) return kFALSE;

  // if( TMath::Abs(fCasc_InvMassLambda-1.115683)>0.005) return kFALSE;
  // if( (part > 4) && TMath::Abs(fCasc_InvMassXi-1.32171)<0.003) return kFALSE;
  // if( (part<5) && (fCasc_CascCtauXi> (4.91*3)) ) return kFALSE;   //4.91 is the ctau of xi in cm
  // if( (part>=5) && (fCasc_CascCtauOmega > (2.461*3)) ) return kFALSE;   //2.461 is the ctau of om in cm

  // if( (part==3) && (TMath::Abs(fCasc_NSigPosPion)>3 || TMath::Abs(fCasc_NSigNegProton)>3 || TMath::Abs(fCasc_NSigBacPion)>3) ) return kFALSE;
  // if( (part==4) && (TMath::Abs(fCasc_NSigNegPion)>3 || TMath::Abs(fCasc_NSigPosProton)>3 || TMath::Abs(fCasc_NSigBacPion)>3) ) return kFALSE;
  // if( (part==5) && (TMath::Abs(fCasc_NSigPosPion)>3 || TMath::Abs(fCasc_NSigNegProton)>3 || TMath::Abs(fCasc_NSigBacKaon)>3) ) return kFALSE;
  // if( (part==6) && (TMath::Abs(fCasc_NSigNegPion)>3 || TMath::Abs(fCasc_NSigPosProton)>3 || TMath::Abs(fCasc_NSigBacKaon)>3) ) return kFALSE;

  void processReconstructedCascade(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::CascDataExt const& Cascades, aod::V0Datas const& fullV0s)
  {
    if (sel8 && !collision.sel8()) {
      return;
    }

    for (auto& casc : Cascades) {
      // histos_Casc.fill(HIST("XiProgSelections"), );
      // histos_Casc.fill(HIST("OmegaProgSelections"), );
      histos_Casc.fill(HIST("CascCosPA"), casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()), casc.sign());
      histos_Casc.fill(HIST("V0CosPA"), casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()), casc.sign());

      // double v0cospatoxi = RecoDecay::CPA(array{casc.x(), casc.y(), casc.z()}, array{casc.xlambda(), casc.ylambda(), casc.zlambda()}, array{v0.px(), v0.py(), v0.pz()});
      double v0cospatoxi = RecoDecay::cpa(array{casc.x(), casc.y(), casc.z()}, array{casc.xlambda(), casc.ylambda(), casc.zlambda()}, array{casc.pxpos() + casc.pxneg(), casc.pypos() + casc.pyneg(), casc.pzpos() + casc.pzneg()});

      histos_Casc.fill(HIST("V0CosPAToXi"), v0cospatoxi, casc.sign());
      histos_Casc.fill(HIST("CascRadius"), casc.cascradius(), casc.sign());
      histos_Casc.fill(HIST("V0Radius"), casc.v0radius(), casc.sign());
      histos_Casc.fill(HIST("CascyXi"), casc.yXi(), casc.sign());
      histos_Casc.fill(HIST("CascyOmega"), casc.yOmega(), casc.sign());

      float cascDecayLength = std::sqrt(std::pow(casc.x() - collision.posX(), 2) + std::pow(casc.y() - collision.posY(), 2) + std::pow(casc.z() - collision.posZ(), 2));
      histos_Casc.fill(HIST("CascDecayLength"), cascDecayLength, casc.sign());
      histos_Casc.fill(HIST("CascDecayLengthXi"), cascDecayLength, casc.sign());
      histos_Casc.fill(HIST("CascDecayLengthOmega"), cascDecayLength, casc.sign());

      float cascTotalMomentum = RecoDecay::sqrtSumOfSquares(casc.px(), casc.py(), casc.pz());
      float CtauXi = cascDecayLength / (cascTotalMomentum + 1E-10) * RecoDecay::getMassPDG(3322); // see O2Physics/Common/Core/MC.h for codes and names accepted
      float CtauOmega = cascDecayLength / (cascTotalMomentum + 1E-10) * RecoDecay::getMassPDG(3334);

      float v0TotalMomentum = RecoDecay::sqrtSumOfSquares(casc.pxpos() + casc.pxneg(), casc.pypos() + casc.pyneg(), casc.pzpos() + casc.pzneg());
      float v0DecayLength = std::sqrt(std::pow(casc.xlambda() - casc.x(), 2) + std::pow(casc.ylambda() - casc.y(), 2) + std::pow(casc.zlambda() - casc.z(), 2));
      float CtauV0 = v0DecayLength / (v0TotalMomentum + 1E-10) * RecoDecay::getMassPDG(kLambda0);

      histos_Casc.fill(HIST("CascCtauXi"), CtauXi, casc.sign());
      histos_Casc.fill(HIST("CascCtauOmega"), CtauOmega, casc.sign());
      histos_Casc.fill(HIST("V0Ctau"), CtauV0, casc.sign());
      histos_Casc.fill(HIST("CascPt"), casc.pt(), casc.sign());
      histos_Casc.fill(HIST("DcaV0Daughters"), casc.dcaV0daughters(), casc.sign());
      histos_Casc.fill(HIST("DcaCascDaughters"), casc.dcacascdaughters(), casc.sign());
      histos_Casc.fill(HIST("DcaV0ToPV"), casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()), casc.sign());
      histos_Casc.fill(HIST("DcaBachToPV"), casc.dcabachtopv(), casc.sign());
      histos_Casc.fill(HIST("DcaPosToPV"), casc.dcapostopv(), casc.sign());
      histos_Casc.fill(HIST("DcaNegToPV"), casc.dcanegtopv(), casc.sign());
      histos_Casc.fill(HIST("InvMassLambdaDaughter"), casc.mLambda(), casc.sign());

      if (casc.v0radius() > Casc_v0radius &&
          casc.cascradius() > Casc_cascradius &&
          casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > Casc_v0cospa &&
          casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) > Casc_casccospa &&
          TMath::Abs(casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ())) > Casc_dcav0topv &&
          TMath::Abs(casc.dcapostopv()) > Casc_dcapostopv && TMath::Abs(casc.dcanegtopv()) > Casc_dcanegtopv && TMath::Abs(casc.dcabachtopv()) > Casc_dcabachtopv &&
          casc.dcaV0daughters() < Casc_dcav0dau && casc.dcacascdaughters() < Casc_dcacascdau) {
        if (casc.sign() < 0) {
          if (TMath::Abs(casc.yXi()) < Casc_rapidity) {
            histos_Casc.fill(HIST("InvMassXiMinus"), casc.pt(), casc.mXi());
            histos_Casc.fill(HIST("InvMassXiMinus_Radius"), casc.cascradius(), casc.mXi());
          }
          if (TMath::Abs(casc.yOmega()) < Casc_rapidity) {
            histos_Casc.fill(HIST("InvMassOmegaMinus"), casc.pt(), casc.mOmega());
          }
        } else {
          if (TMath::Abs(casc.yXi()) < Casc_rapidity) {
            histos_Casc.fill(HIST("InvMassXiPlus"), casc.pt(), casc.mXi());
            histos_Casc.fill(HIST("InvMassXiPlus_Radius"), casc.cascradius(), casc.mXi());
          }
          if (TMath::Abs(casc.yOmega()) < Casc_rapidity) {
            histos_Casc.fill(HIST("InvMassOmegaPlus"), casc.pt(), casc.mOmega());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(v0cascadesQA, processReconstructedCascade, "Process reconstructed level Cascades", true);

  //////////////////////////////////////
  ////////// Cascade QA - MC ///////////
  //////////////////////////////////////

  void processMcCascade(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::CascDataExt const& Cascades, aod::V0sLinked const&, aod::V0Datas const& fullV0s, MyTracksMC const& tracks, aod::McParticles const& mcParticles)
  {
    if (sel8 && !collision.sel8()) {
      return;
    }

    for (auto& casc : Cascades) {

      histos_Casc.fill(HIST("QA_XinusCandidates"), 0.5);

      if (casc.v0radius() > Casc_v0radius &&
          casc.cascradius() > Casc_cascradius &&
          casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > Casc_v0cospa &&
          casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) > Casc_casccospa &&
          TMath::Abs(casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ())) > Casc_dcav0topv &&
          TMath::Abs(casc.dcapostopv()) > Casc_dcapostopv && TMath::Abs(casc.dcanegtopv()) > Casc_dcanegtopv && TMath::Abs(casc.dcabachtopv()) > Casc_dcabachtopv &&
          casc.dcaV0daughters() < Casc_dcav0dau && casc.dcacascdaughters() < Casc_dcacascdau) {

        histos_Casc.fill(HIST("QA_XinusCandidates"), 1.5);

        auto v0index = casc.v0_as<o2::aod::V0sLinked>();
        if (!(v0index.has_v0Data())) {
          continue; // skip those cascades for which V0 doesn't exist
        }
        auto v0 = v0index.v0Data(); // de-reference index to correct v0data in case it exists

        auto reconegtrack = v0.negTrack_as<MyTracksMC>();
        auto recopostrack = v0.posTrack_as<MyTracksMC>();
        auto recobachelor = casc.bachelor_as<MyTracksMC>();
        if (!reconegtrack.has_mcParticle() || !recopostrack.has_mcParticle() || !recobachelor.has_mcParticle()) {
          continue;
        }
        histos_Casc.fill(HIST("QA_XinusCandidates"), 2.5);

        auto mcnegtrack = reconegtrack.mcParticle_as<aod::McParticles>();
        auto mcpostrack = recopostrack.mcParticle_as<aod::McParticles>();
        auto bachelor = recobachelor.mcParticle_as<aod::McParticles>();
        if (!mcnegtrack.has_mothers() || !mcpostrack.has_mothers() || !bachelor.has_mothers()) {
          continue;
        }
        histos_Casc.fill(HIST("QA_XinusCandidates"), 3.5);

        for (auto& particleMotherOfBach : bachelor.mothers_as<aod::McParticles>()) {
          for (auto& particleMotherOfNeg : mcnegtrack.mothers_as<aod::McParticles>()) {
            for (auto& particleMotherOfPos : mcpostrack.mothers_as<aod::McParticles>()) {
              for (auto& particleMotherOfV0 : particleMotherOfNeg.mothers_as<aod::McParticles>()) {

                bool MomOfBachIsPrimary = particleMotherOfBach.isPhysicalPrimary();
                bool MomOfNegIsPrimary = particleMotherOfNeg.isPhysicalPrimary();
                bool MomOfPosIsPrimary = particleMotherOfPos.isPhysicalPrimary();

                bool isXiMinusCascade = MomOfBachIsPrimary && !(MomOfNegIsPrimary) && !(MomOfPosIsPrimary) &&
                                        particleMotherOfNeg == particleMotherOfPos &&
                                        particleMotherOfV0 == particleMotherOfBach &&
                                        particleMotherOfBach.pdgCode() == 3312 &&
                                        bachelor.pdgCode() == -211 &&
                                        particleMotherOfNeg.pdgCode() == 3122 &&
                                        mcnegtrack.pdgCode() == -211 &&
                                        mcpostrack.pdgCode() == 2212;

                bool isOmegaMinusCascade = MomOfBachIsPrimary && !(MomOfNegIsPrimary) && !(MomOfPosIsPrimary) &&
                                           particleMotherOfNeg == particleMotherOfPos &&
                                           particleMotherOfV0 == particleMotherOfBach &&
                                           particleMotherOfBach.pdgCode() == 3334 &&
                                           bachelor.pdgCode() == -321 &&
                                           particleMotherOfNeg.pdgCode() == 3122 &&
                                           mcnegtrack.pdgCode() == -211 &&
                                           mcpostrack.pdgCode() == 2212;

                bool isXiPlusCascade = MomOfBachIsPrimary && !(MomOfNegIsPrimary) && !(MomOfPosIsPrimary) &&
                                       particleMotherOfNeg == particleMotherOfPos &&
                                       particleMotherOfV0 == particleMotherOfBach &&
                                       particleMotherOfBach.pdgCode() == -3312 &&
                                       bachelor.pdgCode() == 211 &&
                                       particleMotherOfNeg.pdgCode() == -3122 &&
                                       mcnegtrack.pdgCode() == -2212 &&
                                       mcpostrack.pdgCode() == 211;

                bool isOmegaPlusCascade = MomOfBachIsPrimary && !(MomOfNegIsPrimary) && !(MomOfPosIsPrimary) &&
                                          particleMotherOfNeg == particleMotherOfPos &&
                                          particleMotherOfV0 == particleMotherOfBach &&
                                          particleMotherOfBach.pdgCode() == -3334 &&
                                          bachelor.pdgCode() == 321 &&
                                          particleMotherOfNeg.pdgCode() == -3122 &&
                                          mcnegtrack.pdgCode() == -2212 &&
                                          mcpostrack.pdgCode() == 211;

                if (isXiMinusCascade) {
                  histos_Casc.fill(HIST("QA_XinusCandidates"), 8.5);
                  if (TMath::Abs(casc.yXi()) < Casc_rapidity) {
                    histos_Casc.fill(HIST("InvMassXiMinusTrue"), casc.pt(), casc.cascradius(), casc.mXi());
                  }
                }

                if (isOmegaMinusCascade) {
                  if (TMath::Abs(casc.yOmega()) < Casc_rapidity) {
                    histos_Casc.fill(HIST("InvMassOmegaMinusTrue"), casc.pt(), casc.cascradius(), casc.mOmega());
                  }
                }

                if (isXiPlusCascade) {
                  if (TMath::Abs(casc.yXi()) < Casc_rapidity) {
                    histos_Casc.fill(HIST("InvMassXiPlusTrue"), casc.pt(), casc.cascradius(), casc.mXi());
                  }
                }

                if (isOmegaPlusCascade) {
                  if (TMath::Abs(casc.yOmega()) < Casc_rapidity) {
                    histos_Casc.fill(HIST("InvMassOmegaPlusTrue"), casc.pt(), casc.cascradius(), casc.mOmega());
                  }
                }

                // QA section - XiMinus
                if (MomOfBachIsPrimary && !(MomOfNegIsPrimary) && !(MomOfPosIsPrimary)) {
                  histos_Casc.fill(HIST("QA_XinusCandidates"), 4.5);
                }
                if ((particleMotherOfNeg.pdgCode() == particleMotherOfPos.pdgCode())) {
                  histos_Casc.fill(HIST("QA_XinusCandidates"), 5.5);
                }
                if ((particleMotherOfV0 == particleMotherOfBach)) {
                  histos_Casc.fill(HIST("QA_XinusCandidates"), 6.5);
                }
                if (particleMotherOfBach.pdgCode() == 3312 &&
                    bachelor.pdgCode() == -211 &&
                    particleMotherOfNeg.pdgCode() == 3122 &&
                    mcnegtrack.pdgCode() == -211 &&
                    mcpostrack.pdgCode() == 2212) {
                  histos_Casc.fill(HIST("QA_XinusCandidates"), 7.5);
                }
              }
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(v0cascadesQA, processMcCascade, "Process MC level Cascades", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<v0cascadesQA>(cfgc)};
}
