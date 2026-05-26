// Copyright 2020-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \author Kotliarov Artem <artem.kotliarov@cern.ch>
/// \file trackJetSpectra.cxx
/// \brief track and jet analysis

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetSubtraction.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Multiplicity.h"

#include "CommonConstants/MathConstants.h"
#include "Framework/ASoA.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

#include "TRandom3.h"
#include <TH1.h>
#include <TMath.h>
#include <TString.h>
#include <TVector2.h>

#include <cmath>
#include <cstdint>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using FilteredColl = soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos>>::iterator;
using FilteredEventMultiplicity = soa::Filtered<soa::Join<aod::JetCollisions, aod::ZDCMults>>::iterator;

using FilteredJetsDetLevel = soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>>;
using FilteredJetsPartLevel = soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents>>;

using FilteredCollDetLevelGetWeight = soa::Filtered<soa::Join<aod::JetCollisionsMCD, aod::BkgChargedRhos, aod::JCollisionOutliers>>::iterator;

using FilteredCollPartLevel = soa::Filtered<soa::Join<aod::JetMcCollisions, aod::BkgChargedMcRhos, aod::JMcCollisionOutliers>>::iterator;
using FilteredParticles = soa::Filtered<aod::JetParticles>;
using FilteredJetsPartLevel = soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents>>;

using FilteredTracks = soa::Filtered<aod::JetTracks>;
using FilteredMatchedJetsDetLevel = soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>>;
using FilteredMatchedJetsPartLevel = soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>>;

struct TrackJetSpectra {

  // List of configurable parameters
  Configurable<std::string> evSel{"evSel", "sel8", "Choose event selection"};
  Configurable<std::string> trkSel{"trkSel", "globalTracks", "Set track selection"};
  Configurable<float> vertexZCut{"vertexZCut", 10., "Accepted z-vertex range"};

  Configurable<float> trkPtMin{"trkPtMin", 0.15, "Minimum pT of acceptanced tracks"};
  Configurable<float> trkPtMax{"trkPtMax", 100., "Maximum pT of acceptanced tracks"};

  Configurable<float> trkEtaCut{"trkEtaCut", 0.9, "Eta acceptance of TPC"};
  Configurable<float> jetR{"jetR", 0.4, "Jet cone radius"};

  Configurable<float> SignalTriggerMin{"STriggerMin", 15., "Minimal Signal Trigger Track pT"};
  Configurable<float> SignalTriggerMax{"STriggerMax", 50., "Maximal Signal Trigger Track pT"};
  Configurable<float> ReferenceTriggerMin{"RTriggerMin", 5., "Minimal Reference Trigger Track pT"};
  Configurable<float> ReferenceTriggerMax{"RTriggerMax", 7., "Maximal Reference Trigger Track pT"};
  Configurable<float> SigToRefFraction{"Fraction", 0.1, ""};
  Configurable<float> DPhiCut{"DPhiCut", 0.6, ""};
  // List of configurable parameters for histograms
  Configurable<uint16_t> histJetPt{"histJetPt", 150, "Maximum value of jet pT shown in histograms"};

  Configurable<float> meanFT0A{"meanFT0A", -1., "Mean value of FT0A signal"};
  Configurable<float> meanFT0C{"meanFT0C", -1., "Mean value of FT0C signal"};

  //------------------------------------------------------------

  TRandom3* rand = new TRandom3(0);
  // Declare filter on collision Z vertex
  Filter collisionFilter = nabs(aod::jcollision::posZ) < vertexZCut;
  Filter collisionFilterMC = nabs(aod::jmccollision::posZ) < vertexZCut;

  // Declare filters on accepted tracks and MC particles (settings for jet reco are provided in the jet finder wagon)
  Filter trackFilter = aod::jtrack::pt > trkPtMin&& aod::jtrack::pt < trkPtMax&& nabs(aod::jtrack::eta) < trkEtaCut;
  Filter partFilter = nabs(aod::jmcparticle::eta) < trkEtaCut;

  // Declare filter on jets
  Filter jetRadiusFilter = aod::jet::r == nround(jetR.node() * 100.);

  HistogramRegistry spectra;

  std::vector<int> eventSelectionBits;
  int trackSelection = -1;

  Service<o2::framework::O2DatabasePDG> pdg;
  Preslice<FilteredMatchedJetsPartLevel> partJetsPerCollision = aod::jet::mcCollisionId;

  void init(InitContext const&)
  {
    // Initialize histogram axes
    AxisSpec pT{histJetPt + 20, -20.0, histJetPt * 1., "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec phiAngle{40, 0.0, constants::math::TwoPI, "#it{#varphi} (rad)"};
    AxisSpec etaTracks{100, -1.0, 1.0, "#it{#eta}_{trk}"};
    AxisSpec etaJets{100, -0.6, 0.6, "#it{#eta}_{jet}"};

    AxisSpec jetArea{50, 0.0, 3., "Area_{jet}"};
    AxisSpec rho{50, 0.0, 50., "#it{#rho}"};
    AxisSpec rhoArea{60, 0.0, 60., "#it{#rho} #times Area_{jet}"};

    // Convert configurable strings to std::string
    std::string evSelToString = static_cast<std::string>(evSel);
    std::string trkSelToString = static_cast<std::string>(trkSel);

    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(evSelToString);
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(trkSelToString);

    /*
    How to add histograms:

    spectra.add("histogram_name", "Histogram title", histogram_type (kTH1F/kTH1D and so on), {axis1, axis2, ...});
    where axis is definened as:

    AxisSpec axis_name{nbins, min_value, max_value, "axis_title"};
    or
    {nbins, min_value, max_value, "axis_title"}
    */

    if (doprocessFilteredCollisions || doprocessCollisions) // FK AKTUALIZOVAT PODLE PROCESS FUNKCI
    {
      // Event selection
      spectra.add("hEventSelectionCount", "Count # of events in the analysis", kTH1F, {{3, 0.0, 3.}});
      spectra.get<TH1>(HIST("hEventSelectionCount"))->GetXaxis()->SetBinLabel(1, "Total # of events w/o cuts");
      spectra.get<TH1>(HIST("hEventSelectionCount"))->GetXaxis()->SetBinLabel(2, Form("# of events after sel. %s", evSelToString.data()));
      spectra.get<TH1>(HIST("hEventSelectionCount"))->GetXaxis()->SetBinLabel(3, "# of events after z cut"); //?

      // Z coordinate of collision vertex
      spectra.add("hVertexZ_NoCut", "z vertex of collisions w/o cut", kTH1F, {{100, -20., 20., "#it{z}_{vertex}"}});
      spectra.add("hVertexZ_Cut", "z vertex of collisions w. cut", kTH1F, {{100, -20., 20., "#it{z}_{vertex}"}});
      spectra.add("hVertexZ_EventFiltering", "z vertex of collisions w. event filtering", kTH1F, {{100, -20., 20., "#it{z}_{vertex}"}});
    }

    // Disitribution of tracks
    if (doprocessTracks || doprocessJets) // FK AKTUALIZOVAT PODLE PROCESS FUNKCI I S OHLEDEM NA PARTICLE LEVEL
    {
      spectra.add("hPocetTracku", "Pocet", kTH1F, {{1000, 0., 1000., "Pocet Eventu"}});
      spectra.add("hPocetSignalTriggeru", "Pocet", kTH1F, {{10, 0., 10., "Pocet Eventu"}});
      spectra.add("hPocetReferenceTriggeru", "Pocet", kTH1F, {{10, 0., 10., "Pocet Eventu"}});
      spectra.add("hRecoilJetRefPt", "#it{p}_{T} distribution of recoil jets", kTH1F, {{200, 0., 200., "Pocet Eventu"}});
      spectra.add("hRecoilJetSigPt", "#it{p}_{T} distribution of recoil jets", kTH1F, {{200, 0., 200., "Pocet Eventu"}});
      spectra.add("hRecoilJetRefCorrPt", "Corrected #it{p}_{T} distribution of recoil jets", kTH1F, {{200, 0., 200., "Pocet Eventu"}});
      spectra.add("hRecoilJetSigCorrPt", "Corrected #it{p}_{T} distribution of recoil jets", kTH1F, {{200, 0., 200., "Pocet Eventu"}});
      spectra.add("hTrackPt", "#it{p}_{T} distribution of tracks", kTH1F, {pT});
      spectra.add("hTrackPhi", "#varphi distribution of tracks", kTH1F, {phiAngle});
      spectra.add("hTrackEta", "#eta distribution of tracks", kTH1F, {etaTracks});

      spectra.add("hTrackPtPhi", "#it{p}_{T} vs. #varphi distribution of tracks", kTH2F, {pT, phiAngle});
      spectra.add("hTrackPtEta", "#it{p}_{T} vs. #eta distribution of tracks", kTH2F, {pT, etaTracks});

      spectra.add("hTrackPtPhiEta", "#it{p}_{T} vs. #varphi vs. #eta distribution of tracks", kTH3F, {pT, phiAngle, etaTracks});

      // Distribution of jets
      spectra.add("hJetPt", "#it{p}_{T} distribution of jets", kTH1F, {pT});
      spectra.add("hJetPhi", "#varphi distribution of jets", kTH1F, {phiAngle});
      spectra.add("hJetEta", "#eta distribution of jets", kTH1F, {etaJets});

      spectra.add("hJetPtPhi", "#it{p}_{T} vs. #varphi distribution of jets", kTH2F, {pT, phiAngle});
      spectra.add("hJetPtEta", "#it{p}_{T} vs. #eta distribution of jets", kTH2F, {pT, etaJets});

      spectra.add("hJetPtPhiEta", "#it{p}_{T} vs. #varphi vs. #eta distribution of jets", kTH3F, {pT, phiAngle, etaJets});

      spectra.add("hMultFT0A", "Mult. signal from FTOA", kTH1F, {{2000, 0.0, 40000., "FT0A"}});
      spectra.add("hMultFT0C", "Mult. signal from FTOC", kTH1F, {{2000, 0.0, 40000., "FT0C"}});

      spectra.add("hJetPtCorr", "#it{p}_{T} distribution of jets", kTH1F, {pT});

      spectra.add("hRho", "distribution of rho", kTH1F, {rho});
      spectra.add("hjetArea", "distribution of jet area", kTH1F, {jetArea});

      spectra.add("JetAreavsPt", "#it{p}_{T} vs. jet area distribution of jets", kTH2F, {pT, jetArea});
      spectra.add("hRhovsNTracks", "Rho vs number of tracks", kTH2F, {{2000, 0.0, 40000., "NTracks"}, rho});

      spectra.add("hSigCorrConePt", "distribution of Signal Corrected perpendicular Cone Pt", kTH1F, {pT});
      spectra.add("hRefCorrConePt", "distribution of Reference Corrected perpendicular Cone Pt", kTH1F, {pT});

      spectra.add("hSigRandCorrConePt", "distribution of Signal Corrected random Cone Pt", kTH1F, {pT});
      spectra.add("hRefRandCorrConePt", "distribution of Reference Corrected random Cone Pt", kTH1F, {pT});

      spectra.add("hEventSelectionCountPartLevel", "Pocet", kTH1F, {{10, 0., 10., "Pocet Eventu"}});
      spectra.add("hMCPTrackPt", "#it{p}_{T} distribution of tracks", kTH1F, {pT});

      spectra.add("hJetPtEtaPhiRhoArea_Part", "Charact. of inclusive part. level jets", kTHnSparseF, {pT, etaJets, phiAngle, rho, jetArea});

      spectra.add("hRecoilJetRefPtMCP", "#it{p}_{T} distribution of recoil jets", kTH1F, {{200, 0., 200., "Pocet Eventu"}});
      spectra.add("hRecoilJetSigPtMCP", "#it{p}_{T} distribution of recoil jets", kTH1F, {{200, 0., 200., "Pocet Eventu"}});
      spectra.add("hRecoilJetRefCorrPtMCP", "Corrected #it{p}_{T} distribution of recoil jets", kTH1F, {{200, 0., 200., "Pocet Eventu"}});
      spectra.add("hRecoilJetSigCorrPtMCP", "Corrected #it{p}_{T} distribution of recoil jets", kTH1F, {{200, 0., 200., "Pocet Eventu"}});

      spectra.add("hTTCountPartLevel", "Pocet", kTH1F, {{2, 0., 2., "Pocet TT MCP"}});
      spectra.add("hTTCountData", "Pocet", kTH1F, {{2, 0., 2., "Pocet TT Data"}});

      spectra.add("hMissedJets_pT", "Part. level jets w/o matched pair", kTH1F, {{200, 0.0, 200.}}, setSumw2);
      spectra.add("hMissedJets_pT_RecoilJets", "Part. level jets w/o matched pair", kTH1F, {{200, 0.0, 200.}}, setSumw2);

      spectra.add("hJetPt_resolution", "Jet p_{T} relative resolution as a func. of jet #it{p}_{T, part}", kTH2F, {{100, -5., 5.}, pT}, setSumw2);
      spectra.add("hJetPhi_resolution", "#varphi resolution as a func. of jet #it{p}_{T, part}", kTH2F, {{40, -1., 1.}, pT}, setSumw2);

      spectra.add("hJetPt_DetLevel_vs_PartLevel_RecoilJets", "Correlation recoil jet pT at part. vs. det. levels", kTH2F, {{200, 0.0, 200.}, {200, 0.0, 200.}}, setSumw2);
      spectra.add("hJetPt_resolution_RecoilJets", "Jet p_{T} relative resolution as a func. of jet #it{p}_{T, part}", kTH2F, {{100, -5., 5.}, pT}, setSumw2);
      spectra.add("hJetPhi_resolution_RecoilJets", "#varphi resolution as a func. of jet #it{p}_{T, part}", kTH2F, {{40, -1., 1.}, pT}, setSumw2);
    }
  }

  //------------------------------------------------------------------------------
  /*
  In jet group, we do not use general information about collisions and tracks. We simply do not
  need all that information for jet finding and jet analysis. Therefore, we use ASoA tables with
  reduced information (aod::JetCollisions and aod::JetTracks). These tables are filled here: https://github.com/AliceO2Group/O2Physics/blob/4e284820d3b2f94327f0b71088c127decd74e8f9/PWGJE/TableProducer/derivedDataProducer.cxx

  Content of the tables can be found here: https://github.com/AliceO2Group/O2Physics/blob/cbc6e1f2eda0de5c8b00f64818e76eaa7e914232/PWGJE/DataModel/JetReducedData.h#L67



  */

  template <typename Collision, typename Jets, typename Tracks>
  void fillHistograms(Collision const& collision, Jets const& jets,
                      Tracks const& tracks)

  {

    if (skipEvent(collision)) {
      return;
    }

    std::vector<float> SigTrackPhi;
    std::vector<float> SigTrackEta;
    std::vector<float> RefTrackPhi;
    std::vector<float> RefTrackEta;

    float multFT0A = collision.multFT0A();
    float multFT0C = collision.multFT0C();
    // float multZNA = collision.multZNA();
    // float multZNC = collision.multZNC();

    spectra.fill(HIST("hMultFT0A"), multFT0A);
    spectra.fill(HIST("hMultFT0C"), multFT0C);

    float jetR2 = jetR * jetR;
    float rho = collision.rho();

    spectra.fill(HIST("hRho"), rho);

    int i = 0;
    for (const auto& track : tracks) {
      // check whether track passes the selection flags
      if (skipTrack(track))
        continue;
      i++;
      auto trackPt = track.pt();
      auto trackPhi = track.phi();
      auto trackEta = track.eta();

      spectra.fill(HIST("hTrackPt"), trackPt);
      spectra.fill(HIST("hTrackPhi"), trackPhi);
      spectra.fill(HIST("hTrackEta"), trackEta);

      spectra.fill(HIST("hTrackPtPhi"), trackPt, trackPhi);
      spectra.fill(HIST("hTrackPtEta"), trackPt, trackEta);

      spectra.fill(HIST("hTrackPtPhiEta"), trackPt, trackPhi, trackEta);

      if ((SignalTriggerMin <= trackPt) && (trackPt <= SignalTriggerMax)) {
        SigTrackPhi.push_back(trackPhi);
        SigTrackEta.push_back(trackEta);
      }

      if ((ReferenceTriggerMin <= trackPt) && (trackPt <= ReferenceTriggerMax)) {
        RefTrackPhi.push_back(trackPhi);
        RefTrackEta.push_back(trackEta);
      }
    }

    spectra.fill(HIST("hRhovsNTracks"), i, rho);

    float PerpConePhi = -99;
    float PhiTT = -999; // FK PSAT DO KODU KOMENTARE
    float EtaTT = -999;
    float rnd = rand->Rndm();
    bool analyzeSignal = 0;

    if (rnd < SigToRefFraction) {
      analyzeSignal = 1;
    }

    if (analyzeSignal == 0 && RefTrackPhi.size() > 0) {

      int ii = rand->Integer(RefTrackPhi.size());
      spectra.fill(HIST("hTTCountData"), 0.5);
      spectra.fill(HIST("hPocetReferenceTriggeru"), (float)RefTrackPhi.size());
      PhiTT = RefTrackPhi[ii];
      EtaTT = RefTrackEta[ii];

    }

    else if (analyzeSignal == 1 && SigTrackPhi.size() > 0) {

      int ii = rand->Integer(SigTrackPhi.size());
      spectra.fill(HIST("hTTCountData"), 1.5);
      spectra.fill(HIST("hPocetSignalTriggeru"), (float)SigTrackPhi.size());
      PhiTT = SigTrackPhi[ii];
      EtaTT = SigTrackEta[ii];
    }

    spectra.fill(HIST("hPocetTracku"), i);

    if (PhiTT < -100) {
      return;
    }

    PerpConePhi = TVector2::Phi_mpi_pi(PhiTT - 0.5 * TMath::Pi());

    float ConePt = 0.;
    for (const auto& track : tracks) {
      // check whether track passes the selection flags
      if (skipTrack(track))
        continue;

      auto trackPt = track.pt();
      auto trackPhi = track.phi();
      auto trackEta = track.eta();

      float DPhiPerp = TVector2::Phi_mpi_pi(trackPhi - PerpConePhi);
      float DEtaPerp = trackEta - EtaTT;

      if ((DPhiPerp) * (DPhiPerp) + (DEtaPerp) * (DEtaPerp) < jetR2) {
        ConePt = ConePt + trackPt; // FK KVULI E SCHEM BY SE MELY SCITAT 4vektory
      }
    }

    float CorrConePt = ConePt - TMath::Pi() * jetR2 * rho;
    if (analyzeSignal == 0) {
      spectra.fill(HIST("hRefCorrConePt"), CorrConePt);
    } else {
      spectra.fill(HIST("hSigCorrConePt"), CorrConePt);
    }

    int jj = 0;
    float Jet1Phi = 0.;
    float Jet1Eta = 0.;
    float Jet2Phi = 0.;
    float Jet2Eta = 0.;

    for (const auto& jet : jets) {

      jj++;

      auto jetPt = jet.pt();

      auto jetPhi = jet.phi();
      auto jetEta = jet.eta();
      float jetArea = jet.area();
      float jetPtCorr = jetPt - rho * jetArea;

      if (jj == 1) {
        Jet1Phi = jetPhi;
        Jet1Eta = jetEta;
      }

      if (jj == 2) {
        Jet2Phi = jetPhi;
        Jet2Eta = jetEta;
      }

      // FK TATO SPEKTRA SE PLNI JEN KDYZ JE TRIGGER. LEPE BY BYLO PLNIT SPEKTRA INKLUSIVNE

      spectra.fill(HIST("hJetPtCorr"), jetPtCorr);
      spectra.fill(HIST("hJetPt"), jetPt);
      spectra.fill(HIST("hJetPhi"), jetPhi);
      spectra.fill(HIST("hJetEta"), jetEta);

      spectra.fill(HIST("hJetPtPhi"), jetPt, jetPhi);
      spectra.fill(HIST("hJetPtEta"), jetPt, jetEta);

      spectra.fill(HIST("hJetPtPhiEta"), jetPt, jetPhi, jetEta);

      spectra.fill(HIST("hjetArea"), jetArea);
      spectra.fill(HIST("JetAreavsPt"), jetPt, jetArea);

      float DPhi = -100;
      DPhi = TMath::Abs(TVector2::Phi_mpi_pi(jetPhi - PhiTT));

      if (DPhi > TMath::Pi() - DPhiCut) {
        if (analyzeSignal == 0) {
          spectra.fill(HIST("hRecoilJetRefPt"), jetPt);
          spectra.fill(HIST("hRecoilJetRefCorrPt"), jetPtCorr);
        } else {
          spectra.fill(HIST("hRecoilJetSigPt"), jetPt);
          spectra.fill(HIST("hRecoilJetSigCorrPt"), jetPtCorr);
        }
      }
    }

    // printf("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
    bool bCloseJet = true;
    float randJetPhi = -999.;
    float randJetEta = -999.;
    while (bCloseJet) { // FK dat nejake pocitadlo, ktere umozni opustit loop
      randJetPhi = rand->Uniform(-TMath::Pi(), TMath::Pi());
      randJetEta = rand->Uniform(-0.5, 0.5);

      float vphi1 = TVector2::Phi_mpi_pi(randJetPhi - Jet1Phi);
      float vphi2 = TVector2::Phi_mpi_pi(randJetPhi - Jet2Phi);
      float veta1 = randJetEta - Jet1Eta;
      float veta2 = randJetEta - Jet2Eta;
      float dist1 = sqrt(vphi1 * vphi1 + veta1 * veta1);
      float dist2 = sqrt(vphi2 * vphi2 + veta2 * veta2);

      if (dist1 > 2 * jetR && dist2 > 2 * jetR)
        bCloseJet = false;
    }

    // printf("yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy\n");

    float RandConePt = 0.;
    for (const auto& track : tracks) {
      if (skipTrack(track))
        continue;

      auto trackPt = track.pt();
      auto trackPhi = track.phi();
      auto trackEta = track.eta();
      float DPhiRand = TVector2::Phi_mpi_pi(trackPhi - randJetPhi);
      float DEtaRand = trackEta - randJetEta;

      if ((DPhiRand) * (DPhiRand) + (DEtaRand) * (DEtaRand) < jetR2) {
        RandConePt = RandConePt + trackPt; // FK scitat 4momenty
      }
    }

    float CorrRandConePt = RandConePt - TMath::Pi() * jetR2 * rho;

    if (analyzeSignal == 0) {
      spectra.fill(HIST("hRefRandCorrConePt"), CorrRandConePt);
    } else {
      spectra.fill(HIST("hSigRandCorrConePt"), CorrRandConePt);
    }
  }

  template <typename Collision, typename Jets, typename Particles>
  void fillMatchedHistograms(JTracksTable const& tracks,
                             JetsBase const& jetsBase,
                             JetsTag const& jetsTag,
                             float weight = 1.)
  {
    std::vector<double> SigTrackPhi;
    double PhiTT = 0.0;

    for (const auto& track : tracks) {
      if (skipTrack(track))
        continue;

      if ((SignalTriggerMin <= trackPt) && (trackPt <= SignalTriggerMax))
        SigTrackPhi.push_back(track.phi());
    }
    bool bIsThereTTSig = SigTrackPhi.size() > 0;

    if (bIsThereTTSig)
      PhiTT = [rand->Integer(SigTrackPhi.size())];

    for (const auto& jetBase : jetsBase) {
      bool bIsBaseJetRecoil = false;

      if (bIsThereTTSig) {
        float DPhi = TMath::Abs(TVector2::Phi_mpi_pi(jetBase.phi() - PhiTT));
        if (DPhi > TMath::Pi() - DPhiCut)
          bIsBaseJetRecoil = true;
      }

      dataForUnfolding(jetBase, jetsTag, bIsBaseJetRecoil, tracks, weight);
    }
  }

  void processCollisions(aod::JetCollisions::iterator const& collision)
  { // no filtering on collisions

    spectra.fill(HIST("hEventSelectionCount"), 0.5); // all events without cuts

    // check whether event contains the required flag for event selection (sel8 in this case)
    if (skipEvent(collision))
      return;

    spectra.fill(HIST("hEventSelectionCount"), 1.5); // number of events after applying event selection
    spectra.fill(HIST("hVertexZ_NoCut"), collision.posZ());

    // Imperative way to apply z-vertex cut (ineffective way, much faster is to apply the filter)
    if (fabs(collision.posZ()) > vertexZCut) {
      return;
    }

    spectra.fill(HIST("hEventSelectionCount"), 2.5); // number of events after applying event selection + z-vertex cut
    spectra.fill(HIST("hVertexZ_Cut"), collision.posZ());
  }
  PROCESS_SWITCH(TrackJetSpectra, processCollisions, "process collisions", true);

  void processFilteredCollisions(soa::Filtered<aod::JetCollisions>::iterator const& collision)
  {
    // check whether event contains the required flag for event selection (sel8 in this case)
    if (skipEvent(collision))
      return;

    spectra.fill(HIST("hVertexZ_EventFiltering"), collision.posZ()); // compare with "hVertexZ_Cut" histogram
  }
  PROCESS_SWITCH(TrackJetSpectra, processFilteredCollisions, "process filtered collisions", false);

  void processTracks(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                     soa::Filtered<aod::JetTracks> const& tracks)
  {
    // check whether event contains the required flag for event selection (sel8 in this case)
    if (skipEvent(collision))
      return;

    int i = 0;
    int PocetSigTrigeru = 0;
    int PocetRefTrigeru = 0;
    for (const auto& track : tracks) {
      // check whether track passes the selection flags
      if (skipTrack(track))
        continue;
      i++;
      auto trackPt = track.pt();
      auto trackPhi = track.phi();
      auto trackEta = track.eta();

      spectra.fill(HIST("hTrackPt"), trackPt);
      spectra.fill(HIST("hTrackPhi"), trackPhi);
      spectra.fill(HIST("hTrackEta"), trackEta);

      spectra.fill(HIST("hTrackPtPhi"), trackPt, trackPhi);
      spectra.fill(HIST("hTrackPtEta"), trackPt, trackEta);

      spectra.fill(HIST("hTrackPtPhiEta"), trackPt, trackPhi, trackEta);

      if ((SignalTriggerMin <= trackPt) && (trackPt <= SignalTriggerMax))
        PocetSigTrigeru++;

      if ((ReferenceTriggerMin <= trackPt) && (trackPt <= ReferenceTriggerMax))
        PocetRefTrigeru++;
    }

    spectra.fill(HIST("hPocetReferenceTriggeru"), PocetRefTrigeru);
    spectra.fill(HIST("hPocetSignalTriggeru"), PocetSigTrigeru);
    spectra.fill(HIST("hPocetTracku"), i);
  }

  PROCESS_SWITCH(TrackJetSpectra, processTracks, "process tracks", false);

  void processJets(FilteredColl const& collision,
                   soa::Filtered<aod::ChargedJets> const& jets,
                   soa::Filtered<aod::JetTracks> const& tracks)
  {
    // check whether event contains the required flag for event selection (sel8 in this case)

    fillHistograms(collision, jets, tracks);
  }

  PROCESS_SWITCH(TrackJetSpectra, processJets, "process inclusive jets", false);

  void processJetsMCDet(FilteredColl const& collision,
                        FilteredJetsDetLevel const& jets,
                        soa::Filtered<aod::JetTracks> const& tracks)
  {
    // check whether event contains the required flag for event selection (sel8 in this case)

    fillHistograms(collision, jets, tracks);
  }

  void processMCPartLevel(FilteredCollPartLevel const& collision,
                          FilteredParticles const& particles,
                          FilteredJetsPartLevel const& jets)
  {
    spectra.fill(HIST("hEventSelectionCountPartLevel"), 0.5);

    std::vector<float> vPhiOfTT;

    float rnd = rand->Rndm();
    bool analyzeSignal = 0;
    float rho = collision.rho();

    if (rnd < SigToRefFraction) {
      analyzeSignal = 1;
    }

    for (const auto& particle : particles) {
      auto pdgParticle = pdg->GetParticle(particle.pdgCode());
      if (!pdgParticle)
        continue;

      // Need charge and physical primary particles
      bool bParticleNeutral = (static_cast<int8_t>(pdgParticle->Charge()) == 0);
      if (bParticleNeutral || !particle.isPhysicalPrimary())
        continue;

      float particlePt = particle.pt();
      float particlePhi = particle.phi();

      if (analyzeSignal && (particlePt > SignalTriggerMin && particlePt < SignalTriggerMax)) {
        vPhiOfTT.push_back(particlePhi);
      }

      if (!analyzeSignal && (particlePt > ReferenceTriggerMin && particlePt < ReferenceTriggerMax)) {
        vPhiOfTT.push_back(particlePhi);
      }
    }
    int ii = -1;
    if (vPhiOfTT.size() > 0) {
      ii = rand->Integer(vPhiOfTT.size());
      if (analyzeSignal == 0) {
        spectra.fill(HIST("hTTCountPartLevel"), 0.5);
      } else {
        spectra.fill(HIST("hTTCountPartLevel"), 1.5);
      }
    }

    for (const auto& jet : jets) {
      float jetPt = jet.pt();
      float jetPhi = jet.phi();
      float jetArea = jet.area();
      float jetPtCorr = jetPt - rho * jetArea;

      spectra.fill(HIST("hJetPtEtaPhiRhoArea_Part"), jetPt, jet.eta(), jet.phi(), rho, jetArea);

      if (ii > -1) {

        float PhiTT = vPhiOfTT[ii];

        float DeltaPhi = TMath::Abs(TVector2::Phi_mpi_pi(PhiTT - jetPhi));
        if (DeltaPhi > TMath::Pi() - DPhiCut) {
          if (analyzeSignal == 0) {
            spectra.fill(HIST("hRecoilJetRefPtMCP"), jetPt);
            spectra.fill(HIST("hRecoilJetRefCorrPtMCP"), jetPtCorr);
          } else {
            spectra.fill(HIST("hRecoilJetSigPtMCP"), jetPt);
            spectra.fill(HIST("hRecoilJetSigCorrPtMCP"), jetPtCorr);
          }
        }
      }
    }
  }

  PROCESS_SWITCH(TrackJetSpectra, processJetsMCDet, "process inclusive jets", false);

  void processJetsMatched(FilteredCollDetLevelGetWeight const& collision,
                          aod::JetMcCollisions const&,
                          FilteredTracks const& tracks,
                          FilteredMatchedJetsDetLevel const& mcdjets,
                          FilteredMatchedJetsPartLevel const& mcpjets)
  {
    if (skipEvent(collision))
      return;

    auto mcpjetsPerMCCollision = mcpjets.sliceBy(partJetsPerCollision, collision.mcCollisionId());

    fillMatchedHistograms(tracks, mcpjetsPerMCCollision, mcdjets);
  }
  PROCESS_SWITCH(TrackJetSpectra, processJetsMatched, "process matching of MC jets (no weight)", false);

  //------------------------------------------------------------------------------
  // Auxiliary functions
  template <typename Collision>
  bool skipEvent(const Collision& coll)
  {
    /// \brief: trigger cut is needed for pp data
    return !jetderiveddatautilities::selectCollision(coll, eventSelectionBits);
  }

  template <typename Track>
  bool skipTrack(const Track& track)
  {
    return !jetderiveddatautilities::selectTrack(track, trackSelection);
  }

  template <typename Jet, typename Tracks>
  bool isJetWithHighPtConstituent(Jet const& jet, Tracks const&)
  {
    bool bIsJetWithHighPtConstituent = false;
    for (const auto& jetConstituent : jet.template tracks_as<Tracks>()) {
      if (jetConstituent.pt() > maxJetConstituentPt) {
        bIsJetWithHighPtConstituent = true;
        break;
      }
    }
    return bIsJetWithHighPtConstituent;
  }

  template <typename PartJet, typename DetJet, typename TracksTable>
  void dataForUnfolding(PartJet const& partJet, DetJet const& detJets, bool bIsBaseJetRecoil, TracksTable const& tracks, float weight = 1.)
  {

    float partJetPt = partJet.pt();
    bool bIsThereMatchedJet = partJet.has_matchedJetGeo();

    // FK mame histogram particle level  pT

    if (bIsThereMatchedJet) {
      const auto& jetsMatched = partJet.template matchedJetGeo_as<std::decay_t<DetJet>>();

      for (const auto& jetMatched : jetsMatched) {
        // skip matches where detector level jets have a constituent with pT above specified cut
        bool skipMatchedDetJet = isJetWithHighPtConstituent(jetMatched, tracks);

        if (skipMatchedDetJet) {
          // Miss jets
          spectra.fill(HIST("hMissedJets_pT"), partJetPt, weight);
          if (bIsBaseJetRecoil)
            spectra.fill(HIST("hMissedJets_pT_RecoilJets"), partJetPt, weight);
        } else {
          float detJetPt = jetMatched.pt();

          spectra.fill(HIST("hJetPt_resolution"), (partJetPt - detJetPt) / partJetPt, partJetPt, weight);
          spectra.fill(HIST("hJetPhi_resolution"), partJet.phi() - jetMatched.phi(), partJetPt, weight);

          if (bIsBaseJetRecoil) {
            spectra.fill(HIST("hJetPt_DetLevel_vs_PartLevel_RecoilJets"), detJetPt, partJetPt, weight);
            spectra.fill(HIST("hJetPt_resolution_RecoilJets"), (partJetPt - detJetPt) / partJetPt, partJetPt, weight);
            spectra.fill(HIST("hJetPhi_resolution_RecoilJets"), partJet.phi() - jetMatched.phi(), partJetPt, weight);
          }
        }
      }
    }
  }

  // template <typename Jet, typename Tracks>
  // bool isJetWithHighPtConstituent(Jet const& jet, Tracks const&)
  // {
  //   bool bIsJetWithHighPtConstituent = false;
  //   for (const auto& jetConstituent : jet.template tracks_as<Tracks>()) {
  //     if (jetConstituent.pt() > maxJetConstituentPt) {
  //       bIsJetWithHighPtConstituent = true;
  //       break;
  //     }
  //   }
  //   return bIsJetWithHighPtConstituent;
  // }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<TrackJetSpectra>(cfgc)};
}
