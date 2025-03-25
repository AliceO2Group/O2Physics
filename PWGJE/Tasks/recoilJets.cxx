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

/// \author Kotliarov Artem <artem.kotliarov@cern.ch>
/// \file recoilJets.cxx
/// \brief hadron-jet correlation analysis

#include <string>
#include <tuple>
#include <vector>

#include "TRandom3.h"
#include "TVector2.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "CommonConstants/MathConstants.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/DataModel/Jet.h"

#include "PWGJE/Core/JetDerivedDataUtilities.h"

#include "EventFiltering/filterTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// Shorthand notations
using FilteredColl = soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos>>::iterator;
using FilteredCollPartLevel = soa::Filtered<soa::Join<aod::JetMcCollisions, aod::BkgChargedMcRhos>>::iterator;
using FilteredCollDetLevelGetWeight = soa::Filtered<soa::Join<aod::JetCollisionsMCD, aod::BkgChargedRhos>>::iterator;

using FilteredJets = soa::Filtered<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>>;
using FilteredJetsDetLevel = soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>>;
using FilteredJetsPartLevel = soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents>>;

using FilteredMatchedJetsDetLevel = soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>>;
using FilteredMatchedJetsPartLevel = soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>>;

using FilteredTracks = soa::Filtered<aod::JetTracks>;

struct RecoilJets {

  // List of configurable parameters
  Configurable<std::string> evSel{"evSel", "sel8", "Choose event selection"};
  Configurable<std::string> trkSel{"trkSel", "globalTracks", "Set track selection"};
  Configurable<float> vertexZCut{"vertexZCut", 10., "Accepted z-vertex range"};
  Configurable<float> fracSig{"fracSig", 0.9, "Fraction of events to use for signal TT"};

  Configurable<float> trkPtMin{"trkPtMin", 0.15, "Minimum pT of acceptanced tracks"};
  Configurable<float> trkPtMax{"trkPtMax", 100., "Maximum pT of acceptanced tracks"};

  Configurable<float> trkEtaCut{"trkEtaCut", 0.9, "Eta acceptance of TPC"};
  Configurable<float> jetR{"jetR", 0.4, "Jet cone radius"};

  Configurable<std::string> triggerMasks{"triggerMasks", "", "Relevant trigger masks: fTrackLowPt,fTrackHighPt"};

  // List of configurable parameters for MC
  Configurable<float> pTHatExponent{"pTHatExponent", 4.0, "Exponent of the event weight for the calculation of pTHat"};
  Configurable<float> pTHatMax{"pTHatMax", 999.0, "Maximum fraction of hard scattering for jet acceptance in MC"};

  // Parameters for recoil jet selection
  Configurable<float> ptTTrefMin{"ptTTrefMin", 5., "Minimum pT of reference TT"};
  Configurable<float> ptTTrefMax{"ptTTrefMax", 7., "Maximum pT of reference TT"};
  Configurable<float> ptTTsigMin{"ptTTsigMin", 20., "Minimum pT of signal TT"};
  Configurable<float> ptTTsigMax{"ptTTsigMax", 50., "Maximum pT of signal TT"};
  Configurable<float> recoilRegion{"recoilRegion", 0.6, "Width of recoil acceptance"};

  // List of configurable parameters for histograms
  Configurable<uint16_t> histJetPt{"histJetPt", 100, "Maximum value of jet pT shown in histograms"};

  // Axes specification
  AxisSpec pT{histJetPt, 0.0, histJetPt * 1.0, "#it{p}_{T} (GeV/#it{c})"};
  AxisSpec jetPTcorr{histJetPt + 20, -20., histJetPt * 1.0, "#it{p}_{T, jet}^{ch, corr} (GeV/#it{c})"};
  AxisSpec phiAngle{40, 0.0, constants::math::TwoPI, "#it{#varphi} (rad)"};
  AxisSpec deltaPhiAngle{52, 0.0, constants::math::PI, "#Delta#it{#varphi} (rad)"};
  AxisSpec pseudorap{40, -1., 1., "#it{#eta}"};
  AxisSpec pseudorapJets{20, -0.5, 0.5, "#it{#eta}_{jet}"};
  AxisSpec jetArea{50, 0.0, 5., "Area_{jet}"};
  AxisSpec rhoArea{60, 0.0, 60., "#it{#rho} #times Area_{jet}"};
  AxisSpec rho{50, 0.0, 50., "#it{#rho}"};

  Preslice<FilteredMatchedJetsPartLevel> partJetsPerCollision = aod::jet::mcCollisionId;

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
  std::vector<int> triggerMaskBits;

  Service<o2::framework::O2DatabasePDG> pdg;

  void init(InitContext const&)
  {
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(evSel));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trkSel));
    triggerMaskBits = jetderiveddatautilities::initialiseTriggerMaskBits(triggerMasks);

    // List of raw and MC det. distributions
    if (doprocessData || doprocessMCDetLevel || doprocessMCDetLevelWeighted) {
      spectra.add("vertexZ", "Z vertex of collisions", kTH1F, {{60, -12., 12.}});

      spectra.add("hTrackPtEtaPhi", "Charact. of tracks", kTH3F, {pT, pseudorap, phiAngle});
      spectra.add("hNtrig", "Total number of selected triggers per class", kTH1F, {{2, 0.0, 2.}});
      spectra.get<TH1>(HIST("hNtrig"))->GetXaxis()->SetBinLabel(1, "TT_{ref}");
      spectra.get<TH1>(HIST("hNtrig"))->GetXaxis()->SetBinLabel(2, "TT_{sig}");

      spectra.add("hTTRef_per_event", "Number of TT_{Ref} per event", kTH1F, {{15, 0.5, 15.5}});
      spectra.add("hTTSig_per_event", "Number of TT_{Sig} per event", kTH1F, {{10, 0.5, 10.5}});

      spectra.add("hJetPtEtaPhiRhoArea", "Charact. of inclusive jets", kTHnSparseF, {pT, pseudorapJets, phiAngle, rho, jetArea});

      spectra.add("hDPhi_JetPt_Corr_TTRef", "Events w. TT_{Ref}: #Delta#varphi & #it{p}_{T, jet}^{ch}", kTH2F, {deltaPhiAngle, jetPTcorr});
      spectra.add("hDPhi_JetPt_Corr_TTSig", "Events w. TT_{Sig}: #Delta#varphi & #it{p}_{T, jet}^{ch}", kTH2F, {deltaPhiAngle, jetPTcorr});
      spectra.add("hDPhi_JetPt_TTRef", "Events w. TT_{Ref}: #Delta#varphi & #it{p}_{T, jet}^{ch}", kTH2F, {deltaPhiAngle, pT});
      spectra.add("hDPhi_JetPt_TTSig", "Events w. TT_{Sig}: #Delta#varphi & #it{p}_{T, jet}^{ch}", kTH2F, {deltaPhiAngle, pT});

      spectra.add("hRecoil_JetPt_Corr_TTRef", "Events w. TT_{Ref}: #it{p}_{T} of recoil jets", kTH1F, {jetPTcorr});
      spectra.add("hRecoil_JetPt_Corr_TTSig", "Events w. TT_{Sig}: #it{p}_{T} of recoil jets", kTH1F, {jetPTcorr});
      spectra.add("hRecoil_JetPt_TTRef", "Events w. TT_{Ref}: #it{p}_{T} of recoil jets", kTH1F, {pT});
      spectra.add("hRecoil_JetPt_TTSig", "Events w. TT_{Sig}: #it{p}_{T} of recoil jets", kTH1F, {pT});

      spectra.add("hJetArea_JetPt_Rho_TTRef", "Events w. TT_{Ref}: A_{jet} & jet pT & #rho", kTH3F, {jetArea, pT, rho});
      spectra.add("hJetArea_JetPt_Rho_TTSig", "Events w. TT_{Sig}: A_{jet} & jet pT & #rho", kTH3F, {jetArea, pT, rho});
    }

    // List of MC particle level distributions
    if (doprocessMCPartLevel || doprocessMCPartLevelWeighted) {
      spectra.add("hPartPtEtaPhi", "Charact. of particles", kTH3F, {pT, pseudorap, phiAngle});
      spectra.add("hNtrig_Part", "Total number of selected triggers per class", kTH1F, {{2, 0.0, 2.}});
      spectra.get<TH1>(HIST("hNtrig_Part"))->GetXaxis()->SetBinLabel(1, "TT_{ref}");
      spectra.get<TH1>(HIST("hNtrig_Part"))->GetXaxis()->SetBinLabel(2, "TT_{sig}");

      spectra.add("hTTRef_per_event_Part", "Number of TT_{Ref} per event", kTH1F, {{15, 0.5, 15.5}});
      spectra.add("hTTSig_per_event_Part", "Number of TT_{Sig} per event", kTH1F, {{10, 0.5, 10.5}});

      spectra.add("hJetPtEtaPhiRhoArea_Part", "Charact. of inclusive part. level jets", kTHnSparseF, {pT, pseudorapJets, phiAngle, rho, jetArea});

      spectra.add("hDPhi_JetPt_Corr_TTRef_Part", "Events w. TT_{Ref}: #Delta#varphi & #it{p}_{T, jet}^{ch}", kTH2F, {deltaPhiAngle, jetPTcorr});
      spectra.add("hDPhi_JetPt_Corr_TTSig_Part", "Events w. TT_{Sig}: #Delta#varphi & #it{p}_{T, jet}^{ch}", kTH2F, {deltaPhiAngle, jetPTcorr});
      spectra.add("hDPhi_JetPt_TTRef_Part", "Events w. TT_{Ref}: #Delta#varphi & #it{p}_{T, jet}^{ch}", kTH2F, {deltaPhiAngle, pT});
      spectra.add("hDPhi_JetPt_TTSig_Part", "Events w. TT_{Sig}: #Delta#varphi & #it{p}_{T, jet}^{ch}", kTH2F, {deltaPhiAngle, pT});

      spectra.add("hRecoil_JetPt_Corr_TTRef_Part", "Events w. TT_{Ref}: #it{p}_{T} of recoil jets", kTH1F, {jetPTcorr});
      spectra.add("hRecoil_JetPt_Corr_TTSig_Part", "Events w. TT_{Sig}: #it{p}_{T} of recoil jets", kTH1F, {jetPTcorr});
      spectra.add("hRecoil_JetPt_TTRef_Part", "Events w. TT_{Ref}: #it{p}_{T} of recoil jets", kTH1F, {pT});
      spectra.add("hRecoil_JetPt_TTSig_Part", "Events w. TT_{Sig}: #it{p}_{T} of recoil jets", kTH1F, {pT});

      spectra.add("hJetArea_JetPt_Rho_TTRef_Part", "Events w. TT_{Ref}: A_{jet} & jet pT & #rho", kTH3F, {jetArea, pT, rho});
      spectra.add("hJetArea_JetPt_Rho_TTSig_Part", "Events w. TT_{Sig}: A_{jet} & jet pT & #rho", kTH3F, {jetArea, pT, rho});
    }

    // Jet matching: part. vs. det.
    if (doprocessJetsMatched || doprocessJetsMatchedWeighted) {
      spectra.add("hJetPt_PartLevel_vs_DetLevel", "Correlation jet pT at part. vs. det. levels", kTH2F, {{200, 0.0, 200.}, {200, 0.0, 200.}});
      // spectra.add("hJetPt_Corr_PartLevel_vs_DetLevel", "Correlation jet pT at part. vs. det. levels", kTH2F, {jetPTcorr, jetPTcorr});
      spectra.add("hJetPt_PartLevel_vs_DetLevel_RecoilJets", "Correlation recoil jet pT at part. vs. det. levels", kTH2F, {{200, 0.0, 200.}, {200, 0.0, 200.}});
      // spectra.add("hJetPt_Corr_PartLevel_vs_DetLevel_RecoilJets", "Correlation recoil jet pT at part. vs. det. levels", kTH2F, {jetPTcorr, jetPTcorr});

      spectra.add("hMissedJets_pT", "Part. level jets w/o matched pair", kTH1F, {{200, 0.0, 200.}});
      // spectra.add("hMissedJets_Corr_pT", "Part. level jets w/o matched pair", kTH1F, {jetPTcorr});
      // spectra.add("hMissedJets_pT_RecoilJets", "Part. level jets w/o matched pair", kTH1F, {{200, 0.0, 200.}});
      // spectra.add("hMissedJets_Corr_pT_RecoilJets", "Part. level jets w/o matched pair", kTH1F, {jetPTcorr});

      spectra.add("hFakeJets_pT", "Det. level jets w/o matched pair", kTH1F, {{200, 0.0, 200.}});
      // spectra.add("hFakeJets_Corr_pT", "Det. level jets w/o matched pair", kTH1F, {jetPTcorr});
      spectra.add("hFakeJets_pT_RecoilJets", "Det. level jets w/o matched pair", kTH1F, {{200, 0.0, 200.}});
      // spectra.add("hFakeJets_Corr_pT_RecoilJets", "Det. level jets w/o matched pair", kTH1F, {jetPTcorr});

      spectra.add("hJetPt_resolution", "Jet p_{T} relative resolution as a func. of jet #it{p}_{T, part}", kTH2F, {{100, -5., 5.}, pT});
      spectra.add("hJetPt_resolution_RecoilJets", "Jet p_{T} relative resolution as a func. of jet #it{p}_{T, part}", kTH2F, {{100, -5., 5.}, pT});

      spectra.add("hJetPhi_resolution", "#varphi resolution as a func. of jet #it{p}_{T, part}", kTH2F, {{40, -1., 1.}, pT});
      spectra.add("hJetPhi_resolution_RecoilJets", "#varphi resolution as a func. of jet #it{p}_{T, part}", kTH2F, {{40, -1., 1.}, pT});
    }
  }

  // Fill histograms with raw or MC det. level data
  template <typename Collision, typename Jets, typename Tracks>
  void fillHistograms(Collision const& collision, Jets const& jets, Tracks const& tracks, bool bIsMC = false, float weight = 1.)
  {

    bool bSigEv = false;
    std::vector<double> vPhiOfTT;
    double phiTT = 0.;
    int nTT = 0;
    float pTHat = 0.;
    if (bIsMC)
      pTHat = getPtHat(weight);

    auto dice = rand->Rndm();
    if (dice < fracSig)
      bSigEv = true;

    for (const auto& track : tracks) {
      if (skipTrack(track))
        continue;

      spectra.fill(HIST("hTrackPtEtaPhi"), track.pt(), track.eta(), track.phi(), weight);

      // Search for TT candidate
      if (bSigEv && (track.pt() > ptTTsigMin && track.pt() < ptTTsigMax)) {
        vPhiOfTT.push_back(track.phi());
        ++nTT;
      }

      if (!bSigEv && (track.pt() > ptTTrefMin && track.pt() < ptTTrefMax)) {
        vPhiOfTT.push_back(track.phi());
        ++nTT;
      }
    }

    if (nTT > 0) { // at least 1 TT

      phiTT = getPhiTT(vPhiOfTT);

      if (bSigEv) {
        spectra.fill(HIST("hNtrig"), 1.5, weight);
        spectra.fill(HIST("hTTSig_per_event"), nTT, weight);
      } else {
        spectra.fill(HIST("hNtrig"), 0.5, weight);
        spectra.fill(HIST("hTTRef_per_event"), nTT, weight);
      }
    }

    for (const auto& jet : jets) {

      if (bIsMC && (jet.pt() > pTHatMax * pTHat))
        continue;

      spectra.fill(HIST("hJetPtEtaPhiRhoArea"), jet.pt(), jet.eta(), jet.phi(), collision.rho(), jet.area(), weight);

      if (nTT > 0) {
        auto [dphi, bRecoilJet] = isRecoilJet(jet, phiTT);

        if (bSigEv) {

          spectra.fill(HIST("hDPhi_JetPt_Corr_TTSig"), dphi, jet.pt() - collision.rho() * jet.area(), weight);
          spectra.fill(HIST("hDPhi_JetPt_TTSig"), dphi, jet.pt(), weight);
          spectra.fill(HIST("hJetArea_JetPt_Rho_TTSig"), jet.area(), jet.pt(), collision.rho(), weight);

          if (bRecoilJet) {
            spectra.fill(HIST("hRecoil_JetPt_Corr_TTSig"), jet.pt() - collision.rho() * jet.area(), weight);
            spectra.fill(HIST("hRecoil_JetPt_TTSig"), jet.pt(), weight);
          }

        } else {
          spectra.fill(HIST("hDPhi_JetPt_Corr_TTRef"), dphi, jet.pt() - collision.rho() * jet.area(), weight);
          spectra.fill(HIST("hDPhi_JetPt_TTRef"), dphi, jet.pt(), weight);
          spectra.fill(HIST("hJetArea_JetPt_Rho_TTRef"), jet.area(), jet.pt(), collision.rho(), weight);

          if (bRecoilJet) {
            spectra.fill(HIST("hRecoil_JetPt_Corr_TTRef"), jet.pt() - collision.rho() * jet.area(), weight);
            spectra.fill(HIST("hRecoil_JetPt_TTRef"), jet.pt(), weight);
          }
        }
      }
    }
  }

  template <typename Collision, typename Jets, typename Particles>
  void fillMCPHistograms(Collision const& collision, Jets const& jets, Particles const& particles, float weight = 1.)
  {
    bool bSigEv = false;
    std::vector<double> vPhiOfTT;
    double phiTT = 0.;
    int nTT = 0;
    float pTHat = getPtHat(weight);

    auto dice = rand->Rndm();
    if (dice < fracSig)
      bSigEv = true;

    for (const auto& particle : particles) {
      auto pdgParticle = pdg->GetParticle(particle.pdgCode());
      if (!pdgParticle)
        continue;

      // Need charge and physical primary particles
      bool bParticleNeutral = (static_cast<int8_t>(pdgParticle->Charge()) == 0);
      if (bParticleNeutral || !particle.isPhysicalPrimary())
        continue;

      spectra.fill(HIST("hPartPtEtaPhi"), particle.pt(), particle.eta(), particle.phi(), weight);

      if (bSigEv && (particle.pt() > ptTTsigMin && particle.pt() < ptTTsigMax)) {
        vPhiOfTT.push_back(particle.phi());
        ++nTT;
      }

      if (!bSigEv && (particle.pt() > ptTTrefMin && particle.pt() < ptTTrefMax)) {
        vPhiOfTT.push_back(particle.phi());
        ++nTT;
      }
    }

    if (nTT > 0) {

      phiTT = getPhiTT(vPhiOfTT);

      if (bSigEv) {
        spectra.fill(HIST("hNtrig_Part"), 1.5, weight);
        spectra.fill(HIST("hTTSig_per_event_Part"), nTT, weight);
      } else {
        spectra.fill(HIST("hNtrig_Part"), 0.5, weight);
        spectra.fill(HIST("hTTRef_per_event_Part"), nTT, weight);
      }
    }

    for (const auto& jet : jets) {

      if (jet.pt() > pTHatMax * pTHat)
        continue;

      spectra.fill(HIST("hJetPtEtaPhiRhoArea_Part"), jet.pt(), jet.eta(), jet.phi(), collision.rho(), jet.area(), weight);

      if (nTT > 0) {

        auto [dphi, bRecoilJet] = isRecoilJet(jet, phiTT);

        if (bSigEv) {

          spectra.fill(HIST("hDPhi_JetPt_Corr_TTSig_Part"), dphi, jet.pt() - collision.rho() * jet.area(), weight);
          spectra.fill(HIST("hDPhi_JetPt_TTSig_Part"), dphi, jet.pt(), weight);
          spectra.fill(HIST("hJetArea_JetPt_Rho_TTSig_Part"), jet.area(), jet.pt(), collision.rho(), weight);

          if (bRecoilJet) {
            spectra.fill(HIST("hRecoil_JetPt_Corr_TTSig_Part"), jet.pt() - collision.rho() * jet.area(), weight);
            spectra.fill(HIST("hRecoil_JetPt_TTSig_Part"), jet.pt(), weight);
          }

        } else {

          spectra.fill(HIST("hDPhi_JetPt_Corr_TTRef_Part"), dphi, jet.pt() - collision.rho() * jet.area(), weight);
          spectra.fill(HIST("hDPhi_JetPt_TTRef_Part"), dphi, jet.pt(), weight);
          spectra.fill(HIST("hJetArea_JetPt_Rho_TTRef_Part"), jet.area(), jet.pt(), collision.rho(), weight);

          if (bRecoilJet) {
            spectra.fill(HIST("hRecoil_JetPt_Corr_TTRef_Part"), jet.pt() - collision.rho() * jet.area(), weight);
            spectra.fill(HIST("hRecoil_JetPt_TTRef_Part"), jet.pt(), weight);
          }
        }
      }
    }
  }

  /// TODO: Add functionality to get rho for particle and detector level
  template <typename Tracks, typename DetLevelJets, typename PartLevelJets>
  void fillMatchedHistograms(Tracks const& tracks, DetLevelJets const& jets_det_level, PartLevelJets const& jets_part_level, float weight = 1.)
  {
    std::vector<double> vPhiOfTT;
    double phiTT = 0.;
    float pTHat = getPtHat(weight);

    for (const auto& track : tracks) {
      if (skipTrack(track))
        continue;

      if (track.pt() > ptTTsigMin && track.pt() < ptTTsigMax) {
        vPhiOfTT.push_back(track.phi());
      }
    }

    bool bTT = vPhiOfTT.size() > 0;
    if (bTT)
      phiTT = getPhiTT(vPhiOfTT);

    for (const auto& jet_det_level : jets_det_level) {
      if (jet_det_level.pt() > pTHatMax * pTHat)
        continue;

      bool bRecoilJet = get<1>(isRecoilJet(jet_det_level, phiTT)) && bTT;

      if (jet_det_level.has_matchedJetGeo()) {

        const auto jetsMatchedPartLevel = jet_det_level.template matchedJetGeo_as<std::decay_t<PartLevelJets>>(); // we can add "matchedJetPt_as" later

        for (const auto& jet_matched_part_level : jetsMatchedPartLevel) {

          /*
          Which histos we want:
          1) det pT vs. part. pT for inclusive jets (corrected for rho*A and not)
          2) det pT vs. part. pT for recoil jets
          3) same as (1) and (2) but 4D with dphi parts
          4) distribution of fake and miss jets
          5) pT and phi resolutions
          */

          spectra.fill(HIST("hJetPt_PartLevel_vs_DetLevel"), jet_det_level.pt(), jet_matched_part_level.pt(), weight);
          spectra.fill(HIST("hJetPt_resolution"), (jet_matched_part_level.pt() - jet_det_level.pt()) / jet_matched_part_level.pt(), jet_matched_part_level.pt(), weight);
          spectra.fill(HIST("hJetPhi_resolution"), jet_matched_part_level.phi() - jet_det_level.phi(), jet_matched_part_level.pt(), weight);

          if (bRecoilJet) {
            spectra.fill(HIST("hJetPt_PartLevel_vs_DetLevel_RecoilJets"), jet_det_level.pt(), jet_matched_part_level.pt(), weight);
            spectra.fill(HIST("hJetPt_resolution_RecoilJets"), (jet_matched_part_level.pt() - jet_det_level.pt()) / jet_matched_part_level.pt(), jet_matched_part_level.pt(), weight);
            spectra.fill(HIST("hJetPhi_resolution_RecoilJets"), jet_matched_part_level.phi() - jet_det_level.phi(), jet_matched_part_level.pt(), weight);
          }
        }
      } else {
        spectra.fill(HIST("hFakeJets_pT"), jet_det_level.pt(), weight);
        if (bRecoilJet)
          spectra.fill(HIST("hFakeJets_pT_RecoilJets"), jet_det_level.pt(), weight);
      }
    }

    // Missed jets
    for (const auto& jet_part_level : jets_part_level) {
      if (!jet_part_level.has_matchedJetGeo()) {
        spectra.fill(HIST("hMissedJets_pT"), jet_part_level.pt(), weight);
      }
    }
  }

  void processData(FilteredColl const& collision,
                   FilteredTracks const& tracks,
                   FilteredJets const& jets)
  {
    if (skipEvent(collision))
      return;

    spectra.fill(HIST("vertexZ"), collision.posZ());
    fillHistograms(collision, jets, tracks);
  }
  PROCESS_SWITCH(RecoilJets, processData, "process data", true);

  void processMCDetLevel(FilteredColl const& collision,
                         FilteredTracks const& tracks,
                         FilteredJetsDetLevel const& jets)
  {
    if (skipEvent(collision))
      return;

    spectra.fill(HIST("vertexZ"), collision.posZ());
    fillHistograms(collision, jets, tracks, true);
  }
  PROCESS_SWITCH(RecoilJets, processMCDetLevel, "process MC detector level", false);

  void processMCDetLevelWeighted(FilteredCollDetLevelGetWeight const& collision,
                                 aod::JetMcCollisions const&,
                                 FilteredTracks const& tracks,
                                 FilteredJetsDetLevel const& jets)
  {
    if (skipEvent(collision))
      return;

    /// \TODO: should we implement function to check whether Collision was reconstructed (has_mcCollision() function)? Example: https://github.com/AliceO2Group/O2Physics/blob/1cba330514ab47c15c0095d8cee9633723d8e2a7/PWGJE/Tasks/v0qa.cxx#L166?
    auto weight = collision.mcCollision().weight();
    spectra.fill(HIST("vertexZ"), collision.posZ(), weight);
    fillHistograms(collision, jets, tracks, true, weight);
  }
  PROCESS_SWITCH(RecoilJets, processMCDetLevelWeighted, "process MC detector level with event weight", false);

  void processMCPartLevel(FilteredCollPartLevel const& collision,
                          aod::JetParticles const& particles,
                          FilteredJetsPartLevel const& jets)
  {
    spectra.fill(HIST("vertexZ"), collision.posZ());
    fillMCPHistograms(collision, jets, particles);
  }
  PROCESS_SWITCH(RecoilJets, processMCPartLevel, "process MC particle level", false);

  void processMCPartLevelWeighted(FilteredCollPartLevel const& collision,
                                  aod::JetParticles const& particles,
                                  FilteredJetsPartLevel const& jets)
  {
    auto weight = collision.weight();
    spectra.fill(HIST("vertexZ"), collision.posZ(), weight);
    fillMCPHistograms(collision, jets, particles, weight);
  }
  PROCESS_SWITCH(RecoilJets, processMCPartLevelWeighted, "process MC particle level with event weight", false);

  void processJetsMatched(FilteredCollDetLevelGetWeight const& collision,
                          aod::JetMcCollisions const&,
                          FilteredTracks const& tracks,
                          FilteredMatchedJetsDetLevel const& mcdjets,
                          FilteredMatchedJetsPartLevel const& mcpjets)
  {
    if (skipEvent(collision))
      return;
    auto mcpjetsPerMCCollision = mcpjets.sliceBy(partJetsPerCollision, collision.mcCollisionId());
    fillMatchedHistograms(tracks, mcdjets, mcpjetsPerMCCollision);
  }
  PROCESS_SWITCH(RecoilJets, processJetsMatched, "process matching of MC jets (no weight)", false);

  void processJetsMatchedWeighted(FilteredCollDetLevelGetWeight const& collision,
                                  aod::JetMcCollisions const&,
                                  FilteredTracks const& tracks,
                                  FilteredMatchedJetsDetLevel const& mcdjets,
                                  FilteredMatchedJetsPartLevel const& mcpjets)
  {
    if (skipEvent(collision))
      return;

    auto mcpjetsPerMCCollision = mcpjets.sliceBy(partJetsPerCollision, collision.mcCollisionId());
    auto weight = collision.mcCollision().weight();

    fillMatchedHistograms(tracks, mcdjets, mcpjetsPerMCCollision, weight);
  }
  PROCESS_SWITCH(RecoilJets, processJetsMatchedWeighted, "process matching of MC jets (weighted)", false);

  //------------------------------------------------------------------------------
  // Auxiliary functions
  template <typename Collision>
  bool skipEvent(const Collision& coll)
  {
    /// \brief: trigger cut is needed for pp data
    return !jetderiveddatautilities::selectCollision(coll, eventSelectionBits) || !jetderiveddatautilities::selectTrigger(coll, triggerMaskBits);
  }

  template <typename Track>
  bool skipTrack(const Track& track)
  {
    return !jetderiveddatautilities::selectTrack(track, trackSelection);
  }

  template <typename Jet>
  std::tuple<double, bool> isRecoilJet(const Jet& jet,
                                       double phiTT)
  {
    double dphi = std::fabs(RecoDecay::constrainAngle(jet.phi() - phiTT, -constants::math::PI));
    return {dphi, (constants::math::PI - recoilRegion) < dphi};
  }

  double getPhiTT(const std::vector<double>& vPhiOfTT)
  {
    auto iTrig = rand->Integer(vPhiOfTT.size());
    return vPhiOfTT[iTrig];
  }

  float getPtHat(float weight)
  {
    return 10. / (std::pow(weight, 1.0 / pTHatExponent));
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<RecoilJets>(cfgc)}; }
