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

// h+jet analysis task
//
// Authors: Daniel Jones

#include <vector>

#include "TRandom3.h"

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

struct hJetAnalysis {

  Configurable<std::string> eventSelections{"eventSelection", "sel8", "choose event selection"};
  Configurable<std::string> trackSelections{"trackSelection", "globalTracks", "set track selections"};
  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum pT acceptance for tracks"};
  Configurable<float> trackPtMax{"trackPtMax", 100.0, "maximum pT acceptance for tracks"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum eta acceptance for tracks"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum eta acceptance for tracks"};
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> pt_TTref_min{"pt_TTref_min", 5, "reference minimum trigger track pt"};
  Configurable<float> pt_TTref_max{"pt_TTref_max", 7, "reference maximum trigger track pt"};
  Configurable<float> pt_TTsig_min{"pt_TTsig_min", 20, "signal minimum trigger track pt"};
  Configurable<float> pt_TTsig_max{"pt_TTsig_max", 50, "signal maximum trigger track pt"};
  Configurable<float> frac_sig{"frac_sig", 0.5, "fraction of events to use for signal"};
  Configurable<float> jetR{"jetR", 0.4, "jet resolution parameter"};
  Configurable<float> pTHatExponent{"pTHatExponent", 6.0, "exponent of the event weight for the calculation of pTHat"};
  Configurable<float> pTHatMaxMCD{"pTHatMaxMCD", 999.0, "maximum fraction of hard scattering for jet acceptance in detector MC"};
  Configurable<float> pTHatMaxMCP{"pTHatMaxMCP", 999.0, "maximum fraction of hard scattering for jet acceptance in particle MC"};

  Preslice<soa::Join<aod::Charged1MCParticleLevelJets, aod::Charged1MCParticleLevelJetConstituents>> PartJetsPerCollision = aod::jet::mcCollisionId;

  TRandom3* rand = new TRandom3(0);

  Filter jetCuts = aod::jet::r == nround(jetR.node() * 100.0f);
  Filter trackCuts = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax);
  Filter eventTrackLevelCuts = nabs(aod::jcollision::posZ) < vertexZCut;

  HistogramRegistry registry{"registry",
                             {{"hNtrig", "number of triggers;trigger type;entries", {HistType::kTH1F, {{2, 0, 2}}}},
                              {"hZvtxSelected", "Z vertex position;Z_{vtx};entries", {HistType::kTH1F, {{80, -20, 20}}}},
                              {"hPtTrack", "Track p_{T};p_{T};entries", {HistType::kTH1F, {{200, 0, 200}}}},
                              {"hEtaTrack", "Track #eta;#eta;entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
                              {"hPhiTrack", "Track #phi;#phi;entries", {HistType::kTH1F, {{160, -1.0, 7.0}}}},
                              {"hReferencePtDPhi", "jet p_{T} vs DPhi;#Delta#phi;p_{T,jet}", {HistType::kTH2F, {{100, 0, 2 * o2::constants::math::PI}, {150, 0, 150}}}},
                              {"hSignalPtDPhi", "jet p_{T} vs DPhi;#Delta#phi;p_{T,jet}", {HistType::kTH2F, {{100, 0, 2 * o2::constants::math::PI}, {150, 0, 150}}}},
                              {"hReferencePt", "jet p_{T};p_{T,jet};entries", {HistType::kTH1F, {{150, 0, 150}}}},
                              {"hSignalPt", "jet p_{T};p_{T,jet};entries", {HistType::kTH1F, {{150, 0, 150}}}},
                              {"hSignalLeadingTrack", "leading track p_{T};p_{T,jet};#Delta#phi;leading track p_{T}", {HistType::kTH3F, {{150, 0, 150}, {100, 0, 2 * o2::constants::math::PI}, {150, 0, 150}}}},
                              {"hReferenceLeadingTrack", "leading track p_{T};p_{T,jet};#Delta#phi;leading track p_{T}", {HistType::kTH3F, {{150, 0, 150}, {100, 0, 2 * o2::constants::math::PI}, {150, 0, 150}}}},
                              {"hJetSignalMultiplicity", "jet multiplicity;N_{jets};entries", {HistType::kTH1F, {{10, 0, 10}}}},
                              {"hJetReferenceMultiplicity", "jet multiplicity;N_{jets};entries", {HistType::kTH1F, {{10, 0, 10}}}},
                              {"hJetSignalConstituentMultiplicity", "jet constituent multiplicity;p_{T,jet};#Delta#phi;N_{constituents}", {HistType::kTH3F, {{150, 0, 150}, {100, 0, 2 * o2::constants::math::PI}, {50, 0, 50}}}},
                              {"hJetReferenceConstituentMultiplicity", "jet constituent multiplicity;p_{T,jet};#Delta#phi;N_{constituents}", {HistType::kTH3F, {{150, 0, 150}, {100, 0, 2 * o2::constants::math::PI}, {50, 0, 50}}}},
                              {"hJetSignalConstituentPt", "jet constituent p_{T};p_{T,jet};#Delta#phi;p_{T,constituent}", {HistType::kTH3F, {{150, 0, 150}, {100, 0, 2 * o2::constants::math::PI}, {150, 0, 150}}}},
                              {"hJetReferenceConstituentPt", "jet constituent p_{T};p_{T,jet};#Delta#phi;p_{T,constituent}", {HistType::kTH3F, {{150, 0, 150}, {100, 0, 2 * o2::constants::math::PI}, {150, 0, 150}}}},
                              {"hSigEventTriggers", "N_{triggers};events", {HistType::kTH1F, {{10, 0, 10}}}},
                              {"hRefEventTriggers", "N_{triggers};events", {HistType::kTH1F, {{10, 0, 10}}}},
                              {"hJetPt", "jet p_{T};p_{T,jet};entries", {HistType::kTH1F, {{200, 0, 200}}}},
                              {"hJetEta", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
                              {"hJetPhi", "jet #phi;#phi_{jet};entries", {HistType::kTH1F, {{160, -1.0, 7.0}}}},
                              {"hPtPart", "Particle p_{T};p_{T};entries", {HistType::kTH1F, {{200, 0, 200}}}},
                              {"hEtaPart", "Particle #eta;#eta;entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
                              {"hPhiPart", "Particle #phi;#phi;entries", {HistType::kTH1F, {{160, -1.0, 7.0}}}},
                              {"hDeltaR", "#DeltaR;#DeltaR;#frac{dN_{jets}}{d#DeltaR}", {HistType::kTH1F, {{50, 0.0, 0.15}}}},
                              {"hDeltaRPart", "Particle #DeltaR;#DeltaR;#frac{1}{N_{jets}}#frac{dN_{jets}}{d#DeltaR}", {HistType::kTH1F, {{50, 0.0, 0.15}}}},
                              {"hDeltaRpT", "jet p_{T} vs #DeltaR;p_{T,jet};#DeltaR", {HistType::kTH2F, {{200, 0, 200}, {50, 0.0, 0.15}}}},
                              {"hDeltaRpTPart", "Particle jet p_{T} vs #DeltaR;p_{T,jet};#DeltaR", {HistType::kTH2F, {{200, 0, 200}, {50, 0.0, 0.15}}}},
                              {"hDeltaRSignal", "#DeltaR;#DeltaR;#frac{dN_{jets}}{d#DeltaR}", {HistType::kTH1F, {{50, 0.0, 0.15}}}},
                              {"hDeltaRSignalPart", "Particle #DeltaR;#DeltaR;#frac{1}{N_{jets}}#frac{dN_{jets}}{d#DeltaR}", {HistType::kTH1F, {{50, 0.0, 0.15}}}},
                              {"hDeltaRpTSignal", "jet p_{T} vs #DeltaR;p_{T,jet};#DeltaR", {HistType::kTH2F, {{200, 0, 200}, {50, 0.0, 0.15}}}},
                              {"hDeltaRpTSignalPart", "Particle jet p_{T} vs #DeltaR;p_{T,jet};#DeltaR", {HistType::kTH2F, {{200, 0, 200}, {50, 0.0, 0.15}}}},
                              {"hDeltaRpTDPhiSignal", "jet p_{T} vs #DeltaR vs #Delta#phi;p_{T,jet};#Delta#phi;#DeltaR", {HistType::kTH3F, {{200, 0, 200}, {100, 0, 2 * o2::constants::math::PI}, {50, 0.0, 0.15}}}},
                              {"hDeltaRpTDPhiSignalPart", "Particle jet p_{T} vs #DeltaR vs #Delta#phi;p_{T,jet};#Delta#phi;#DeltaR", {HistType::kTH3F, {{200, 0, 200}, {100, 0, 2 * o2::constants::math::PI}, {50, 0.0, 0.15}}}},
                              {"hDeltaRReference", "#DeltaR;#DeltaR;#frac{dN_{jets}}{d#DeltaR}", {HistType::kTH1F, {{50, 0.0, 0.15}}}},
                              {"hDeltaRPartReference", "Particle #DeltaR;#DeltaR;#frac{1}{N_{jets}}#frac{dN_{jets}}{d#DeltaR}", {HistType::kTH1F, {{50, 0.0, 0.15}}}},
                              {"hDeltaRpTReference", "jet p_{T} vs #DeltaR;p_{T,jet};#DeltaR", {HistType::kTH2F, {{200, 0, 200}, {50, 0.0, 0.15}}}},
                              {"hDeltaRpTPartReference", "Particle jet p_{T} vs #DeltaR;p_{T,jet};#DeltaR", {HistType::kTH2F, {{200, 0, 200}, {50, 0.0, 0.15}}}},
                              {"hDeltaRpTDPhiReference", "jet p_{T} vs #DeltaR vs #Delta#phi;p_{T,jet};#Delta#phi;#DeltaR", {HistType::kTH3F, {{200, 0, 200}, {100, 0, 2 * o2::constants::math::PI}, {50, 0.0, 0.15}}}},
                              {"hDeltaRpTDPhiReferencePart", "jet p_{T} vs #DeltaR vs #Delta#phi;p_{T,jet};#Delta#phi;#DeltaR", {HistType::kTH3F, {{200, 0, 200}, {100, 0, 2 * o2::constants::math::PI}, {50, 0.0, 0.15}}}},
                              {"hPtMatched", "p_{T} matching;p_{T,det};p_{T,part}", {HistType::kTH2F, {{200, 0, 200}, {200, 0, 200}}}},
                              {"hPhiMatched", "#phi matching;#phi_{det};#phi_{part}", {HistType::kTH2F, {{160, -1.0, 7.0}, {160, -1.0, 7.0}}}},
                              {"hDeltaRMatched", "#DeltaR matching;#DeltaR_{det};#DeltaR_{part}", {HistType::kTH2F, {{50, 0.0, 0.15}, {50, 0.0, 0.15}}}},
                              {"hPtResolution", "p_{T} resolution;p_{T,part};Relative Resolution", {HistType::kTH2F, {{200, 0, 200}, {1000, -5.0, 5.0}}}},
                              {"hPhiResolution", "#phi resolution;#p{T,part};Resolution", {HistType::kTH2F, {{200, 0, 200}, {1000, -7.0, 7.0}}}},
                              {"hDeltaRResolution", "#DeltaR Resolution;p_{T,part};Resolution", {HistType::kTH2F, {{200, 0, 200}, {1000, -0.15, 0.15}}}},
                              {"hFullMatching", "Full 6D matching;p_{T,det};p_{T,part};#phi_{det};#phi_{part};#DeltaR_{det};#DeltaR_{part}", {HistType::kTHnSparseD, {{200, 0, 200}, {200, 0, 200}, {160, -1.0, 7.0}, {160, -1.0, 7.0}, {50, 0.0, 0.15}, {50, 0.0, 0.15}}}}}};

  int eventSelection = -1;
  int trackSelection = -1;

  Service<o2::framework::O2DatabasePDG> pdg;

  void init(InitContext const&)
  {
    eventSelection = jetderiveddatautilities::initialiseEventSelection(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

    Filter jetCuts = aod::jet::r == nround(jetR.node() * 100.0f);
    Filter trackCuts = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax);
    Filter eventTrackLevelCuts = nabs(aod::jcollision::posZ) < vertexZCut;
  }

  template <typename T, typename U, typename W>
  void fillHistograms(T const& jets, W const& /*jetsWTA*/, U const& tracks, float weight = 1.0)
  {
    bool is_sig_col;
    std::vector<double> phi_TT_ar;
    double phi_TT = 0;
    int trig_number = 0;
    int n_TT = 0;
    double leadingPT = 0;
    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));

    float dice = rand->Rndm();
    if (dice < frac_sig)
      is_sig_col = true;
    else
      is_sig_col = false;

    for (auto& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
      if (is_sig_col && track.pt() < pt_TTsig_max && track.pt() > pt_TTsig_min) {
        phi_TT_ar.push_back(track.phi());
        n_TT++;
      }
      if (!is_sig_col && track.pt() < pt_TTref_max && track.pt() > pt_TTref_min) {
        phi_TT_ar.push_back(track.phi());
        n_TT++;
      }
      registry.fill(HIST("hPtTrack"), track.pt(), weight);
      registry.fill(HIST("hEtaTrack"), track.eta(), weight);
      registry.fill(HIST("hPhiTrack"), track.phi(), weight);
    }

    if (n_TT > 0) {
      trig_number = rand->Integer(n_TT);
      phi_TT = phi_TT_ar[trig_number];
      if (is_sig_col) {
        registry.fill(HIST("hNtrig"), 1.5, weight);
        registry.fill(HIST("hJetSignalMultiplicity"), jets.size(), weight);
        registry.fill(HIST("hSigEventTriggers"), n_TT, weight);
      }
      if (!is_sig_col) {
        registry.fill(HIST("hNtrig"), 0.5, weight);
        registry.fill(HIST("hJetReferenceMultiplicity"), jets.size(), weight);
        registry.fill(HIST("hRefEventTriggers"), n_TT, weight);
      }
    }

    for (auto& jet : jets) {
      if (jet.pt() > pTHatMaxMCD * pTHat) {
        continue;
      }
      registry.fill(HIST("hJetPt"), jet.pt(), weight);
      registry.fill(HIST("hJetEta"), jet.eta(), weight);
      registry.fill(HIST("hJetPhi"), jet.phi(), weight);
      for (auto& jetWTA : jet.template matchedJetGeo_as<std::decay_t<W>>()) {
        double deltaPhi = RecoDecay::constrainAngle(jetWTA.phi() - jet.phi(), -o2::constants::math::PI);
        double deltaEta = jetWTA.eta() - jet.eta();
        double dR = RecoDecay::sqrtSumOfSquares(deltaPhi, deltaEta);
        registry.fill(HIST("hDeltaR"), dR, weight);
        registry.fill(HIST("hDeltaRpT"), jet.pt(), dR, weight);
      }
      if (n_TT > 0) {
        float dphi = RecoDecay::constrainAngle(jet.phi() - phi_TT);
        if (is_sig_col) {
          for (auto& jetWTA : jet.template matchedJetGeo_as<std::decay_t<W>>()) {
            double deltaPhi = RecoDecay::constrainAngle(jetWTA.phi() - jet.phi(), -o2::constants::math::PI);
            double deltaEta = jetWTA.eta() - jet.eta();
            double dR = RecoDecay::sqrtSumOfSquares(deltaPhi, deltaEta);
            if (std::abs(dphi - o2::constants::math::PI) < 0.6) {
              registry.fill(HIST("hDeltaRpTSignal"), jet.pt(), dR, weight);
              registry.fill(HIST("hDeltaRSignal"), dR, weight);
            }
            registry.fill(HIST("hDeltaRpTDPhiSignal"), jet.pt(), dphi, dR, weight);
          }
          registry.fill(HIST("hSignalPtDPhi"), dphi, jet.pt(), weight);
          if (std::abs(dphi - o2::constants::math::PI) < 0.6) {
            registry.fill(HIST("hSignalPt"), jet.pt(), weight);
          }
          registry.fill(HIST("hJetSignalConstituentMultiplicity"), jet.pt(), dphi, jet.tracksIds().size(), weight);
          for (auto& constituent : jet.template tracks_as<U>()) {
            if (constituent.pt() > leadingPT) {
              leadingPT = constituent.pt();
            }
            registry.fill(HIST("hJetSignalConstituentPt"), jet.pt(), dphi, constituent.pt(), weight);
          }
          registry.fill(HIST("hSignalLeadingTrack"), jet.pt(), dphi, leadingPT, weight);
        }
        if (!is_sig_col) {
          for (auto& jetWTA : jet.template matchedJetGeo_as<std::decay_t<W>>()) {
            double deltaPhi = RecoDecay::constrainAngle(jetWTA.phi() - jet.phi(), -o2::constants::math::PI);
            double deltaEta = jetWTA.eta() - jet.eta();
            double dR = RecoDecay::sqrtSumOfSquares(deltaPhi, deltaEta);
            if (std::abs(dphi - o2::constants::math::PI) < 0.6) {
              registry.fill(HIST("hDeltaRpTReference"), jet.pt(), dR, weight);
              registry.fill(HIST("hDeltaRReference"), dR, weight);
            }
            registry.fill(HIST("hDeltaRpTDPhiReference"), jet.pt(), dphi, dR, weight);
          }
          registry.fill(HIST("hReferencePtDPhi"), dphi, jet.pt(), weight);
          if (std::abs(dphi - o2::constants::math::PI) < 0.6) {
            registry.fill(HIST("hReferencePt"), jet.pt(), weight);
          }
          registry.fill(HIST("hJetReferenceConstituentMultiplicity"), jet.pt(), dphi, jet.tracksIds().size(), weight);
          for (auto& constituent : jet.template tracks_as<U>()) {
            if (constituent.pt() > leadingPT) {
              leadingPT = constituent.pt();
            }
            registry.fill(HIST("hJetReferenceConstituentPt"), jet.pt(), dphi, constituent.pt(), weight);
          }
          registry.fill(HIST("hReferenceLeadingTrack"), jet.pt(), dphi, leadingPT, weight);
        }
      }
    }
  }

  template <typename T, typename W, typename U>
  void fillMCPHistograms(T const& jets, W const& /*jetsWTA*/, U const& particles, float weight = 1.0)
  {
    bool is_sig_col;
    std::vector<double> phi_TT_ar;
    double phi_TT = 0;
    int trig_number = 0;
    int n_TT = 0;
    double leadingPT = 0;
    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));

    float dice = rand->Rndm();
    if (dice < frac_sig)
      is_sig_col = true;
    else
      is_sig_col = false;

    for (auto& particle : particles) {
      auto pdgParticle = pdg->GetParticle(particle.pdgCode());
      if (!pdgParticle) {
        continue;
      }
      if ((pdgParticle->Charge() == 0.0) || (!particle.isPhysicalPrimary())) {
        continue;
      }
      if (is_sig_col && particle.pt() < pt_TTsig_max && particle.pt() > pt_TTsig_min) {
        phi_TT_ar.push_back(particle.phi());
        n_TT++;
      }
      if (!is_sig_col && particle.pt() < pt_TTref_max && particle.pt() > pt_TTref_min) {
        phi_TT_ar.push_back(particle.phi());
        n_TT++;
      }
      registry.fill(HIST("hPtPart"), particle.pt(), weight);
      registry.fill(HIST("hEtaPart"), particle.eta(), weight);
      registry.fill(HIST("hPhiPart"), particle.phi(), weight);
    }

    if (n_TT > 0) {
      trig_number = rand->Integer(n_TT);
      phi_TT = phi_TT_ar[trig_number];
      if (is_sig_col) {
        registry.fill(HIST("hNtrig"), 1.5, weight);
        registry.fill(HIST("hJetSignalMultiplicity"), jets.size(), weight);
        registry.fill(HIST("hSigEventTriggers"), n_TT, weight);
      }
      if (!is_sig_col) {
        registry.fill(HIST("hNtrig"), 0.5, weight);
        registry.fill(HIST("hJetReferenceMultiplicity"), jets.size(), weight);
        registry.fill(HIST("hRefEventTriggers"), n_TT, weight);
      }
    }

    for (auto& jet : jets) {
      if (jet.pt() > pTHatMaxMCP * pTHat) {
        continue;
      }
      registry.fill(HIST("hJetPt"), jet.pt(), weight);
      registry.fill(HIST("hJetEta"), jet.eta(), weight);
      registry.fill(HIST("hJetPhi"), jet.phi(), weight);
      for (auto& jetWTA : jet.template matchedJetGeo_as<std::decay_t<W>>()) {
        double deltaPhi = RecoDecay::constrainAngle(jetWTA.phi() - jet.phi(), -o2::constants::math::PI);
        double deltaEta = jetWTA.eta() - jet.eta();
        double dR = RecoDecay::sqrtSumOfSquares(deltaPhi, deltaEta);
        registry.fill(HIST("hDeltaRPart"), dR, weight);
        registry.fill(HIST("hDeltaRpTPart"), jet.pt(), dR, weight);
      }
      if (n_TT > 0) {
        float dphi = RecoDecay::constrainAngle(jet.phi() - phi_TT);
        if (is_sig_col) {
          for (auto& jetWTA : jet.template matchedJetGeo_as<std::decay_t<W>>()) {
            double deltaPhi = RecoDecay::constrainAngle(jetWTA.phi() - jet.phi(), -o2::constants::math::PI);
            double deltaEta = jetWTA.eta() - jet.eta();
            double dR = RecoDecay::sqrtSumOfSquares(deltaPhi, deltaEta);
            if (std::abs(dphi - o2::constants::math::PI) < 0.6) {
              registry.fill(HIST("hDeltaRpTSignalPart"), jet.pt(), dR, weight);
              registry.fill(HIST("hDeltaRSignalPart"), dR, weight);
            }
            registry.fill(HIST("hDeltaRpTDPhiSignalPart"), jet.pt(), dphi, dR, weight);
          }
          registry.fill(HIST("hSignalPtDPhi"), dphi, jet.pt(), weight);
          if (std::abs(dphi - o2::constants::math::PI) < 0.6) {
            registry.fill(HIST("hSignalPt"), jet.pt(), weight);
          }
          registry.fill(HIST("hJetSignalConstituentMultiplicity"), jet.pt(), dphi, jet.tracksIds().size(), weight);
          for (auto& constituent : jet.template tracks_as<U>()) {
            if (constituent.pt() > leadingPT) {
              leadingPT = constituent.pt();
            }
            registry.fill(HIST("hJetSignalConstituentPt"), jet.pt(), dphi, constituent.pt(), weight);
          }
          registry.fill(HIST("hSignalLeadingTrack"), jet.pt(), dphi, leadingPT, weight);
        }
        if (!is_sig_col) {
          for (auto& jetWTA : jet.template matchedJetGeo_as<std::decay_t<W>>()) {
            double deltaPhi = RecoDecay::constrainAngle(jetWTA.phi() - jet.phi(), -o2::constants::math::PI);
            double deltaEta = jetWTA.eta() - jet.eta();
            double dR = RecoDecay::sqrtSumOfSquares(deltaPhi, deltaEta);
            if (std::abs(dphi - o2::constants::math::PI) < 0.6) {
              registry.fill(HIST("hDeltaRpTPartReference"), jet.pt(), dR, weight);
              registry.fill(HIST("hDeltaRPartReference"), dR, weight);
            }
            registry.fill(HIST("hDeltaRpTDPhiReferencePart"), jet.pt(), dphi, dR, weight);
          }
          registry.fill(HIST("hReferencePtDPhi"), dphi, jet.pt(), weight);
          if (std::abs(dphi - o2::constants::math::PI) < 0.6) {
            registry.fill(HIST("hReferencePt"), jet.pt(), weight);
          }
          registry.fill(HIST("hJetReferenceConstituentMultiplicity"), jet.pt(), dphi, jet.tracksIds().size(), weight);
          for (auto& constituent : jet.template tracks_as<U>()) {
            if (constituent.pt() > leadingPT) {
              leadingPT = constituent.pt();
            }
            registry.fill(HIST("hJetReferenceConstituentPt"), jet.pt(), dphi, constituent.pt(), weight);
          }
          registry.fill(HIST("hReferenceLeadingTrack"), jet.pt(), dphi, leadingPT, weight);
        }
      }
    }
  }

  template <typename T, typename V, typename W, typename U>
  void fillMatchedHistograms(T const& jetBase, V const& mcdjetsWTA, W const& mcpjetsWTA, U const&, float weight = 1.0)
  {
    double dR = 0;
    double dRp = 0;

    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    if (jetBase.pt() > pTHatMaxMCD * pTHat) {
      return;
    }

    for (const auto& mcdjetWTA : mcdjetsWTA) {
      double djet = RecoDecay::sqrtSumOfSquares(RecoDecay::constrainAngle(jetBase.phi() - mcdjetWTA.phi(), -o2::constants::math::PI), jetBase.eta() - mcdjetWTA.eta());
      if (mcdjetWTA.pt() > pTHatMaxMCD * pTHat) {
        continue;
      }
      if (djet < 0.6 * jetR) {
        dR = djet;
        break;
      }
    }

    if (jetBase.has_matchedJetGeo()) {
      for (auto& jetTag : jetBase.template matchedJetGeo_as<std::decay_t<U>>()) {
        if (jetTag.pt() > pTHatMaxMCP * pTHat) {
          continue;
        }
        for (const auto& mcpjetWTA : mcpjetsWTA) {
          double djetp = RecoDecay::sqrtSumOfSquares(RecoDecay::constrainAngle(jetTag.phi() - mcpjetWTA.phi(), -o2::constants::math::PI), jetTag.eta() - mcpjetWTA.eta());
          if (mcpjetWTA.pt() > pTHatMaxMCP * pTHat) {
            continue;
          }
          if (djetp < 0.6 * jetR) {
            dRp = djetp;
            break;
          }
        }
        registry.fill(HIST("hPtMatched"), jetBase.pt(), jetTag.pt(), weight);
        registry.fill(HIST("hPhiMatched"), jetBase.phi(), jetTag.phi(), weight);
        registry.fill(HIST("hPtResolution"), jetTag.pt(), (jetTag.pt() - jetBase.pt()) / jetTag.pt(), weight);
        registry.fill(HIST("hPhiResolution"), jetTag.pt(), jetTag.phi() - jetBase.phi(), weight);
        registry.fill(HIST("hDeltaRMatched"), dR, dRp, weight);
        registry.fill(HIST("hDeltaRResolution"), jetTag.pt(), dRp - dR, weight);
        registry.fill(HIST("hFullMatching"), jetBase.pt(), jetTag.pt(), jetBase.phi(), jetTag.phi(), dR, dRp, weight);
      }
    }
  }

  void processData(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                   soa::Filtered<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents, aod::ChargedJetsMatchedToCharged1Jets>> const& jets,
                   soa::Filtered<soa::Join<aod::Charged1Jets, aod::Charged1JetConstituents, aod::Charged1JetsMatchedToChargedJets>> const& jetsWTA,
                   soa::Filtered<aod::JetTracks> const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }
    registry.fill(HIST("hZvtxSelected"), collision.posZ());
    fillHistograms(jets, jetsWTA, tracks);
  }
  PROCESS_SWITCH(hJetAnalysis, processData, "process data", true);

  void processMCD(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                  soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToCharged1MCDetectorLevelJets>> const& jets,
                  soa::Filtered<soa::Join<aod::Charged1MCDetectorLevelJets, aod::Charged1MCDetectorLevelJetConstituents, aod::Charged1MCDetectorLevelJetsMatchedToChargedMCDetectorLevelJets>> const& jetsWTA,
                  soa::Filtered<aod::JetTracks> const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }
    registry.fill(HIST("hZvtxSelected"), collision.posZ());
    fillHistograms(jets, jetsWTA, tracks);
  }
  PROCESS_SWITCH(hJetAnalysis, processMCD, "process MC detector level", false);

  void processMCDWeighted(soa::Filtered<soa::Join<aod::JetCollisions, aod::JMcCollisionLbs>>::iterator const& collision,
                          aod::JMcCollisions const&,
                          soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToCharged1MCDetectorLevelJets>> const& jets,
                          soa::Filtered<soa::Join<aod::Charged1MCDetectorLevelJets, aod::Charged1MCDetectorLevelJetConstituents, aod::Charged1MCDetectorLevelJetsMatchedToChargedMCDetectorLevelJets>> const& jetsWTA,
                          soa::Filtered<aod::JetTracks> const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }
    registry.fill(HIST("hZvtxSelected"), collision.posZ(), collision.mcCollision().weight());
    fillHistograms(jets, jetsWTA, tracks, collision.mcCollision().weight());
  }
  PROCESS_SWITCH(hJetAnalysis, processMCDWeighted, "process MC detector level with event weights", false);

  void processMCP(aod::JetMcCollision const& collision,
                  soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToCharged1MCParticleLevelJets>> const& jets,
                  soa::Filtered<soa::Join<aod::Charged1MCParticleLevelJets, aod::Charged1MCParticleLevelJetConstituents, aod::Charged1MCParticleLevelJetsMatchedToChargedMCParticleLevelJets>> const& jetsWTA,
                  aod::JetParticles const& particles)
  {
    if (std::abs(collision.posZ()) > vertexZCut) {
      return;
    }
    registry.fill(HIST("hZvtxSelected"), collision.posZ());
    fillMCPHistograms(jets, jetsWTA, particles);
  }
  PROCESS_SWITCH(hJetAnalysis, processMCP, "process MC particle level", false);

  void processMCPWeighted(aod::JetMcCollision const& collision,
                          soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToCharged1MCParticleLevelJets>> const& jets,
                          soa::Filtered<soa::Join<aod::Charged1MCParticleLevelJets, aod::Charged1MCParticleLevelJetConstituents, aod::Charged1MCParticleLevelJetsMatchedToChargedMCParticleLevelJets>> const& jetsWTA,
                          aod::JetParticles const& particles)
  {
    if (std::abs(collision.posZ()) > vertexZCut) {
      return;
    }
    registry.fill(HIST("hZvtxSelected"), collision.posZ(), collision.weight());
    fillMCPHistograms(jets, jetsWTA, particles, collision.weight());
  }
  PROCESS_SWITCH(hJetAnalysis, processMCPWeighted, "process MC particle level with event weights", false);

  void processJetsMCPMCDMatched(soa::Filtered<aod::JetCollisionsMCD>::iterator const& collision,
                                soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>> const& mcdjets,
                                soa::Filtered<soa::Join<aod::Charged1MCDetectorLevelJets, aod::Charged1MCDetectorLevelJetConstituents>> const& mcdjetsWTA,
                                soa::Filtered<soa::Join<aod::Charged1MCParticleLevelJets, aod::Charged1MCParticleLevelJetConstituents>> const& mcpjetsWTA,
                                aod::JetTracks const&,
                                aod::JetParticles const&,
                                aod::JetMcCollisions const&,
                                soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>> const& mcpjets)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }
    registry.fill(HIST("hZvtxSelected"), collision.posZ());
    const auto& mcpjetsWTACut = mcpjetsWTA.sliceBy(PartJetsPerCollision, collision.mcCollisionId());
    for (const auto& mcdjet : mcdjets) {
      fillMatchedHistograms(mcdjet, mcdjetsWTA, mcpjetsWTACut, mcpjets);
    }
  }
  PROCESS_SWITCH(hJetAnalysis, processJetsMCPMCDMatched, "process MC matched (inc jets)", false);

  void processJetsMCPMCDMatchedWeighted(soa::Filtered<aod::JetCollisionsMCD>::iterator const& collision,
                                        soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets, aod::ChargedMCDetectorLevelJetEventWeights>> const& mcdjets,
                                        soa::Filtered<soa::Join<aod::Charged1MCDetectorLevelJets, aod::Charged1MCDetectorLevelJetConstituents>> const& mcdjetsWTA,
                                        soa::Filtered<soa::Join<aod::Charged1MCParticleLevelJets, aod::Charged1MCParticleLevelJetConstituents>> const& mcpjetsWTA,
                                        aod::JetTracks const&,
                                        aod::JetParticles const&,
                                        aod::JetMcCollisions const&,
                                        soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets, aod::ChargedMCParticleLevelJetEventWeights>> const& mcpjets)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }
    registry.fill(HIST("hZvtxSelected"), collision.posZ());
    const auto& mcpjetsWTACut = mcpjetsWTA.sliceBy(PartJetsPerCollision, collision.mcCollisionId());
    for (const auto& mcdjet : mcdjets) {
      fillMatchedHistograms(mcdjet, mcdjetsWTA, mcpjetsWTACut, mcpjets, mcdjet.eventWeight());
    }
  }
  PROCESS_SWITCH(hJetAnalysis, processJetsMCPMCDMatchedWeighted, "process MC matched with event weights (inc jets)", false);

  void processRecoilJetsMCPMCDMatched(soa::Filtered<aod::JetCollisionsMCD>::iterator const& collision,
                                      soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>> const& mcdjets,
                                      soa::Filtered<soa::Join<aod::Charged1MCDetectorLevelJets, aod::Charged1MCDetectorLevelJetConstituents>> const& mcdjetsWTA,
                                      soa::Filtered<soa::Join<aod::Charged1MCParticleLevelJets, aod::Charged1MCParticleLevelJetConstituents>> const& mcpjetsWTA,
                                      aod::JetTracks const& tracks,
                                      aod::JetParticles const&,
                                      aod::JetMcCollisions const&,
                                      soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>> const& mcpjets)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }
    registry.fill(HIST("hZvtxSelected"), collision.posZ());
    const auto& mcpjetsWTACut = mcpjetsWTA.sliceBy(PartJetsPerCollision, collision.mcCollisionId());
    bool ishJetEvent = false;
    for (auto& track : tracks) {
      if (track.pt() < pt_TTsig_max && track.pt() > pt_TTsig_min) {
        ishJetEvent = true;
        break;
      }
    }
    if (ishJetEvent) {
      for (const auto& mcdjet : mcdjets) {
        fillMatchedHistograms(mcdjet, mcdjetsWTA, mcpjetsWTACut, mcpjets);
      }
    }
  }
  PROCESS_SWITCH(hJetAnalysis, processRecoilJetsMCPMCDMatched, "process MC matched (recoil jets)", false);

  void processRecoilJetsMCPMCDMatchedWeighted(soa::Filtered<aod::JetCollisionsMCD>::iterator const& collision,
                                              soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets, aod::ChargedMCDetectorLevelJetEventWeights>> const& mcdjets,
                                              soa::Filtered<soa::Join<aod::Charged1MCDetectorLevelJets, aod::Charged1MCDetectorLevelJetConstituents>> const& mcdjetsWTA,
                                              soa::Filtered<soa::Join<aod::Charged1MCParticleLevelJets, aod::Charged1MCParticleLevelJetConstituents>> const& mcpjetsWTA,
                                              aod::JetTracks const& tracks,
                                              aod::JetParticles const&,
                                              aod::JetMcCollisions const&,
                                              soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets, aod::ChargedMCParticleLevelJetEventWeights>> const& mcpjets)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }
    registry.fill(HIST("hZvtxSelected"), collision.posZ());
    const auto& mcpjetsWTACut = mcpjetsWTA.sliceBy(PartJetsPerCollision, collision.mcCollisionId());
    bool ishJetEvent = false;
    for (auto& track : tracks) {
      if (track.pt() < pt_TTsig_max && track.pt() > pt_TTsig_min) {
        ishJetEvent = true;
        break;
      }
    }
    if (ishJetEvent) {
      for (const auto& mcdjet : mcdjets) {
        fillMatchedHistograms(mcdjet, mcdjetsWTA, mcpjetsWTACut, mcpjets, mcdjet.eventWeight());
      }
    }
  }
  PROCESS_SWITCH(hJetAnalysis, processRecoilJetsMCPMCDMatchedWeighted, "process MC matched with event weights (recoil jets)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<hJetAnalysis>(cfgc, TaskName{"hJetAnalysis"})}; }
