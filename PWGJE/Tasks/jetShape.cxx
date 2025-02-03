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

/// \file   jetShape.cxx
/// \author Yuto Nishida <yuto.nishida@cern.ch>
/// \brief Task for measuring the dependence of the jet shape function rho(r) on the distance r from the jet axis.

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"

#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetShapeTask {
  HistogramRegistry registry{"registry",
                             {{"tpcPi", "tpcPi", {HistType::kTH2F, {{1000, 0, 5}, {401, -10.025f, 10.025f}}}},
                              {"tofPi", "tofPi", {HistType::kTH2F, {{1000, 0, 5}, {401, -10.025f, 10.025f}}}},
                              {"tpcPr", "tpcPr", {HistType::kTH2F, {{1000, 0, 5}, {401, -10.025f, 10.025f}}}},
                              {"tofPr", "tofPr", {HistType::kTH2F, {{1000, 0, 5}, {401, -10.025f, 10.025f}}}},
                              {"tofTpc", "tofTpc", {HistType::kTH2F, {{401, -10.025f, 10.025f}, {401, -10.025f, 10.025f}}}},
                              {"tpcDedx", "tpcDedx", {HistType::kTH2F, {{1000, 0, 5}, {1000, 0, 1000}}}},
                              {"tofBeta", "tofBeta", {HistType::kTH2F, {{1000, 0, 5}, {900, 0.2, 1.1}}}},
                              {"tofMass", "tofMass", {HistType::kTH1F, {{3000, 0, 3}}}},
                              {"jetPt", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"jetEta", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
                              {"jetPhi", "jet #phi;#phi_{jet};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
                              {"area", "area", {HistType::kTH1F, {{200, 0, 4}}}},
                              {"rho", "rho", {HistType::kTH1F, {{200, -1, 119}}}},
                              {"ptCorr", "Corrected jet pT; p_{T}^{corr} (GeV/c); Counts", {HistType::kTH1F, {{200, 0, 200}}}},
                              {"ptCorrVsDistance", "ptcorr_vs_distance", {HistType::kTH2F, {{70, 0, 0.7}, {100, 0, 100}}}},
                              {"trackPtVsDistance", "trackpt_vs_distance", {HistType::kTH2F, {{70, 0, 0.7}, {100, 0, 100}}}},
                              {"ptSum", "ptSum", {HistType::kTH2F, {{14, 0, 0.7}, {300, 0, 300}}}},
                              {"ptSumBg1", "ptSumBg1", {HistType::kTH2F, {{14, 0, 0.7}, {300, 0, 300}}}},
                              {"ptSumBg2", "ptSumBg2", {HistType::kTH2F, {{14, 0, 0.7}, {300, 0, 300}}}},
                              {"ptVsCentrality", "ptvscentrality", {HistType::kTH2F, {{100, 0, 100}, {300, 0, 300}}}}}};

  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};

  Configurable<float> jetPtMin{"jetPtMin", 5.0, "minimum jet pT cut"};
  Configurable<float> jetR{"jetR", 0.4, "jet resolution parameter"};

  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

  Configurable<float> jetAreaFractionMin{"jetAreaFractionMin", -99.0, "used to make a cut on the jet areas"};
  Configurable<float> leadingConstituentPtMin{"leadingConstituentPtMin", 5.0, "minimum pT selection on jet constituent"};
  Configurable<float> leadingConstituentPtMax{"leadingConstituentPtMax", 9999.0, "maximum pT selection on jet constituent"};

  Configurable<float> kappa{"kappa", 1.0, "angularity kappa"};
  Configurable<float> alpha{"alpha", 1.0, "angularity alpha"};

  Configurable<std::string> triggerMasks{"triggerMasks", "", "possible JE Trigger masks: fJetChLowPt,fJetChHighPt,fTrackLowPt,fTrackHighPt,fJetD0ChLowPt,fJetD0ChHighPt,fJetLcChLowPt,fJetLcChHighPt,fEMCALReadout,fJetFullHighPt,fJetFullLowPt,fJetNeutralHighPt,fJetNeutralLowPt,fGammaVeryHighPtEMCAL,fGammaVeryHighPtDCAL,fGammaHighPtEMCAL,fGammaHighPtDCAL,fGammaLowPtEMCAL,fGammaLowPtDCAL,fGammaVeryLowPtEMCAL,fGammaVeryLowPtDCAL"};

  int eventSelection = -1;
  int trackSelection = -1;
  std::vector<int> triggerMaskBits;

  void init(o2::framework::InitContext&)
  {
    eventSelection = jetderiveddatautilities::initialiseEventSelection(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
    triggerMaskBits = jetderiveddatautilities::initialiseTriggerMaskBits(triggerMasks);
  }

  template <typename T, typename U>
  bool isAcceptedJet(U const& jet)
  {
    if (jetAreaFractionMin > -98.0) {
      if (jet.area() < jetAreaFractionMin * o2::constants::math::PI * (jet.r() / 100.0) * (jet.r() / 100.0)) {
        return false;
      }
      if (jet.area() < o2::constants::math::PIHALF * (jet.r() / 100.0) * (jet.r() / 100.0)) {
        return false;
      }
    }
    bool checkConstituentPt = true;
    bool checkConstituentMinPt = (leadingConstituentPtMin > 5);
    bool checkConstituentMaxPt = (leadingConstituentPtMax < 9998.0);
    if (!checkConstituentMinPt && !checkConstituentMaxPt) {
      checkConstituentPt = false;
    }

    if (checkConstituentPt) {
      bool isMinLeadingConstituent = !checkConstituentMinPt;
      bool isMaxLeadingConstituent = true;

      for (const auto& constituent : jet.template tracks_as<T>()) {
        double pt = constituent.pt();

        if (checkConstituentMinPt && pt >= leadingConstituentPtMin) {
          isMinLeadingConstituent = true;
        }
        if (checkConstituentMaxPt && pt > leadingConstituentPtMax) {
          isMaxLeadingConstituent = false;
        }
      }
      return isMinLeadingConstituent && isMaxLeadingConstituent;
    }

    return true;
  }

  Filter jetCuts = aod::jet::pt > jetPtMin&& aod::jet::r == nround(jetR.node() * 100.0f);
  Filter collisionFilter = nabs(aod::jcollision::posZ) < vertexZCut;
  Filter mcCollisionFilter = nabs(aod::jmccollision::posZ) < vertexZCut;

  Preslice<soa::Filtered<aod::ChargedMCParticleLevelJets>> perMcCollisionJets = aod::jet::mcCollisionId;

  void processCharged(soa::Filtered<soa::Join<JetCollisions, aod::BkgChargedRhos>>::iterator const& collision, soa::Join<JetTracks, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::TracksExtra, aod::TracksDCA, aod::pidTOFbeta, aod::pidTOFmass> const& tracks, soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }
    double ptDensity[14] = {0};
    double ptDensityBg1[14] = {0};
    double ptDensityBg2[14] = {0};

    // std::cout << collision.centrality() << std::endl;

    for (auto const& jet : jets) {
      if (!isAcceptedJet<JetTracks>(jet)) {
        continue;
      }

      float ptCorr = jet.pt() - collision.rho() * jet.area();

      int nJet = 0;
      nJet += jets.size();

      for (const auto& track : tracks) {
        double preDeltaPhi1 = 0;
        double deltaPhi1 = 0;
        double deltaEta = 0;

        preDeltaPhi1 = track.phi() - jet.phi();

        if (preDeltaPhi1 > o2::constants::math::PI) {
          deltaPhi1 = preDeltaPhi1 - o2::constants::math::TWO_PI;
        } else if (preDeltaPhi1 < -o2::constants::math::PI) {
          deltaPhi1 = preDeltaPhi1 + o2::constants::math::TWO_PI;
        } else {
          deltaPhi1 = preDeltaPhi1;
        }

        float distance = std::sqrt(deltaEta * deltaEta + deltaPhi1 * deltaPhi1);
        registry.fill(HIST("trackPtVsDistance"), distance, track.pt());
        registry.fill(HIST("ptCorrVsDistance"), distance, ptCorr);

        registry.fill(HIST("ptVsCentrality"), collision.centrality(), track.pt());

        double trackPtSum[14] = {0};
        double trackPtSumBg1[14] = {0};
        double trackPtSumBg2[14] = {0};
        double phi2 = jet.phi() + (o2::constants::math::PIHALF);
        double phi3 = jet.phi() - (o2::constants::math::PIHALF);
        double preDeltaPhi2 = 0;
        double deltaPhi2 = 0;
        double preDeltaPhi3 = 0;
        double deltaPhi3 = 0;

        preDeltaPhi2 = track.phi() - phi2;
        preDeltaPhi3 = track.phi() - phi3;

        if (preDeltaPhi2 > o2::constants::math::PI) {
          deltaPhi2 = preDeltaPhi2 - o2::constants::math::TWO_PI;
        } else if (preDeltaPhi2 < -o2::constants::math::PI) {
          deltaPhi2 = preDeltaPhi2 + o2::constants::math::TWO_PI;
        } else {
          deltaPhi2 = preDeltaPhi2;
        }

        if (preDeltaPhi3 > o2::constants::math::PI) {
          deltaPhi3 = preDeltaPhi3 - o2::constants::math::TWO_PI;
        } else if (preDeltaPhi3 < -o2::constants::math::PI) {
          deltaPhi3 = preDeltaPhi3 + o2::constants::math::TWO_PI;
        } else {
          deltaPhi3 = preDeltaPhi3;
        }

        deltaEta = track.eta() - jet.eta();

        float distanceBg1 = std::sqrt(deltaEta * deltaEta + deltaPhi2 * deltaPhi2);
        float distanceBg2 = std::sqrt(deltaEta * deltaEta + deltaPhi3 * deltaPhi3);

        float trackMomentumCategory[] = {0.05f, 0.1f, 0.15f, 0.2f, 0.25f, 0.3f, 0.35f, 0.4f, 0.45f, 0.5f, 0.55f, 0.6f, 0.65f, 0.7f};

        for (int i = 0; i < 14; i++) {
          if (distance < trackMomentumCategory[i])
            trackPtSum[i] += track.pt();
        }

        for (int8_t i = 0; i < 14; i++) {
          ptDensity[i] += trackPtSum[i] / (0.05 * ptCorr);
        }

        for (int i = 0; i < 14; i++) {
          if (distanceBg1 < trackMomentumCategory[i])
            trackPtSumBg1[i] += track.pt();
        }

        for (int i = 0; i < 14; i++) {
          ptDensityBg1[i] += trackPtSumBg1[i] / (0.05 * ptCorr);
        }

        for (int i = 0; i < 14; i++) {
          if (distanceBg2 < trackMomentumCategory[i])
            trackPtSumBg2[i] += track.pt();
        }

        for (int i = 0; i < 14; i++) {
          ptDensityBg2[i] += trackPtSumBg2[i] / (0.05 * ptCorr);
        }

        registry.fill(HIST("tpcDedx"), track.pt(), track.tpcSignal());
        registry.fill(HIST("tofBeta"), track.pt(), track.beta());
        registry.fill(HIST("tofMass"), track.mass());
        registry.fill(HIST("tpcPi"), track.pt(), track.tpcNSigmaPi());
        registry.fill(HIST("tofPi"), track.pt(), track.tofNSigmaPi());
        registry.fill(HIST("tpcPr"), track.pt(), track.tpcNSigmaPr());
        registry.fill(HIST("tofPr"), track.pt(), track.tofNSigmaPr());
      }

      registry.fill(HIST("jetPt"), jet.pt());
      registry.fill(HIST("jetEta"), jet.eta());
      registry.fill(HIST("jetPhi"), jet.phi());
      registry.fill(HIST("area"), jet.area());
      registry.fill(HIST("rho"), collision.rho());
      registry.fill(HIST("ptCorr"), ptCorr);
      // std::cout << nJet << std::endl;

      for (int i = 0; i < 14; i++) {
        double jetX = 0.05 * i + 0.025;
        double jetShapeFunction = ptDensity[i] / nJet;
        double jetShapeFunctionBg1 = ptDensityBg1[i] / nJet;
        double jetShapeFunctionBg2 = ptDensityBg2[i] / nJet;
        registry.fill(HIST("ptSum"), jetX, jetShapeFunction);
        registry.fill(HIST("ptSumBg1"), jetX, jetShapeFunctionBg1);
        registry.fill(HIST("ptSumBg2"), jetX, jetShapeFunctionBg2);
      }
    }
  }
  PROCESS_SWITCH(JetShapeTask, processCharged, "charged jets", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<JetShapeTask>(cfgc)}; }
