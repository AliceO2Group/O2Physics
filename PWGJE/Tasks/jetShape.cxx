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

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetSubtraction.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include <CommonConstants/MathConstants.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetShapeTask {
  HistogramRegistry registry{"registry",
                             {{"tpcTofPi", "tpcTofPi", {HistType::kTHnSparseD, {{101, -10.1f, 10.1f}, {20, -10, 10}, {25, 0, 5}, {14, 0, 7}}}},
                              {"tpcPi", "tpcPi", {HistType::kTH2F, {{100, 0, 5}, {401, -10.025f, 10.025f}}}},
                              {"tofPi", "tofPi", {HistType::kTH2F, {{100, 0, 5}, {401, -10.025f, 10.025f}}}},
                              {"tpcTofPr", "tpcTofPr", {HistType::kTHnSparseD, {{101, -10.1f, 10.1f}, {20, -10, 10}, {25, 0, 5}, {14, 0, 7}}}},
                              {"tpcPr", "tpcPr", {HistType::kTH2F, {{100, 0, 5}, {401, -10.025f, 10.025f}}}},
                              {"tofPr", "tofPr", {HistType::kTH2F, {{100, 0, 5}, {401, -10.025f, 10.025f}}}},
                              {"tpcDedx", "tpcDedx", {HistType::kTH2F, {{500, 0, 5}, {1000, 0, 1000}}}},
                              {"tofBeta", "tofBeta", {HistType::kTH2F, {{500, 0, 5}, {450, 0.2, 1.1}}}},
                              {"tofMass", "tofMass", {HistType::kTH1F, {{3000, 0, 3}}}},
                              {"jetPt", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"jetEta", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
                              {"jetPhi", "jet #phi;#phi_{jet};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
                              {"area", "area", {HistType::kTH1F, {{200, 0, 4}}}},
                              {"rho", "rho", {HistType::kTH1F, {{300, 0, 300}}}},
                              {"ptCorr", "Corrected jet pT; p_{T}^{corr} (GeV/c); Counts", {HistType::kTH1F, {{200, 0, 200}}}},
                              {"ptCorrVsDistance", "ptcorr_vs_distance", {HistType::kTH2F, {{70, 0, 0.7}, {100, 0, 100}}}},
                              {"distanceVsTrackpt", "trackpt_vs_distance", {HistType::kTH2F, {{70, 0, 0.7}, {100, 0, 100}}}},
                              {"ptSum", "ptSum", {HistType::kTH2F, {{14, 0, 0.7}, {300, 0, 300}}}},
                              {"ptSumBg1", "ptSumBg1", {HistType::kTH2F, {{14, 0, 0.7}, {300, 0, 300}}}},
                              {"ptSumBg2", "ptSumBg2", {HistType::kTH2F, {{14, 0, 0.7}, {300, 0, 300}}}},
                              {"event/vertexz", ";Vtx_{z} (cm);Entries", {HistType::kTH1F, {{100, -20, 20}}}},
                              {"ptVsCentrality", "ptvscentrality", {HistType::kTH2F, {{100, 0, 100}, {300, 0, 300}}}}}};

  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};

  Configurable<float> jetPtMin{"jetPtMin", 5.0, "minimum jet pT cut"};
  Configurable<float> jetR{"jetR", 0.4, "jet resolution parameter"};

  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

  Configurable<float> jetAreaFractionMin{"jetAreaFractionMin", -99.0, "used to make a cut on the jet areas"};
  Configurable<float> leadingConstituentPtMin{"leadingConstituentPtMin", 5.0, "minimum pT selection on jet constituent"};
  Configurable<float> leadingConstituentPtMax{"leadingConstituentPtMax", 9999.0, "maximum pT selection on jet constituent"};

  // for jet shape
  Configurable<std::vector<float>> distanceCategory{"distanceCategory", {0.00f, 0.05f, 0.10f, 0.15f, 0.20f, 0.25f, 0.30f, 0.35f, 0.40f, 0.45f, 0.50f, 0.55f, 0.60f, 0.65f, 0.70f}, "distance of category"};

  // for ppi production
  Configurable<float> etaTrUp{"etaTrUp", 0.7f, "maximum track eta"};
  Configurable<float> dcaxyMax{"dcaxyMax", 2.0f, "mximum DCA xy"};
  Configurable<float> chi2ItsMax{"chi2ItsMax", 15.0f, "its chi2 cut"};
  Configurable<float> chi2TpcMax{"chi2TpcMax", 4.0f, "tpc chi2 cut"};
  Configurable<float> nclItsMin{"nclItsMin", 2.0f, "its # of cluster cut"};
  Configurable<float> nclTpcMin{"nclTpcMin", 100.0f, "tpc # if cluster cut"};
  Configurable<float> nclcrossTpcMin{"nclcrossTpcMin", 70.0f, "tpc # of crossedRows cut"};

  Configurable<std::string> triggerMasks{"triggerMasks", "", "possible JE Trigger masks: fJetChLowPt,fJetChHighPt,fTrackLowPt,fTrackHighPt,fJetD0ChLowPt,fJetD0ChHighPt,fJetLcChLowPt,fJetLcChHighPt,fEMCALReadout,fJetFullHighPt,fJetFullLowPt,fJetNeutralHighPt,fJetNeutralLowPt,fGammaVeryHighPtEMCAL,fGammaVeryHighPtDCAL,fGammaHighPtEMCAL,fGammaHighPtDCAL,fGammaLowPtEMCAL,fGammaLowPtDCAL,fGammaVeryLowPtEMCAL,fGammaVeryLowPtDCAL"};

  std::vector<int> eventSelectionBits;
  int trackSelection = -1;
  std::vector<int> triggerMaskBits;

  void init(o2::framework::InitContext&)
  {
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
    triggerMaskBits = jetderiveddatautilities::initialiseTriggerMaskBits(triggerMasks);
  }

  template <typename T, typename U>
  bool isAcceptedJet(U const& jet)
  {
    static constexpr double kJetAreaFractionMinValue = -98.0;
    if (jetAreaFractionMin > kJetAreaFractionMinValue) {
      if (jet.area() < jetAreaFractionMin * o2::constants::math::PI * (jet.r() / 100.0) * (jet.r() / 100.0)) {
        return false;
      }
      if (jet.area() < o2::constants::math::PIHalf * (jet.r() / 100.0) * (jet.r() / 100.0)) {
        return false;
      }
    }
    static constexpr double kLeadingConstituentPtMinValue = 5.0;
    static constexpr double kLeadingConstituentPtMaxValue = 9998.0;
    bool checkConstituentPt = true;
    bool checkConstituentMinPt = (leadingConstituentPtMin > kLeadingConstituentPtMinValue);
    bool checkConstituentMaxPt = (leadingConstituentPtMax < kLeadingConstituentPtMaxValue);
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

  void processJetShape(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos>>::iterator const& collision, aod::JetTracks const& tracks, soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets)
  {

    std::vector<float> ptDensity(distanceCategory->size() - 1, 0.f);
    std::vector<float> ptDensityBg1(distanceCategory->size() - 1, 0.f);
    std::vector<float> ptDensityBg2(distanceCategory->size() - 1, 0.f);

    for (auto const& jet : jets) {
      if (!isAcceptedJet<aod::JetTracks>(jet)) {
        continue;
      }

      // Get underlying event subtracted jet.pt() as ptCorr
      float ptCorr = jet.pt() - collision.rho() * jet.area();

      for (const auto& track : tracks) {
        float preDeltaPhi1 = track.phi() - jet.phi();
        float deltaPhi1 = RecoDecay::constrainAngle(preDeltaPhi1);
        float deltaEta = track.eta() - jet.eta();

        // calculate distance from jet axis
        float distance = std::sqrt(deltaEta * deltaEta + deltaPhi1 * deltaPhi1);

        registry.fill(HIST("ptCorrVsDistance"), distance, ptCorr);
        registry.fill(HIST("ptVsCentrality"), collision.centrality(), track.pt());

        // calculate compornents of jetshapefunction rho(r)
        std::vector<float> trackPtSum(distanceCategory->size() - 1, 0.f);
        std::vector<float> trackPtSumBg1(distanceCategory->size() - 1, 0.f);
        std::vector<float> trackPtSumBg2(distanceCategory->size() - 1, 0.f);

        float phiBg1 = jet.phi() + (o2::constants::math::PIHalf);
        float phiBg2 = jet.phi() - (o2::constants::math::PIHalf);

        float preDeltaPhiBg1 = track.phi() - phiBg1;
        float preDeltaPhiBg2 = track.phi() - phiBg2;

        float deltaPhiBg1 = RecoDecay::constrainAngle(preDeltaPhiBg1);
        float deltaPhiBg2 = RecoDecay::constrainAngle(preDeltaPhiBg2);

        float distanceBg1 = std::sqrt(deltaEta * deltaEta + deltaPhiBg1 * deltaPhiBg1);
        float distanceBg2 = std::sqrt(deltaEta * deltaEta + deltaPhiBg2 * deltaPhiBg2);

        for (size_t i = 0; i < distanceCategory->size() - 1; i++) {
          if (distanceCategory->at(i) <= distance && distance < distanceCategory->at(i + 1))
            trackPtSum[i] += track.pt();
          if (distanceCategory->at(i) <= distanceBg1 && distanceBg1 < distanceCategory->at(i + 1))
            trackPtSumBg1[i] += track.pt();
          if (distanceCategory->at(i) <= distanceBg2 && distanceBg2 < distanceCategory->at(i + 1))
            trackPtSumBg2[i] += track.pt();
        }

        for (size_t i = 0; i < distanceCategory->size() - 1; i++) {
          ptDensity[i] += trackPtSum[i] / ((distanceCategory->at(i + 1) - distanceCategory->at(i)) * ptCorr);
          ptDensityBg1[i] += trackPtSumBg1[i] / ((distanceCategory->at(i + 1) - distanceCategory->at(i)) * ptCorr);
          ptDensityBg2[i] += trackPtSumBg2[i] / ((distanceCategory->at(i + 1) - distanceCategory->at(i)) * ptCorr);
        }
      }

      registry.fill(HIST("jetPt"), jet.pt());
      registry.fill(HIST("jetEta"), jet.eta());
      registry.fill(HIST("jetPhi"), jet.phi());
      registry.fill(HIST("area"), jet.area());
      registry.fill(HIST("rho"), collision.rho());
      registry.fill(HIST("ptCorr"), ptCorr);

      for (size_t i = 0; i < distanceCategory->size() - 1; i++) {
        double jetX = (distanceCategory->at(i + 1) - distanceCategory->at(i)) * i + (distanceCategory->at(i + 1) - distanceCategory->at(i)) / 2;
        double jetShapeFunction = ptDensity[i];
        double jetShapeFunctionBg1 = ptDensityBg1[i];
        double jetShapeFunctionBg2 = ptDensityBg2[i];
        registry.fill(HIST("ptSum"), jetX, jetShapeFunction);
        registry.fill(HIST("ptSumBg1"), jetX, jetShapeFunctionBg1);
        registry.fill(HIST("ptSumBg2"), jetX, jetShapeFunctionBg2);
      }
    }
  }
  PROCESS_SWITCH(JetShapeTask, processJetShape, "JetShape", true);

  void processProductionRatio(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos>>::iterator const& collision, soa::Join<aod::JetTracks, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::TracksExtra, aod::TracksDCA, aod::pidTOFbeta, aod::pidTOFmass> const& tracks, soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }

    registry.fill(HIST("event/vertexz"), collision.posZ());

    for (auto const& jet : jets) {
      if (!isAcceptedJet<aod::JetTracks>(jet)) {
        continue;
      }

      // tracks conditions
      for (const auto& track : tracks) {
        if (std::abs(track.eta()) > etaTrUp)
          continue;
        if (track.tpcNClsCrossedRows() < nclcrossTpcMin)
          continue;
        if (std::abs(track.dcaXY()) > dcaxyMax)
          continue;
        if (track.itsChi2NCl() > chi2ItsMax)
          continue;
        if (track.tpcChi2NCl() > chi2TpcMax)
          continue;
        if (track.tpcNClsFound() < nclTpcMin)
          continue;
        if (track.itsNCls() < nclItsMin)
          continue;

        // PID check
        registry.fill(HIST("tpcDedx"), track.pt(), track.tpcSignal());
        registry.fill(HIST("tofBeta"), track.pt(), track.beta());
        registry.fill(HIST("tofMass"), track.mass());

        // for calculate purity
        registry.fill(HIST("tpcPi"), track.pt(), track.tpcNSigmaPi());
        registry.fill(HIST("tofPi"), track.pt(), track.tofNSigmaPi());
        registry.fill(HIST("tpcPr"), track.pt(), track.tpcNSigmaPr());
        registry.fill(HIST("tofPr"), track.pt(), track.tofNSigmaPr());

        // for calculate distance
        float preDeltaPhi1 = track.phi() - jet.phi();
        float deltaPhi1 = RecoDecay::constrainAngle(preDeltaPhi1);
        float deltaEta = track.eta() - jet.eta();

        // calculate distance from jet axis
        float distance = std::sqrt(deltaEta * deltaEta + deltaPhi1 * deltaPhi1);

        registry.fill(HIST("distanceVsTrackpt"), distance, track.pt());
        registry.fill(HIST("tpcTofPi"), track.tpcNSigmaPi(), track.tofNSigmaPi(), track.pt(), distance);
        registry.fill(HIST("tpcTofPr"), track.tpcNSigmaPr(), track.tofNSigmaPr(), track.pt(), distance);
      }
    }
  }
  PROCESS_SWITCH(JetShapeTask, processProductionRatio, "production ratio", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<JetShapeTask>(cfgc)}; }
