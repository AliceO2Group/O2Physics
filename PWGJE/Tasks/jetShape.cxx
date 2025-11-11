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

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/DataModel/Jet.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include <TPDGCode.h>

#include <cmath>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetShapeTask {

  Configurable<int> nBinsNSigma{"nBinsNSigma", 101, "Number of nsigma bins"};
  Configurable<float> nSigmaMin{"nSigmaMin", -10.1f, "Min value of nsigma"};
  Configurable<float> nSigmaMax{"nSigmaMax", 10.1f, "Max value of nsigma"};
  Configurable<int> nBinsPForDedx{"nBinsPForDedx", 700, "Number of p bins"};
  Configurable<int> nBinsPForBeta{"nBinsPForBeta", 500, "Number of pT bins"};
  Configurable<int> nBinsTpcDedx{"nBinsTpcDedx", 500, "Number of DEdx bins"};
  Configurable<int> nBinsTofBeta{"nBinsTofBeta", 350, "Number of Beta bins"};
  Configurable<float> pMax{"pMax", 7.0f, "Max value of p"};
  Configurable<float> ptMax{"ptMax", 5.0f, "Max value of pT"};
  Configurable<float> jetPtMinForCut{"jetPtMinForCut", 0.0f, "Minimum value of jet pT cut"};
  Configurable<float> jetPtMaxForCut{"jetPtMaxForCut", 200.0f, "Maximum value of the jet pT cut"};
  Configurable<int> nBinsP{"nBinsP", 70, "Number of p bins"};
  Configurable<int> nBinsPt{"nBinsPt", 50, "Number of pT bins"};
  Configurable<int> nBinsJetPt{"nBinsJetPt", 10, "Number of jet pT bins"};
  Configurable<int> nBinsDistance{"nBinsDistance", 7, "Number of distance bins"};
  Configurable<float> distanceMax{"distanceMax", 0.7f, "Max value of distance"};
  Configurable<float> nSigmaTofCut{"nSigmaTofCut", 2.0f, "Number of sigma cut for TOF PID"};
  Configurable<float> tpcNSigmaPrMin{"tpcNSigmaPrMin", -3.5f, "Min value of tpcNsigmaProton"};
  Configurable<float> tpcNSigmaPrMax{"tpcNSigmaPrMax", 0.5f, "Max value of tpcNsigmaProton"};
  Configurable<float> tpcNSigmaPiMin{"tpcNSigmaPiMin", -0.5f, "Min value of tpcNsigmaPion"};
  Configurable<float> tpcNSigmaPiMax{"tpcNSigmaPiMax", 3.5f, "Max value of tpcNsigmaPion"};

  HistogramRegistry registry{"registry",
                             {{"tpcTofPi", "tpcTofPi", {HistType::kTHnSparseD, {{35, 0, pMax}, {nBinsNSigma, nSigmaMin, nSigmaMax}, {nBinsDistance, 0, distanceMax}, {nBinsJetPt, jetPtMinForCut, jetPtMaxForCut}}}},
                              {"tpcTofPr", "tpcTofPr", {HistType::kTHnSparseD, {{35, 0, pMax}, {nBinsNSigma, nSigmaMin, nSigmaMax}, {nBinsDistance, 0, distanceMax}, {nBinsJetPt, jetPtMinForCut, jetPtMaxForCut}}}},
                              {"tpcTofPiOutOfJet", "tpcTofPiOutOfJet", {HistType::kTHnSparseD, {{35, 0, pMax}, {nBinsNSigma, nSigmaMin, nSigmaMax}, {nBinsJetPt, jetPtMinForCut, jetPtMaxForCut}}}},
                              {"tpcTofPrOutOfJet", "tpcTofPrOutOfJet", {HistType::kTHnSparseD, {{35, 0, pMax}, {nBinsNSigma, nSigmaMin, nSigmaMax}, {nBinsJetPt, jetPtMinForCut, jetPtMaxForCut}}}},
                              {"tpcPi", "tpcPi", {HistType::kTH2F, {{nBinsP, 0, pMax}, {nBinsNSigma, nSigmaMin, nSigmaMax}}}},
                              {"tofPi", "tofPi", {HistType::kTH2F, {{nBinsPt, 0, ptMax}, {nBinsNSigma, nSigmaMin, nSigmaMax}}}},
                              {"tpcPr", "tpcPr", {HistType::kTH2F, {{nBinsP, 0, pMax}, {nBinsNSigma, nSigmaMin, nSigmaMax}}}},
                              {"tofPr", "tofPr", {HistType::kTH2F, {{nBinsPt, 0, ptMax}, {nBinsNSigma, nSigmaMin, nSigmaMax}}}},
                              {"tpcDedx", "tpcDedx", {HistType::kTHnSparseD, {{nBinsPForDedx, 0, pMax}, {nBinsTpcDedx, 0, 1000}, {nBinsDistance, 0, distanceMax}}}},
                              {"tpcDedxOutOfJet", "tpcDedxOutOfJet", {HistType::kTH2F, {{nBinsPForDedx, 0, pMax}, {nBinsTpcDedx, 0, 1000}}}},
                              {"tofBeta", "tofBeta", {HistType::kTH2F, {{nBinsPForBeta, 0, pMax}, {nBinsTofBeta, 0.4, 1.1}}}},
                              {"pVsPtForPr", "pVsPtForPr", {HistType::kTHnSparseD, {{nBinsP, 0, pMax}, {nBinsPt, 0, ptMax}, {nBinsDistance, 0, distanceMax}, {nBinsJetPt, jetPtMinForCut, jetPtMaxForCut}}}},
                              {"pVsPtForPi", "pVsPtPi", {HistType::kTHnSparseD, {{nBinsP, 0, pMax}, {nBinsPt, 0, ptMax}, {nBinsDistance, 0, distanceMax}, {nBinsJetPt, jetPtMinForCut, jetPtMaxForCut}}}},
                              {"pVsPtForPrOutOfJet", "pVsPtForPrOutOfJet", {HistType::kTHnSparseD, {{nBinsP, 0, pMax}, {nBinsPt, 0, ptMax}, {nBinsJetPt, jetPtMinForCut, jetPtMaxForCut}}}},
                              {"pVsPtForPiOutOfJet", "pVsPtPionOutOfJet", {HistType::kTHnSparseD, {{nBinsP, 0, pMax}, {nBinsPt, 0, ptMax}, {nBinsJetPt, jetPtMinForCut, jetPtMaxForCut}}}},
                              {"tofMass", "tofMass", {HistType::kTH1F, {{300, 0, 3}}}},
                              {"trackPhi", "trackPhi", {HistType::kTH1F, {{80, -1, 7}}}},
                              {"trackEta", "trackEta", {HistType::kTH1F, {{100, -1, 1}}}},
                              {"trackTpcNClsCrossedRows", "trackTpcNClsCrossedRows", {HistType::kTH1F, {{50, 0, 200}}}},
                              {"trackDcaXY", "trackDcaXY", {HistType::kTH1F, {{40, -10, 10}}}},
                              {"trackItsChi2NCl", "trackItsChi2NCl", {HistType::kTH1F, {{60, 0, 30}}}},
                              {"trackTpcChi2NCl", "trackTpcChi2NCl", {HistType::kTH1F, {{100, 0, 50}}}},
                              {"trackTpcNClsFound", "trackTpcNClsFound", {HistType::kTH1F, {{100, 0, 200}}}},
                              {"trackItsNCls", "trackItsNCls", {HistType::kTH1F, {{10, 0, 10}}}},
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
                              {"eventCounter", "eventCounter", {HistType::kTH1F, {{1, 0, +1, ""}}}},
                              {"ptVsCentrality", "ptvscentrality", {HistType::kTH2F, {{100, 0, 100}, {300, 0, 300}}}},
                              {"ptResolution", "ptResolution", {HistType::kTH2F, {{nBinsPt, 0, ptMax}, {100, -1.0, +1.0}}}},
                              {"ptHistogramPion", "ptHistogramPion", {HistType::kTH1F, {{nBinsPt, 0, ptMax}}}},
                              {"ptHistogramKaon", "ptHistogramKaon", {HistType::kTH1F, {{nBinsPt, 0, ptMax}}}},
                              {"ptHistogramProton", "ptHistogramProton", {HistType::kTH1F, {{nBinsPt, 0, ptMax}}}},
                              {"ptHistogramPionTof", "ptHistogramPionTof", {HistType::kTH1F, {{nBinsPt, 0, ptMax}}}},
                              {"ptHistogramKaonTof", "ptHistogramKaonTof", {HistType::kTH1F, {{nBinsPt, 0, ptMax}}}},
                              {"ptHistogramProtonTof", "ptHistogramProtonTof", {HistType::kTH1F, {{nBinsPt, 0, ptMax}}}},
                              {"ptGeneratedPion", "ptGeneratedPion", {HistType::kTH1F, {{nBinsPt, 0, ptMax}}}},
                              {"ptGeneratedKaon", "ptGeneratedKaon", {HistType::kTH1F, {{nBinsPt, 0, ptMax}}}},
                              {"ptGeneratedProton", "ptGeneratedProton", {HistType::kTH1F, {{nBinsPt, 0, ptMax}}}}}};

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
  Configurable<float> mcRapidityMax{"mcRapidityMax", 0.5f, "maximum mctrack y"};

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
    static constexpr double JetAreaFractionMinValue = -98.0;
    if (jetAreaFractionMin > JetAreaFractionMinValue) {
      if (jet.area() < jetAreaFractionMin * o2::constants::math::PI * (jet.r() / 100.0) * (jet.r() / 100.0)) {
        return false;
      }
      if (jet.area() < o2::constants::math::PIHalf * (jet.r() / 100.0) * (jet.r() / 100.0)) {
        return false;
      }
    }
    static constexpr double LeadingConstituentPtMinValue = 5.0;
    static constexpr double LeadingConstituentPtMaxValue = 9998.0;
    bool checkConstituentPt = true;
    bool checkConstituentMinPt = (leadingConstituentPtMin > LeadingConstituentPtMinValue);
    bool checkConstituentMaxPt = (leadingConstituentPtMax < LeadingConstituentPtMaxValue);
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
        registry.fill(HIST("ptVsCentrality"), collision.centFT0M(), track.pt());

        // calculate compornents of jetshapefunction rho(r)
        std::vector<float> trackPtSum(distanceCategory->size() - 1, 0.f);
        std::vector<float> trackPtSumBg1(distanceCategory->size() - 1, 0.f);
        std::vector<float> trackPtSumBg2(distanceCategory->size() - 1, 0.f);

        float phiBg1 = jet.phi() + (o2::constants::math::PIHalf);
        float phiBg2 = jet.phi() - (o2::constants::math::PIHalf);

        float preDeltaPhiBg1 = track.phi() - phiBg1;
        float preDeltaPhiBg2 = track.phi() - phiBg2;

        float deltaPhiBg1 = RecoDecay::constrainAngle(preDeltaPhiBg1, -o2::constants::math::PI);
        float deltaPhiBg2 = RecoDecay::constrainAngle(preDeltaPhiBg2, -o2::constants::math::PI);

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
  PROCESS_SWITCH(JetShapeTask, processJetShape, "JetShape", false);

  void processProductionRatio(soa::Filtered<aod::JetCollisions>::iterator const& collision, soa::Join<aod::JetTracks, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::TracksExtra, aod::TracksDCA, aod::pidTOFbeta, aod::pidTOFmass> const& tracks, soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets)
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
        registry.fill(HIST("trackTpcNClsCrossedRows"), track.tpcNClsCrossedRows());
        registry.fill(HIST("trackDcaXY"), track.dcaXY());
        registry.fill(HIST("trackItsChi2NCl"), track.itsChi2NCl());
        registry.fill(HIST("trackTpcChi2NCl"), track.tpcChi2NCl());
        registry.fill(HIST("trackTpcNClsFound"), track.tpcNClsFound());
        registry.fill(HIST("trackItsNCls"), track.itsNCls());
        registry.fill(HIST("trackEta"), track.eta());
        registry.fill(HIST("trackPhi"), track.phi());

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
        registry.fill(HIST("tofMass"), track.mass());
        registry.fill(HIST("tpcPi"), track.p(), track.tpcNSigmaPi());
        registry.fill(HIST("tofPi"), track.pt(), track.tofNSigmaPi());
        registry.fill(HIST("tpcPr"), track.p(), track.tpcNSigmaPr());
        registry.fill(HIST("tofPr"), track.pt(), track.tofNSigmaPr());

        // for calculate distance
        float preDeltaPhi1 = track.phi() - jet.phi();
        float deltaPhi1 = RecoDecay::constrainAngle(preDeltaPhi1);
        float deltaEta = track.eta() - jet.eta();

        // calculate distance from jet axis
        float distance = std::sqrt(deltaEta * deltaEta + deltaPhi1 * deltaPhi1);

        // Define perpendicular cone axes in phi
        float phiBg1 = jet.phi() + (o2::constants::math::PIHalf);
        float phiBg2 = jet.phi() - (o2::constants::math::PIHalf);

        // Calculate delta phi for background cones
        float preDeltaPhiBg1 = track.phi() - phiBg1;
        float preDeltaPhiBg2 = track.phi() - phiBg2;
        float deltaPhiBg1 = RecoDecay::constrainAngle(preDeltaPhiBg1, -o2::constants::math::PI);
        float deltaPhiBg2 = RecoDecay::constrainAngle(preDeltaPhiBg2, -o2::constants::math::PI);

        // Calculate distance to background cone axes
        float distanceBg1 = std::sqrt(deltaEta * deltaEta + deltaPhiBg1 * deltaPhiBg1);
        float distanceBg2 = std::sqrt(deltaEta * deltaEta + deltaPhiBg2 * deltaPhiBg2);

        // Fill histogram if track is inside one of the perpendicular cones
        if (distanceBg1 < jetR || distanceBg2 < jetR) {
          registry.fill(HIST("tpcDedxOutOfJet"), track.p(), track.tpcSignal());

          if (std::abs(track.tofNSigmaPi()) < nSigmaTofCut) {
            registry.fill(HIST("tpcTofPiOutOfJet"), track.p(), track.tpcNSigmaPi(), jet.pt());
            if (track.tpcNSigmaPi() > tpcNSigmaPiMin && track.tpcNSigmaPi() < tpcNSigmaPiMax) {
              registry.fill(HIST("pVsPtForPiOutOfJet"), track.p(), track.pt(), jet.pt());
            }
          }
          if (std::abs(track.tofNSigmaPr()) < nSigmaTofCut) {
            registry.fill(HIST("tpcTofPrOutOfJet"), track.p(), track.tpcNSigmaPr(), jet.pt());
            if (track.tpcNSigmaPr() > tpcNSigmaPrMin && track.tpcNSigmaPr() < tpcNSigmaPrMax) {
              registry.fill(HIST("pVsPtForPrOutOfJet"), track.p(), track.pt(), jet.pt());
            }
          }
        }

        registry.fill(HIST("distanceVsTrackpt"), distance, track.pt());
        registry.fill(HIST("tpcDedx"), track.p(), track.tpcSignal(), distance);
        registry.fill(HIST("tofBeta"), track.p(), track.beta());

        if (std::abs(track.tofNSigmaPr()) < nSigmaTofCut) {
          registry.fill(HIST("tpcTofPr"), track.p(), track.tpcNSigmaPr(), distance, jet.pt());
          if (track.tpcNSigmaPr() > tpcNSigmaPrMin && track.tpcNSigmaPr() < tpcNSigmaPrMax) {
            registry.fill(HIST("pVsPtForPr"), track.p(), track.pt(), distance, jet.pt());
          }
        }

        if (std::abs(track.tofNSigmaPi()) < nSigmaTofCut) {
          registry.fill(HIST("tpcTofPi"), track.p(), track.tpcNSigmaPi(), distance, jet.pt());
          if (track.tpcNSigmaPi() > tpcNSigmaPiMin && track.tpcNSigmaPi() < tpcNSigmaPiMax) {
            registry.fill(HIST("pVsPtForPi"), track.p(), track.pt(), distance, jet.pt());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(JetShapeTask, processProductionRatio, "production ratio", false);

  void processReco(soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels> const& tracks, aod::McParticles const&)
  {
    registry.fill(HIST("eventCounter"), 0.5);

    for (const auto& track : tracks) {
      if (track.has_mcParticle()) {
        auto mcParticle = track.mcParticle();
        registry.fill(HIST("ptResolution"), track.pt(), track.pt() - mcParticle.pt());

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

        if (mcParticle.isPhysicalPrimary() && std::fabs(mcParticle.y()) < mcRapidityMax) { // do this in the context of the track ! (context matters!!!)
          if (std::abs(mcParticle.pdgCode()) == PDG_t::kPiPlus)
            registry.fill(HIST("ptHistogramPion"), mcParticle.pt());
          if (std::abs(mcParticle.pdgCode()) == PDG_t::kKPlus)
            registry.fill(HIST("ptHistogramKaon"), mcParticle.pt());
          if (std::abs(mcParticle.pdgCode()) == PDG_t::kProton)
            registry.fill(HIST("ptHistogramProton"), mcParticle.pt());
        }

        if (track.hasTOF()) {
          if (mcParticle.isPhysicalPrimary() && std::fabs(mcParticle.y()) < mcRapidityMax) {
            if (std::abs(mcParticle.pdgCode()) == PDG_t::kPiPlus)
              registry.fill(HIST("ptHistogramPionTof"), mcParticle.pt());
            if (std::abs(mcParticle.pdgCode()) == PDG_t::kKPlus)
              registry.fill(HIST("ptHistogramKaonTof"), mcParticle.pt());
            if (std::abs(mcParticle.pdgCode()) == PDG_t::kProton)
              registry.fill(HIST("ptHistogramProtonTof"), mcParticle.pt());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(JetShapeTask, processReco, "process reconstructed information", true);

  void processSim(aod::McParticles const& mcParticles)
  {
    for (const auto& mcParticle : mcParticles) {
      if (mcParticle.isPhysicalPrimary() && std::fabs(mcParticle.y()) < mcRapidityMax) {
        if (std::abs(mcParticle.pdgCode()) == PDG_t::kPiPlus)
          registry.fill(HIST("ptGeneratedPion"), mcParticle.pt());
        if (std::abs(mcParticle.pdgCode()) == PDG_t::kKPlus)
          registry.fill(HIST("ptGeneratedKaon"), mcParticle.pt());
        if (std::abs(mcParticle.pdgCode()) == PDG_t::kProton)
          registry.fill(HIST("ptGeneratedProton"), mcParticle.pt());
      }
    }
  }
  PROCESS_SWITCH(JetShapeTask, processSim, "process pure simulation information", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<JetShapeTask>(cfgc)}; }
