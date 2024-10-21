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

// Task to produce a table joinable to the jet tables for hf jet tagging
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>
/// \author Hanseo Park <hanseo.park@cern.ch>

#include <TF1.h>
#include <TH1.h>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "Common/Core/trackUtilities.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetTagging.h"
#include "PWGJE/Core/JetTaggingUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

template <typename JetTableData, typename JetTableMCD, typename JetTableMCP, typename JetTaggingTableData, typename JetTaggingTableMCD, typename JetTaggingTableMCP>
struct JetTaggerHFTask {

  Produces<JetTaggingTableData> taggingTableData;
  Produces<JetTaggingTableMCD> taggingTableMCD;
  Produces<JetTaggingTableMCP> taggingTableMCP;

  // configuration topological cut for track and sv
  Configurable<float> trackDcaXYMax{"trackDcaXYMax", 1, "minimum DCA xy acceptance for tracks [cm]"};
  Configurable<float> trackDcaZMax{"trackDcaZMax", 2, "minimum DCA z acceptance for tracks [cm]"};
  Configurable<float> prongsigmaLxyMax{"prongsigmaLxyMax", 100, "maximum sigma of decay length of prongs on xy plane"};
  Configurable<float> prongsigmaLxyzMax{"prongsigmaLxyzMax", 100, "maximum sigma of decay length of prongs on xyz plane"};
  Configurable<float> prongIPxyMin{"prongIPxyMin", 0.008, "maximum impact paramter of prongs on xy plane [cm]"};
  Configurable<float> prongIPxyMax{"prongIpxyMax", 1, "minimum impact parmeter of prongs on xy plane [cm]"};
  Configurable<float> prongChi2PCAMin{"prongChi2PCAMin", 4, "minimum Chi2 PCA of decay length of prongs"};
  Configurable<float> prongChi2PCAMax{"prongChi2PCAMax", 100, "maximum Chi2 PCA of decay length of prongs"};
  Configurable<float> svDispersionMax{"svDispersionMax", 1, "maximum dispersion of sv"};

  // jet flavour definition
  Configurable<float> maxDeltaR{"maxDeltaR", 0.25, "maximum distance of jet axis from flavour initiating parton"};
  Configurable<bool> removeGluonShower{"removeGluonShower", true, "find jet origin removed gluon spliting"}; // true:: remove gluon spliting
  Configurable<bool> searchUpToQuark{"searchUpToQuark", true, "Finding first mother in particles to quark"};

  // configuration about IP method
  Configurable<bool> useJetProb{"useJetProb", false, "fill table for track counting algorithm"};
  Configurable<bool> trackProbQA{"trackProbQA", false, "fill track probability histograms separately for geometric positive and negative tracks for QA"};
  Configurable<int> numCount{"numCount", 3, "number of track counting"};
  Configurable<int> resoFuncMatching{"resoFuncMatching", 0, "matching parameters of resolution function as MC samble (0: custom, 1: custom & inc, 2: MB, 3: MB & inc, 4: JJ, 5: JJ & inc)"};
  Configurable<std::vector<float>> paramsResoFuncData{"paramsResoFuncData", std::vector<float>{1306800, -0.1049, 0.861425, 13.7547, 0.977967, 8.96823, 0.151595, 6.94499, 0.0250301}, "parameters of gaus(0)+expo(3)+expo(5)+expo(7))"};
  Configurable<std::vector<float>> paramsResoFuncIncJetMC{"paramsResoFuncIncJetMC", std::vector<float>{1908803.027, -0.059, 0.895, 13.467, 1.005, 8.867, 0.098, 6.929, 0.011}, "parameters of gaus(0)+expo(3)+expo(5)+expo(7)))"};
  Configurable<std::vector<float>> paramsResoFuncCharmJetMC{"paramsResoFuncCharmJetMC", std::vector<float>{282119.753, -0.065, 0.893, 11.608, 0.945, 8.029, 0.131, 6.244, 0.027}, "parameters of gaus(0)+expo(3)+expo(5)+expo(7)))"};
  Configurable<std::vector<float>> paramsResoFuncBeautyJetMC{"paramsResoFuncBeautyJetMC", std::vector<float>{74901.583, -0.082, 0.874, 10.332, 0.941, 7.352, 0.097, 6.220, 0.022}, "parameters of gaus(0)+expo(3)+expo(5)+expo(7)))"};
  Configurable<std::vector<float>> paramsResoFuncLfJetMC{"paramsResoFuncLfJetMC", std::vector<float>{1539435.343, -0.061, 0.896, 13.272, 1.034, 5.884, 0.004, 7.843, 0.090}, "parameters of gaus(0)+expo(3)+expo(5)+expo(7)))"};
  Configurable<float> minSignImpXYSig{"minsIPs", -40.0, "minimum of signed impact parameter significance"};
  Configurable<int> minIPCount{"minIPCount", 2, "Select at least N signed impact parameter significance in jets"}; // default 2
  Configurable<float> tagPointForIP{"tagPointForIP", 2.5, "tagging working point for IP"};
  Configurable<float> tagPointForIPxyz{"tagPointForIPxyz", 2.5, "tagging working point for IP xyz"};
  // configuration about SV method
  Configurable<float> tagPointForSV{"tagPointForSV", 40, "tagging working point for SV"};
  Configurable<float> tagPointForSVxyz{"tagPointForSVxyz", 40, "tagging working point for SV xyz"};

  // axis spec
  ConfigurableAxis binTrackProbability{"binTrackProbability", {100, 0.f, 1.f}, ""};
  ConfigurableAxis binJetFlavour{"binJetFlavour", {6, -0.5, 5.5}, ""};

  using JetTagTracksData = soa::Join<aod::JetTracks, aod::JTrackExtras, aod::JTrackPIs>;
  using JetTagTracksMCD = soa::Join<aod::JetTracksMCD, aod::JTrackExtras, aod::JTrackPIs>;

  std::vector<float> vecParamsData;
  std::vector<float> vecParamsIncJetMC;
  std::vector<float> vecParamsCharmJetMC;
  std::vector<float> vecParamsBeautyJetMC;
  std::vector<float> vecParamsLfJetMC;
  std::vector<float> jetProb;
  bool useResoFuncFromIncJet = false;
  int maxOrder = -1;
  int resoFuncMatch = 0;
  std::unique_ptr<TF1> fSignImpXYSigData = nullptr;
  std::unique_ptr<TF1> fSignImpXYSigIncJetMC = nullptr;
  std::unique_ptr<TF1> fSignImpXYSigCharmJetMC = nullptr;
  std::unique_ptr<TF1> fSignImpXYSigBeautyJetMC = nullptr;
  std::unique_ptr<TF1> fSignImpXYSigLfJetMC = nullptr;

  template <typename T, typename U>
  void calculateJetProbability(int origin, T const& jet, U const& jtracks, std::vector<float>& jetProb, bool const& isMC = true)
  {
    jetProb.clear();
    jetProb.reserve(maxOrder);
    for (int order = 0; order < maxOrder; order++) {
      if (!isMC) {
        jetProb.push_back(jettaggingutilities::getJetProbability(fSignImpXYSigData, jet, jtracks, trackDcaXYMax, trackDcaZMax, order, tagPointForIP, minSignImpXYSig));
      } else {
        if (useResoFuncFromIncJet) {
          jetProb.push_back(jettaggingutilities::getJetProbability(fSignImpXYSigIncJetMC, jet, jtracks, trackDcaXYMax, trackDcaZMax, order, tagPointForIP, minSignImpXYSig));
        } else {
          if (origin == JetTaggingSpecies::charm) {
            jetProb.push_back(jettaggingutilities::getJetProbability(fSignImpXYSigCharmJetMC, jet, jtracks, trackDcaXYMax, trackDcaZMax, order, tagPointForIP, minSignImpXYSig));
          }
          if (origin == JetTaggingSpecies::beauty) {
            jetProb.push_back(jettaggingutilities::getJetProbability(fSignImpXYSigBeautyJetMC, jet, jtracks, trackDcaXYMax, trackDcaZMax, order, tagPointForIP, minSignImpXYSig));
          }
          if (origin == JetTaggingSpecies::lightflavour) {
            jetProb.push_back(jettaggingutilities::getJetProbability(fSignImpXYSigLfJetMC, jet, jtracks, trackDcaXYMax, trackDcaZMax, order, tagPointForIP, minSignImpXYSig));
          }
          if (origin != JetTaggingSpecies::charm && origin != JetTaggingSpecies::beauty && origin != JetTaggingSpecies::lightflavour) {
            jetProb.push_back(-1);
          }
        }
      }
    }
  }

  template <typename T, typename U>
  void evaluateTrackProbQA(int origin, T const& jet, U const& /*jtracks*/, bool const& isMC = true)
  {
    for (auto& jtrack : jet.template tracks_as<U>()) {
      if (!jettaggingutilities::trackAcceptanceWithDca(jtrack, trackDcaXYMax, trackDcaZMax))
        continue;
      auto geoSign = jettaggingutilities::getGeoSign(jet, jtrack);
      float probTrack = -1;
      if (!isMC) {
        probTrack = jettaggingutilities::getTrackProbability(fSignImpXYSigData, jtrack, minSignImpXYSig);
        if (geoSign > 0)
          registry.fill(HIST("h_pos_track_probability"), probTrack);
        else
          registry.fill(HIST("h_neg_track_probability"), probTrack);
      } else {
        if (useResoFuncFromIncJet) {
          probTrack = jettaggingutilities::getTrackProbability(fSignImpXYSigIncJetMC, jtrack, minSignImpXYSig);
        } else {
          if (origin == JetTaggingSpecies::charm) {
            probTrack = jettaggingutilities::getTrackProbability(fSignImpXYSigCharmJetMC, jtrack, minSignImpXYSig);
          }
          if (origin == JetTaggingSpecies::beauty) {
            probTrack = jettaggingutilities::getTrackProbability(fSignImpXYSigBeautyJetMC, jtrack, minSignImpXYSig);
          }
          if (origin == JetTaggingSpecies::lightflavour) {
            probTrack = jettaggingutilities::getTrackProbability(fSignImpXYSigLfJetMC, jtrack, minSignImpXYSig);
          }
        }
        if (geoSign > 0)
          registry.fill(HIST("h2_pos_track_probability_flavour"), probTrack, origin);
        else
          registry.fill(HIST("h2_neg_track_probability_flavour"), probTrack, origin);
      }
    }
  }

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};
  void init(InitContext const&)
  {
    maxOrder = numCount + 1; // 0: untagged, >1 : N ordering

    // Set up the resolution function
    resoFuncMatch = resoFuncMatching;
    switch (resoFuncMatch) {
      case 0:
        vecParamsData = (std::vector<float>)paramsResoFuncData;
        vecParamsCharmJetMC = (std::vector<float>)paramsResoFuncCharmJetMC;
        vecParamsBeautyJetMC = (std::vector<float>)paramsResoFuncBeautyJetMC;
        vecParamsLfJetMC = (std::vector<float>)paramsResoFuncLfJetMC;
        LOG(info) << "defined parameters of resolution function: custom";
        break;
      case 1:
        vecParamsData = (std::vector<float>)paramsResoFuncData;
        vecParamsIncJetMC = (std::vector<float>)paramsResoFuncIncJetMC;
        useResoFuncFromIncJet = true;
        LOG(info) << "defined parameters of resolution function: custom & use inclusive distribution";
        break;
      case 2: // TODO
        vecParamsData = (std::vector<float>)paramsResoFuncData;
        vecParamsCharmJetMC = {282119.753, -0.065, 0.893, 11.608, 0.945, 8.029, 0.131, 6.244, 0.027};
        vecParamsBeautyJetMC = {74901.583, -0.082, 0.874, 10.332, 0.941, 7.352, 0.097, 6.220, 0.022};
        vecParamsLfJetMC = {1539435.343, -0.061, 0.896, 13.272, 1.034, 5.884, 0.004, 7.843, 0.090};
        LOG(info) << "defined parameters of resolution function: PYTHIA8, MB, LHC23d1k";
        break;
      case 3: // TODO
        vecParamsData = (std::vector<float>)paramsResoFuncData;
        vecParamsIncJetMC = {1908803.027, -0.059, 0.895, 13.467, 1.005, 8.867, 0.098, 6.929, 0.011};
        LOG(info) << "defined parameters of resolution function: PYTHIA8, MB, LHC23d1k & use inclusive distribution";
        useResoFuncFromIncJet = true;
        break;
      case 4: // TODO
        vecParamsData = (std::vector<float>)paramsResoFuncData;
        vecParamsCharmJetMC = {743719.121, -0.960, -0.240, 13.765, 1.314, 10.761, 0.293, 8.538, 0.052};
        vecParamsBeautyJetMC = {88888.418, 0.256, 1.003, 10.185, 0.740, 8.216, 0.147, 7.228, 0.040};
        vecParamsLfJetMC = {414860.372, -1.000, 0.285, 14.561, 1.464, 11.693, 0.339, 9.183, 0.052};
        LOG(info) << "defined parameters of resolution function: PYTHIA8, JJ, weighted, LHC24g4";
        break;
      case 5: // TODO
        vecParamsData = (std::vector<float>)paramsResoFuncData;
        vecParamsIncJetMC = {2211391.862, 0.360, 1.028, 13.019, 0.650, 11.151, 0.215, 9.462, 0.044};
        LOG(info) << "defined parameters of resolution function: PYTHIA8, JJ, weighted, LHC24g4 & use inclusive distribution";
        useResoFuncFromIncJet = true;
        break;
      default:
        LOG(fatal) << "undefined parameters of resolution function. Fix it!";
        break;
    }

    fSignImpXYSigData = jettaggingutilities::setResolutionFunction(vecParamsData);
    fSignImpXYSigIncJetMC = jettaggingutilities::setResolutionFunction(vecParamsIncJetMC);
    fSignImpXYSigCharmJetMC = jettaggingutilities::setResolutionFunction(vecParamsCharmJetMC);
    fSignImpXYSigBeautyJetMC = jettaggingutilities::setResolutionFunction(vecParamsBeautyJetMC);
    fSignImpXYSigLfJetMC = jettaggingutilities::setResolutionFunction(vecParamsLfJetMC);

    // Use QA for effectivness of track probability
    if (trackProbQA) {
      AxisSpec trackProbabilityAxis = {binTrackProbability, "Track proability"};
      AxisSpec jetFlavourAxis = {binJetFlavour, "Jet flavour"};
      if (doprocessData || doprocessDataWithSV) {
        registry.add("h_pos_track_probability", "positive track probability", {HistType::kTH1F, {{trackProbabilityAxis}}});
        registry.add("h_neg_track_probability", "negative track probability", {HistType::kTH1F, {{trackProbabilityAxis}}});
      }
      if (doprocessMCD || doprocessMCDWithSV) {
        registry.add("h2_pos_track_probability_flavour", "positive track probability", {HistType::kTH2F, {{trackProbabilityAxis}, {jetFlavourAxis}}});
        registry.add("h2_neg_track_probability_flavour", "negative track probability", {HistType::kTH2F, {{trackProbabilityAxis}, {jetFlavourAxis}}});
      }
    }
  }

  void processDummy(aod::JetCollisions const&)
  {
  }
  PROCESS_SWITCH(JetTaggerHFTask, processDummy, "Dummy process", true);

  void processData(aod::JetCollision const& /*collision*/, JetTableData const& jets, JetTagTracksData const& jtracks)
  {
    for (auto& jet : jets) {
      bool flagtaggedjetIP = false;
      bool flagtaggedjetIPxyz = false;
      bool flagtaggedjetSV = false;
      bool flagtaggedjetSVxyz = false;
      if (useJetProb) {
        calculateJetProbability(0, jet, jtracks, jetProb, false);
        if (trackProbQA) {
          evaluateTrackProbQA(0, jet, jtracks, false);
        }
      }
      if (jettaggingutilities::isGreaterThanTaggingPoint(jet, jtracks, trackDcaXYMax, trackDcaZMax, tagPointForIP, minIPCount, false))
        flagtaggedjetIP = true;
      if (jettaggingutilities::isGreaterThanTaggingPoint(jet, jtracks, trackDcaXYMax, trackDcaZMax, tagPointForIP, minIPCount, true))
        flagtaggedjetIPxyz = true;

      taggingTableData(0, jetProb, flagtaggedjetIP, flagtaggedjetIPxyz, flagtaggedjetSV, flagtaggedjetSVxyz);
    }
  }
  PROCESS_SWITCH(JetTaggerHFTask, processData, "Fill tagging decision for data jets", false);

  void processDataWithSV(aod::JetCollision const& /*collision*/, soa::Join<JetTableData, aod::DataSecondaryVertex3ProngIndices> const& jets, JetTagTracksData const& jtracks, aod::DataSecondaryVertex3Prongs const& prongs)
  {
    for (auto& jet : jets) {
      bool flagtaggedjetIP = false;
      bool flagtaggedjetIPxyz = false;
      bool flagtaggedjetSV = false;
      bool flagtaggedjetSVxyz = false;
      if (useJetProb) {
        calculateJetProbability(0, jet, jtracks, jetProb, false);
        if (trackProbQA) {
          evaluateTrackProbQA(0, jet, jtracks, false);
        }
      }
      if (jettaggingutilities::isGreaterThanTaggingPoint(jet, jtracks, trackDcaXYMax, trackDcaZMax, tagPointForIP, minIPCount, false))
        flagtaggedjetIP = true;
      if (jettaggingutilities::isGreaterThanTaggingPoint(jet, jtracks, trackDcaXYMax, trackDcaZMax, tagPointForIP, minIPCount, true))
        flagtaggedjetIPxyz = true;
      flagtaggedjetSV = jettaggingutilities::isTaggedJetSV(jet, prongs, prongChi2PCAMin, prongChi2PCAMax, prongsigmaLxyMax, svDispersionMax, false, tagPointForSV);
      flagtaggedjetSVxyz = jettaggingutilities::isTaggedJetSV(jet, prongs, prongChi2PCAMin, prongChi2PCAMax, prongsigmaLxyzMax, svDispersionMax, true, tagPointForSV);
      taggingTableData(0, jetProb, flagtaggedjetIP, flagtaggedjetIPxyz, flagtaggedjetSV, flagtaggedjetSVxyz);
    }
  }
  PROCESS_SWITCH(JetTaggerHFTask, processDataWithSV, "Fill tagging decision for data jets", false);

  void processMCD(aod::JetCollision const& /*collision*/, JetTableMCD const& mcdjets, JetTagTracksMCD const& jtracks, aod::JetParticles const& particles)
  {
    for (auto& mcdjet : mcdjets) {
      bool flagtaggedjetIP = false;
      bool flagtaggedjetIPxyz = false;
      bool flagtaggedjetSV = false;
      bool flagtaggedjetSVxyz = false;
      typename JetTagTracksMCD::iterator hftrack;
      int origin = 0;
      if (removeGluonShower)
        origin = jettaggingutilities::mcdJetFromHFShower(mcdjet, jtracks, particles, maxDeltaR, searchUpToQuark);
      else
        origin = jettaggingutilities::jetTrackFromHFShower(mcdjet, jtracks, particles, hftrack, searchUpToQuark);
      if (useJetProb) {
        calculateJetProbability(origin, mcdjet, jtracks, jetProb);
        if (trackProbQA) {
          evaluateTrackProbQA(origin, mcdjet, jtracks);
        }
      }
      if (jettaggingutilities::isGreaterThanTaggingPoint(mcdjet, jtracks, trackDcaXYMax, trackDcaZMax, tagPointForIP, minIPCount, false))
        flagtaggedjetIP = true;
      if (jettaggingutilities::isGreaterThanTaggingPoint(mcdjet, jtracks, trackDcaXYMax, trackDcaZMax, tagPointForIP, minIPCount, true))
        flagtaggedjetIPxyz = true;
      taggingTableMCD(origin, jetProb, flagtaggedjetIP, flagtaggedjetIPxyz, flagtaggedjetSV, flagtaggedjetSVxyz);
    }
  }
  PROCESS_SWITCH(JetTaggerHFTask, processMCD, "Fill tagging decision for mcd jets", false);

  void processMCDWithSV(aod::JetCollision const& /*collision*/, soa::Join<JetTableMCD, aod::MCDSecondaryVertex3ProngIndices> const& mcdjets, JetTagTracksMCD const& jtracks, aod::MCDSecondaryVertex3Prongs const& prongs, aod::JetParticles const& particles)
  {
    for (auto& mcdjet : mcdjets) {
      bool flagtaggedjetIP = false;
      bool flagtaggedjetIPxyz = false;
      bool flagtaggedjetSV = false;
      bool flagtaggedjetSVxyz = false;
      typename JetTagTracksMCD::iterator hftrack;
      int origin = 0;
      if (removeGluonShower)
        origin = jettaggingutilities::mcdJetFromHFShower(mcdjet, jtracks, particles, maxDeltaR, searchUpToQuark);
      else
        origin = jettaggingutilities::jetTrackFromHFShower(mcdjet, jtracks, particles, hftrack, searchUpToQuark);
      if (useJetProb) {
        calculateJetProbability(origin, mcdjet, jtracks, jetProb);
        if (trackProbQA) {
          evaluateTrackProbQA(origin, mcdjet, jtracks);
        }
      }
      if (jettaggingutilities::isGreaterThanTaggingPoint(mcdjet, jtracks, trackDcaXYMax, trackDcaZMax, tagPointForIP, minIPCount, false))
        flagtaggedjetIP = true;
      if (jettaggingutilities::isGreaterThanTaggingPoint(mcdjet, jtracks, trackDcaXYMax, trackDcaZMax, tagPointForIP, minIPCount, true))
        flagtaggedjetIPxyz = true;
      flagtaggedjetSV = jettaggingutilities::isTaggedJetSV(mcdjet, prongs, prongChi2PCAMin, prongChi2PCAMax, prongsigmaLxyMax, prongIPxyMin, prongIPxyMax, svDispersionMax, false, tagPointForSV);
      flagtaggedjetSVxyz = jettaggingutilities::isTaggedJetSV(mcdjet, prongs, prongChi2PCAMin, prongChi2PCAMax, prongsigmaLxyzMax, prongIPxyMin, prongIPxyMax, svDispersionMax, true, tagPointForSV);
      taggingTableMCD(origin, jetProb, flagtaggedjetIP, flagtaggedjetIPxyz, flagtaggedjetSV, flagtaggedjetSVxyz);
    }
  }
  PROCESS_SWITCH(JetTaggerHFTask, processMCDWithSV, "Fill tagging decision for mcd jets with sv", false);

  void processMCP(aod::JetMcCollision const& /*collision*/, JetTableMCP const& mcpjets, aod::JetParticles const& particles)
  {
    for (auto& mcpjet : mcpjets) {
      bool flagtaggedjetIP = false;
      bool flagtaggedjetIPxyz = false;
      bool flagtaggedjetSV = false;
      bool flagtaggedjetSVxyz = false;
      typename aod::JetParticles::iterator hfparticle;
      int origin = 0;
      // TODO
      if (removeGluonShower) {
        if (jettaggingutilities::mcpJetFromHFShower(mcpjet, particles, maxDeltaR, searchUpToQuark))
          origin = jettaggingutilities::mcpJetFromHFShower(mcpjet, particles, maxDeltaR, searchUpToQuark);
        else
          origin = 0;
      } else {
        if (jettaggingutilities::jetParticleFromHFShower(mcpjet, particles, hfparticle, searchUpToQuark))
          origin = jettaggingutilities::jetParticleFromHFShower(mcpjet, particles, hfparticle, searchUpToQuark);
        else
          origin = 0;
      }
      jetProb.clear();
      jetProb.reserve(maxOrder);
      jetProb.push_back(-1);
      taggingTableMCP(origin, jetProb, flagtaggedjetIP, flagtaggedjetIPxyz, flagtaggedjetSV, flagtaggedjetSVxyz);
    }
  }
  PROCESS_SWITCH(JetTaggerHFTask, processMCP, "Fill tagging decision for mcp jets with sv", false);

  void processTraining(aod::JetCollision const& /*collision*/, JetTableMCD const& /*mcdjets*/, JetTagTracksMCD const& /*tracks*/)
  {
    // To create table for ML
  }
  PROCESS_SWITCH(JetTaggerHFTask, processTraining, "Fill tagging decision for mcd jets", false);
};

using JetTaggerChargedJets = JetTaggerHFTask<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>, soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>, soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents>, aod::ChargedJetTags, aod::ChargedMCDetectorLevelJetTags, aod::ChargedMCParticleLevelJetTags>;
using JetTaggerFullJets = JetTaggerHFTask<soa::Join<aod::FullJets, aod::FullJetConstituents>, soa::Join<aod::FullMCDetectorLevelJets, aod::FullMCDetectorLevelJetConstituents>, soa::Join<aod::FullMCParticleLevelJets, aod::FullMCParticleLevelJetConstituents>, aod::FullJetTags, aod::FullMCDetectorLevelJetTags, aod::FullMCParticleLevelJetTags>;
// using JetTaggerNeutralJets = JetTaggerHFTask<soa::Join<aod::NeutralJets, aod::NeutralJetConstituents>,soa::Join<aod::NeutralMCDetectorLevelJets, aod::NeutralMCDetectorLevelJetConstituents>, aod::NeutralJetTags, aod::NeutralMCDetectorLevelJetTags>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{

  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(
    adaptAnalysisTask<JetTaggerChargedJets>(cfgc,
                                            SetDefaultProcesses{}, TaskName{"jet-taggerhf-charged"}));

  tasks.emplace_back(
    adaptAnalysisTask<JetTaggerFullJets>(cfgc,
                                         SetDefaultProcesses{}, TaskName{"jet-taggerhf-full"}));
  /*
    tasks.emplace_back(
      adaptAnalysisTask<JetTaggerNeutralJets>(cfgc,
                                                  SetDefaultProcesses{}, TaskName{"jet-taggerhf-neutral"}));
  */
  return WorkflowSpec{tasks};
}
