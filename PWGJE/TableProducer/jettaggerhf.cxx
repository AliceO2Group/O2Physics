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

template <typename JetTableData, typename JetTableMCD, typename JetTaggingTableData, typename JetTaggingTableMCD>
struct JetTaggerHFTask {

  Produces<JetTaggingTableData> taggingTableData;
  Produces<JetTaggingTableMCD> taggingTableMCD;

  Configurable<float> maxDeltaR{"maxDeltaR", 0.25, "maximum distance of jet axis from flavour initiating parton"};
  Configurable<bool> removeGluonShower{"removeGluonShower", true, "find jet origin removed gluon spliting"}; // true:: remove gluon spliting
  Configurable<bool> searchUpToQuark{"searchUpToQuark", true, "Finding first mother in particles to quark"};
  Configurable<bool> useJetProb{"useJetProb", false, "fill table for track counting algorithm"};
  Configurable<bool> trackProbQA{"trackProbQA", false, "fill track probability histograms separately for geometric positive and negative tracks for QA"};
  Configurable<bool> doSV{"doSV", false, "fill table for secondary vertex algorithm"};
  Configurable<int> numCount{"numCount", 3, "number of track counting"};
  Configurable<int> resoFuncMatching{"resoFuncMatching", 0, "matching parameters of resolution function as MC samble (0: custom, 1: custom & inc, 2: MB, 3: MB & inc, 4: JJ, 5: JJ & inc)"};
  Configurable<std::vector<float>> paramsResoFuncData{"paramsResoFuncData", std::vector<float>{1306800, -0.1049, 0.861425, 13.7547, 0.977967, 8.96823, 0.151595, 6.94499, 0.0250301}, "parameters of gaus(0)+expo(3)+expo(5)+expo(7))"};
  Configurable<std::vector<float>> paramsResoFuncIncJetMC{"paramsResoFuncIncJetMC", std::vector<float>{1908803.027, -0.059, 0.895, 13.467, 1.005, 8.867, 0.098, 6.929, 0.011}, "parameters of gaus(0)+expo(3)+expo(5)+expo(7)))"};
  Configurable<std::vector<float>> paramsResoFuncCharmJetMC{"paramsResoFuncCharmJetMC", std::vector<float>{282119.753, -0.065, 0.893, 11.608, 0.945, 8.029, 0.131, 6.244, 0.027}, "parameters of gaus(0)+expo(3)+expo(5)+expo(7)))"};
  Configurable<std::vector<float>> paramsResoFuncBeautyJetMC{"paramsResoFuncBeautyJetMC", std::vector<float>{74901.583, -0.082, 0.874, 10.332, 0.941, 7.352, 0.097, 6.220, 0.022}, "parameters of gaus(0)+expo(3)+expo(5)+expo(7)))"};
  Configurable<std::vector<float>> paramsResoFuncLfJetMC{"paramsResoFuncLfJetMC", std::vector<float>{1539435.343, -0.061, 0.896, 13.272, 1.034, 5.884, 0.004, 7.843, 0.090}, "parameters of gaus(0)+expo(3)+expo(5)+expo(7)))"};
  Configurable<float> minSignImpXYSig{"minsIPs", -40.0, "minimum of signed impact parameter significance"};
  Configurable<float> tagPoint{"tagPoint", 2.5, "tagging working point"};

  ConfigurableAxis binTrackProbability{"binTrackProbability", {100, 0.f, 1.f}, ""};
  ConfigurableAxis binJetFlavour{"binJetFlavour", {6, -0.5, 5.5}, ""};

  using JetTagTracksData = soa::Join<JetTracks, aod::JTrackPIs, aod::JTracksTag>;
  using JetTagTracksMCD = soa::Join<JetTracksMCD, aod::JTrackPIs, aod::JTracksTag>;
  using OriTracksData = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TracksDCACov, aod::TrackSelection>;
  using OriTracksMCD = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TracksDCACov, aod::TrackSelection, aod::McTrackLabels>;

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

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};
  void init(InitContext const&)
  {
    maxOrder = numCount + 1; // 0: untagged, >1 : N ordering

    // Set up the resolution function
    resoFuncMatch = resoFuncMatching;
    switch (resoFuncMatch) {
      case 0:
        vecParamsData = (std::vector<float>)paramsResoFuncData;
        vecParamsIncJetMC = (std::vector<float>)paramsResoFuncIncJetMC;
        LOG(info) << "defined parameters of resolution function: custom";
        break;
      case 1:
        vecParamsData = (std::vector<float>)paramsResoFuncData;
        vecParamsCharmJetMC = (std::vector<float>)paramsResoFuncCharmJetMC;
        vecParamsBeautyJetMC = (std::vector<float>)paramsResoFuncBeautyJetMC;
        vecParamsLfJetMC = (std::vector<float>)paramsResoFuncLfJetMC;
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
        vecParamsCharmJetMC = {281446.003, -0.063, 0.894, 11.598, 0.943, 8.025, 0.130, 6.227, 0.027};
        vecParamsBeautyJetMC = {74839.065, -0.081, 0.875, 10.314, 0.939, 7.326, 0.101, 6.309, 0.024};
        vecParamsLfJetMC = {1531580.038, -0.062, -0.896, 13.267, 1.034, 5.866, 0.004, 7.836, 0.090};
        LOG(info) << "defined parameters of resolution function: PYTHIA8, JJ, LHC23d4";
        break;
      case 5: // TODO
        vecParamsData = (std::vector<float>)paramsResoFuncData;
        vecParamsIncJetMC = {1900387.527, -0.059, 0.895, 13.461, 1.004, 8.860, 0.098, 6.931, 0.011};
        LOG(info) << "defined parameters of resolution function: PYTHIA8, JJ, LHC23d4 & use inclusive distribution";
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
      registry.add("h2_pos_track_probability_flavour", "positive track probability", {HistType::kTH2F, {{trackProbabilityAxis}, {jetFlavourAxis}}});
      registry.add("h2_neg_track_probability_flavour", "negative track probability", {HistType::kTH2F, {{trackProbabilityAxis}, {jetFlavourAxis}}});
    }
  }

  void processDummy(JetCollisions const&)
  {
  }
  PROCESS_SWITCH(JetTaggerHFTask, processDummy, "Dummy process", true);

  void processData(JetCollision const& collision, JetTableData const& jets, JetTagTracksData const& jtracks, OriTracksData const& tracks)
  {
    for (auto& jet : jets) {
      int algorithm2 = 0;
      int algorithm3 = 0;
      if (useJetProb) {
        jetProb.clear();
        jetProb.reserve(maxOrder);
        for (int order = 0; order < maxOrder; order++) {
          jetProb.push_back(jettaggingutilities::getJetProbability(fSignImpXYSigData, collision, jet, jtracks, tracks, order, tagPoint, minSignImpXYSig));
        }
      }
      // if (doSV) algorithm2 = jettaggingutilities::Algorithm2((mcdjet, tracks);
      taggingTableData(0, jetProb, algorithm2, algorithm3);
    }
  }
  PROCESS_SWITCH(JetTaggerHFTask, processData, "Fill tagging decision for data jets", false);

  void processMCD(JetCollision const& collision, JetTableMCD const& mcdjets, JetTagTracksMCD const& jtracks, OriTracksMCD const& tracks, JetParticles const& particles)
  {
    for (auto& mcdjet : mcdjets) {
      typename JetTagTracksMCD::iterator hftrack;
      int origin = 0;
      if (removeGluonShower)
        origin = jettaggingutilities::mcdJetFromHFShower(mcdjet, jtracks, particles, maxDeltaR, searchUpToQuark);
      else
        origin = jettaggingutilities::jetTrackFromHFShower(mcdjet, jtracks, particles, hftrack, searchUpToQuark);
      int algorithm2 = 0;
      int algorithm3 = 0;
      if (useJetProb) {
        jetProb.clear();
        jetProb.reserve(maxOrder);
        for (int order = 0; order < maxOrder; order++) {
          if (useResoFuncFromIncJet) {
            jetProb.push_back(jettaggingutilities::getJetProbability(fSignImpXYSigIncJetMC, collision, mcdjet, jtracks, tracks, order, tagPoint, minSignImpXYSig));
          } else {
            if (origin == JetTaggingSpecies::charm) {
              jetProb.push_back(jettaggingutilities::getJetProbability(fSignImpXYSigCharmJetMC, collision, mcdjet, jtracks, tracks, order, tagPoint, minSignImpXYSig));
            }
            if (origin == JetTaggingSpecies::beauty) {
              jetProb.push_back(jettaggingutilities::getJetProbability(fSignImpXYSigBeautyJetMC, collision, mcdjet, jtracks, tracks, order, tagPoint, minSignImpXYSig));
            }
            if (origin == JetTaggingSpecies::lightflavour) {
              jetProb.push_back(jettaggingutilities::getJetProbability(fSignImpXYSigLfJetMC, collision, mcdjet, jtracks, tracks, order, tagPoint, minSignImpXYSig));
            }
          }
        }
        if (trackProbQA) {
          for (auto& jtrack : mcdjet.template tracks_as<JetTagTracksMCD>()) {
            auto track = jtrack.template track_as<OriTracksMCD>();
            auto geoSign = jettaggingutilities::getGeoSign(collision, mcdjet, track);
            float probTrack = -1;
            if (useResoFuncFromIncJet) {
              probTrack = jettaggingutilities::getTrackProbability(fSignImpXYSigIncJetMC, track, minSignImpXYSig);
            } else {
              if (origin == JetTaggingSpecies::charm) {
                probTrack = jettaggingutilities::getTrackProbability(fSignImpXYSigCharmJetMC, track, minSignImpXYSig);
              }
              if (origin == JetTaggingSpecies::beauty) {
                probTrack = jettaggingutilities::getTrackProbability(fSignImpXYSigBeautyJetMC, track, minSignImpXYSig);
              }
              if (origin == JetTaggingSpecies::lightflavour) {
                probTrack = jettaggingutilities::getTrackProbability(fSignImpXYSigLfJetMC, track, minSignImpXYSig);
              }
            }
            if (geoSign > 0) {
              registry.fill(HIST("h2_pos_track_probability_flavour"), probTrack, origin);
            } else {
              registry.fill(HIST("h2_neg_track_probability_flavour"), probTrack, origin);
            }
          }
        }
      }
      // if (doSV) algorithm2 = jettaggingutilities::Algorithm2((mcdjet, tracks);
      taggingTableMCD(origin, jetProb, algorithm2, algorithm3);
    }
  }
  PROCESS_SWITCH(JetTaggerHFTask, processMCD, "Fill tagging decision for mcd jets", false);

  void processTraining(JetCollision const& /*collision*/, JetTableMCD const& /*mcdjets*/, JetTagTracksMCD const& /*tracks*/)
  {
    // To create table for ML
  }
  PROCESS_SWITCH(JetTaggerHFTask, processTraining, "Fill tagging decision for mcd jets", false);
};

struct JetTaggerHFExtTask {

  Produces<aod::JTracksTag> jTracksTagTable;

  void init(InitContext const&)
  {
  }

  void processTracks(soa::Join<aod::Tracks, aod::TracksCov, aod::TrackSelection, aod::TracksDCA, aod::TracksDCACov>::iterator const& track)
  {
    float sigmaDcaXYZ2 = 0;
    float dcaXYZ = getDcaXYZ(track, &sigmaDcaXYZ2);

    jTracksTagTable(dcaXYZ, sigmaDcaXYZ2);
  }
  PROCESS_SWITCH(JetTaggerHFExtTask, processTracks, "produces derived track table for tagging", true);
};

using JetTaggerChargedJets = JetTaggerHFTask<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>, soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>, aod::ChargedJetTags, aod::ChargedMCDetectorLevelJetTags>;
using JetTaggerFullJets = JetTaggerHFTask<soa::Join<aod::FullJets, aod::FullJetConstituents>, soa::Join<aod::FullMCDetectorLevelJets, aod::FullMCDetectorLevelJetConstituents>, aod::FullJetTags, aod::FullMCDetectorLevelJetTags>;
// using JetTaggerNeutralJets = JetTaggerHFTask<soa::Join<aod::NeutralJets, aod::NeutralJetConstituents>,soa::Join<aod::NeutralMCDetectorLevelJets, aod::NeutralMCDetectorLevelJetConstituents>, aod::NeutralJetTags, aod::NeutralMCDetectorLevelJetTags>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{

  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(
    adaptAnalysisTask<JetTaggerHFExtTask>(cfgc,
                                          SetDefaultProcesses{}, TaskName{"jet-taggerhf-extension"}));

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
