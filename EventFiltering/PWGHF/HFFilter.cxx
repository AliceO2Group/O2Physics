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
// O2 includes

/// \file HFFilter.cxx
/// \brief task for selection of events with HF signals
///
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN
/// \author Marcel Lesch <marcel.lesch@tum.de>, TUM
/// \author Alexandre Bigot <alexandre.bigot@cern.ch>, Strasbourg University
/// \author Biao Zhang <biao.zhang@cern.ch>, CCNU
/// \author Federica Zanone <federica.zanone@cern.ch>, Heidelberg University

#include <onnxruntime/core/session/experimental_onnxruntime_cxx_api.h> // needed for HFFilterHelpers, to be fixed

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "EventFiltering/filterTables.h"
#include "EventFiltering/PWGHF/HFFilterHelpers.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::aod::hffilters;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct HfFilter { // Main struct for HF triggers

  Produces<aod::HfFilters> tags;
  Produces<aod::HFOptimisationTreeBeauty> optimisationTreeBeauty;
  Produces<aod::HFOptimisationTreeCharm> optimisationTreeCharm;
  Produces<aod::HFOptimisationTreeFemto> optimisationTreeFemto;
  Produces<aod::HFOptimisationTreeCollisions> optimisationTreeCollisions;

  Configurable<int> activateQA{"activateQA", 0, "flag to enable QA histos (0 no QA, 1 basic QA, 2 extended QA, 3 very extended QA)"};
  Configurable<bool> applyEventSelection{"applyEventSelection", false, "flag to enable event selection (sel8 + Zvt)"};

  // parameters for all triggers
  // nsigma PID (except for V0 and cascades)
  Configurable<LabeledArray<float>> nSigmaPidCuts{"nSigmaPidCuts", {cutsNsigma[0], 3, 5, labelsRowsNsigma, labelsColumnsNsigma}, "Nsigma cuts for TPC/TOF PID (except for V0 and cascades)"};
  // min pts for tracks and bachelors (except for V0 and cascades)
  Configurable<LabeledArray<float>> minPtCuts{"minPtCuts", {cutsMinPt[0], 1, 4, labelsEmpty, labelsColumnsMinPt}, "minimum pT for bachelor tracks (except for V0 and cascades)"};

  // parameters for high-pT triggers
  Configurable<LabeledArray<float>> ptThresholds{"ptThresholds", {cutsHighPtThresholds[0], 1, 2, labelsEmpty, labelsColumnsHighPtThresholds}, "pT treshold for high pT charm hadron candidates for kHighPt triggers in GeV/c"};

  // parameters for beauty triggers
  Configurable<LabeledArray<float>> deltaMassBeauty{"deltaMassBeauty", {cutsDeltaMassB[0], 1, kNBeautyParticles + 1, labelsEmpty, labelsColumnsDeltaMassB}, "invariant-mass delta with respect to the b-hadron masses in GeV/c2"};
  Configurable<std::vector<double>> pTBinsTrack{"pTBinsTrack", std::vector<double>{hf_cuts_single_track::vecBinsPtTrack}, "track pT bin limits for DCAXY pT-dependent cut"};
  Configurable<LabeledArray<double>> cutsTrackBeauty3Prong{"cutsTrackBeauty3Prong", {hf_cuts_single_track::cutsTrack[0], hf_cuts_single_track::nBinsPtTrack, hf_cuts_single_track::nCutVarsTrack, hf_cuts_single_track::labelsPtTrack, hf_cuts_single_track::labelsCutVarTrack}, "Single-track selections per pT bin for 3-prong beauty candidates"};
  Configurable<LabeledArray<double>> cutsTrackBeauty4Prong{"cutsTrackBeauty4Prong", {hf_cuts_single_track::cutsTrack[0], hf_cuts_single_track::nBinsPtTrack, hf_cuts_single_track::nCutVarsTrack, hf_cuts_single_track::labelsPtTrack, hf_cuts_single_track::labelsCutVarTrack}, "Single-track selections per pT bin for 4-prong beauty candidates"};

  // parameters for femto triggers
  Configurable<float> femtoMaxRelativeMomentum{"femtoMaxRelativeMomentum", 2., "Maximal allowed value for relative momentum between charm-proton pairs in GeV/c"};
  Configurable<LabeledArray<int>> enableFemtoChannels{"enableFemtoChannels", {activeFemtoChannels[0], 1, 5, labelsEmpty, labelsColumnsFemtoChannels}, "Flags to enable/disable femto channels"};
  Configurable<bool> requireCharmMassForFemto{"requireCharmMassForFemto", false, "Flags to enable/disable cut on charm-hadron invariant-mass window for femto"};
  Configurable<float> ptThresholdForFemtoPid{"ptThresholdForFemtoPid", 8., "pT threshold for changing strategy of proton PID in femto triggers"};

  // double charm
  Configurable<LabeledArray<int>> enableDoubleCharmChannels{"enableDoubleCharmChannels", {activeDoubleCharmChannels[0], 1, 3, labelsEmpty, labelsColumnsDoubleCharmChannels}, "Flags to enable/disable double charm channels"};

  // parameters for V0 + charm triggers
  Configurable<LabeledArray<float>> cutsGammaK0sLambda{"cutsGammaK0sLambda", {cutsV0s[0], 1, 6, labelsEmpty, labelsColumnsV0s}, "Selections for V0s (gamma, K0s, Lambda) for D+V0 triggers"};
  Configurable<LabeledArray<float>> maxDeltaMassCharmReso{"maxDeltaMassCharmReso", {cutsMassCharmReso[0], 1, 6, labelsEmpty, labelsColumnsDeltaMasseCharmReso}, "maximum invariant-mass delta for charm hadron resonances in GeV/c2"};

  // parameters for charm baryons to Xi bachelor
  Configurable<LabeledArray<float>> cutsXiCascades{"cutsXiCascades", {cutsCascades[0], 1, 7, labelsEmpty, labelsColumnsCascades}, "Selections for cascades (Xi) for Xi+bachelor triggers"};
  Configurable<LabeledArray<float>> cutsXiBachelor{"cutsXiBachelor", {cutsCharmBaryons[0], 1, 4, labelsEmpty, labelsColumnsCharmBaryons}, "Selections for charm baryons (Xi+Pi and Xi+Ka)"};
  Configurable<LabeledArray<double>> cutsTrackCharmBaryonBachelor{"cutsTrackCharmBaryonBachelor", {hf_cuts_single_track::cutsTrack[0], hf_cuts_single_track::nBinsPtTrack, hf_cuts_single_track::nCutVarsTrack, hf_cuts_single_track::labelsPtTrack, hf_cuts_single_track::labelsCutVarTrack}, "Single-track selections per pT bin for charm-baryon bachelor candidates"};

  // parameters for ML application with ONNX
  Configurable<bool> applyML{"applyML", false, "Flag to enable or disable ML application"};
  Configurable<std::vector<double>> pTBinsBDT{"pTBinsBDT", std::vector<double>{hf_cuts_bdt_multiclass::vecBinsPt}, "track pT bin limits for BDT cut"};

  Configurable<std::string> onnxFileD0ToKPiConf{"onnxFileD0ToKPiConf", "XGBoostModel.onnx", "ONNX file for ML model for D0 candidates"};
  Configurable<LabeledArray<double>> thresholdBDTScoreD0ToKPi{"thresholdBDTScoreD0ToKPi", {hf_cuts_bdt_multiclass::cuts[0], hf_cuts_bdt_multiclass::nBinsPt, hf_cuts_bdt_multiclass::nCutBdtScores, hf_cuts_bdt_multiclass::labelsPt, hf_cuts_bdt_multiclass::labelsCutBdt}, "Threshold values for BDT output scores of D0 candidates"};
  Configurable<std::string> onnxFileDPlusToPiKPiConf{"onnxFileDPlusToPiKPiConf", "", "ONNX file for ML model for D+ candidates"};
  Configurable<LabeledArray<double>> thresholdBDTScoreDPlusToPiKPi{"thresholdBDTScoreDPlusToPiKPi", {hf_cuts_bdt_multiclass::cuts[0], hf_cuts_bdt_multiclass::nBinsPt, hf_cuts_bdt_multiclass::nCutBdtScores, hf_cuts_bdt_multiclass::labelsPt, hf_cuts_bdt_multiclass::labelsCutBdt}, "Threshold values for BDT output scores of D+ candidates"};
  Configurable<std::string> onnxFileDSToPiKKConf{"onnxFileDSToPiKKConf", "", "ONNX file for ML model for Ds+ candidates"};
  Configurable<LabeledArray<double>> thresholdBDTScoreDSToPiKK{"thresholdBDTScoreDSToPiKK", {hf_cuts_bdt_multiclass::cuts[0], hf_cuts_bdt_multiclass::nBinsPt, hf_cuts_bdt_multiclass::nCutBdtScores, hf_cuts_bdt_multiclass::labelsPt, hf_cuts_bdt_multiclass::labelsCutBdt}, "Threshold values for BDT output scores of Ds+ candidates"};
  Configurable<std::string> onnxFileLcToPiKPConf{"onnxFileLcToPiKPConf", "", "ONNX file for ML model for Lc+ candidates"};
  Configurable<LabeledArray<double>> thresholdBDTScoreLcToPiKP{"thresholdBDTScoreLcToPiKP", {hf_cuts_bdt_multiclass::cuts[0], hf_cuts_bdt_multiclass::nBinsPt, hf_cuts_bdt_multiclass::nCutBdtScores, hf_cuts_bdt_multiclass::labelsPt, hf_cuts_bdt_multiclass::labelsCutBdt}, "Threshold values for BDT output scores of Lc+ candidates"};
  Configurable<std::string> onnxFileXicToPiKPConf{"onnxFileXicToPiKPConf", "", "ONNX file for ML model for Xic+ candidates"};
  Configurable<LabeledArray<double>> thresholdBDTScoreXicToPiKP{"thresholdBDTScoreXicToPiKP", {hf_cuts_bdt_multiclass::cuts[0], hf_cuts_bdt_multiclass::nBinsPt, hf_cuts_bdt_multiclass::nCutBdtScores, hf_cuts_bdt_multiclass::labelsPt, hf_cuts_bdt_multiclass::labelsCutBdt}, "Threshold values for BDT output scores of Xic+ candidates"};

  // CCDB configuration
  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> mlModelPathCCDB{"mlModelPathCCDB", "Analysis/PWGHF/ML/HFTrigger/", "Path on CCDB"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB. Exceptions: > 0 for the specific timestamp, 0 gets the run dependent timestamp"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};
  int currentRun{0}; // needed to detect if the run changed and trigger update of calibrations etc.

  // TPC PID calibrations
  Configurable<int> setTPCCalib{"setTPCCalib", 0, "0 is not use re-calibrations, 1 is compute TPC post-calibrated n-sigmas, 2 is using TPC Spline"};
  Configurable<std::string> ccdbBBProton{"ccdbBBProton", "Users/l/lserksny/PIDProton", "Path to the CCDB ocject for proton BB param"};
  Configurable<std::string> ccdbBBAntiProton{"ccdbBBAntiProton", "Users/l/lserksny/PIDAntiProton", "Path to the CCDB ocject for antiproton BB param"};
  Configurable<std::string> ccdbBBPion{"ccdbBBPion", "Users/l/lserksny/PIDPion", "Path to the CCDB ocject for Pion BB param"};
  Configurable<std::string> ccdbBBAntiPion{"ccdbBBAntiPion", "Users/l/lserksny/PIDAntiPion", "Path to the CCDB ocject for antiPion BB param"};
  Configurable<std::string> ccdbBBKaon{"ccdbBBKaon", "Users/l/lserksny/PIDPion", "Path to the CCDB ocject for Kaon BB param"};
  Configurable<std::string> ccdbBBAntiKaon{"ccdbBBAntiKaon", "Users/l/lserksny/PIDAntiPion", "Path to the CCDB ocject for antiKaon BB param"};
  Configurable<string> ccdbPathTPC{"ccdbPathTPC", "Users/i/iarsene/Calib/TPCpostCalib", "base path to the CCDB object"};

  // parameter for Optimisation Tree
  Configurable<bool> applyOptimisation{"applyOptimisation", false, "Flag to enable or disable optimisation"};

  // array of ONNX config and BDT thresholds
  std::array<std::string, kNCharmParticles> onnxFiles;
  std::array<LabeledArray<double>, kNCharmParticles> thresholdBDTScores;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  std::shared_ptr<TH1> hProcessedEvents;

  // QA histos
  std::shared_ptr<TH1> hN2ProngCharmCand, hN3ProngCharmCand;
  std::array<std::shared_ptr<TH1>, kNCharmParticles> hCharmHighPt{};
  std::array<std::shared_ptr<TH1>, kNCharmParticles> hCharmProtonKstarDistr{};
  std::array<std::shared_ptr<TH2>, kNBeautyParticles> hMassVsPtB{};
  std::array<std::shared_ptr<TH2>, kNCharmParticles + 8> hMassVsPtC{}; // +6 for resonances (D*+, D*0, Ds*+, Ds1+, Ds2*+, Xic*) +2 for charm baryons (Xi+Pi, Xi+Ka)
  std::shared_ptr<TH2> hProtonTPCPID, hProtonTOFPID;
  std::array<std::shared_ptr<TH1>, kNCharmParticles> hBDTScoreBkg{};
  std::array<std::shared_ptr<TH1>, kNCharmParticles> hBDTScorePrompt{};
  std::array<std::shared_ptr<TH1>, kNCharmParticles> hBDTScoreNonPrompt{};
  std::array<std::shared_ptr<TH2>, kNV0> hArmPod{};
  std::shared_ptr<TH2> hV0Selected;
  std::shared_ptr<TH1> hMassXi;

  // ONNX
  std::array<std::shared_ptr<Ort::Experimental::Session>, kNCharmParticles> sessionML = {nullptr, nullptr, nullptr, nullptr, nullptr};
  std::array<std::vector<std::vector<int64_t>>, kNCharmParticles> inputShapesML{};
  std::array<Ort::Env, kNCharmParticles> envML = {
    Ort::Env{ORT_LOGGING_LEVEL_ERROR, "ml-model-d0-triggers"},
    Ort::Env{ORT_LOGGING_LEVEL_ERROR, "ml-model-dplus-triggers"},
    Ort::Env{ORT_LOGGING_LEVEL_ERROR, "ml-model-ds-triggers"},
    Ort::Env{ORT_LOGGING_LEVEL_ERROR, "ml-model-lc-triggers"},
    Ort::Env{ORT_LOGGING_LEVEL_ERROR, "ml-model-xic-triggers"}};
  std::array<Ort::SessionOptions, kNCharmParticles> sessionOptions{Ort::SessionOptions(), Ort::SessionOptions(), Ort::SessionOptions(), Ort::SessionOptions(), Ort::SessionOptions()};
  std::array<int, kNCharmParticles> dataTypeML{};

  // material correction for track propagation
  o2::base::MatLayerCylSet* lut;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  o2::base::Propagator::MatCorrType noMatCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  // helper object
  HfFilterHelper helper;

  void init(InitContext&)
  {
    helper.setPtBinsSingleTracks(pTBinsTrack);
    helper.setMinPtBeautyBachelor(minPtCuts->get(0u, 0u));
    helper.setMinPtDstarSoftPion(minPtCuts->get(0u, 1u));
    helper.setCutsSingleTrackBeauty(cutsTrackBeauty3Prong, cutsTrackBeauty4Prong);
    helper.setMinPtProtonForFemto(minPtCuts->get(0u, 2u));
    helper.setPtThresholdPidStrategyForFemto(ptThresholdForFemtoPid);
    helper.setNsigmaProtonCutsForFemto(std::array{nSigmaPidCuts->get(0u, 3u), nSigmaPidCuts->get(1u, 3u), nSigmaPidCuts->get(2u, 3u)});
    helper.setNsigmaProtonCutsForCharmBaryons(nSigmaPidCuts->get(0u, 0u), nSigmaPidCuts->get(1u, 0u));
    helper.setNsigmaPionKaonCutsForDzero(nSigmaPidCuts->get(0u, 1u), nSigmaPidCuts->get(1u, 1u));
    helper.setNsigmaKaonCutsFor3Prongs(nSigmaPidCuts->get(0u, 2u), nSigmaPidCuts->get(1u, 2u));
    helper.setDeltaMassCharmHadForBeauty(deltaMassBeauty->get(0u, kNBeautyParticles));
    helper.setV0Selections(cutsGammaK0sLambda->get(0u, 0u), cutsGammaK0sLambda->get(0u, 1u), cutsGammaK0sLambda->get(0u, 2u), cutsGammaK0sLambda->get(0u, 3u), cutsGammaK0sLambda->get(0u, 4u), cutsGammaK0sLambda->get(0u, 5u));
    helper.setXiSelections(cutsXiCascades->get(0u, 0u), cutsXiCascades->get(0u, 1u), cutsXiCascades->get(0u, 2u), cutsXiCascades->get(0u, 3u), cutsXiCascades->get(0u, 4u), cutsXiCascades->get(0u, 5u), cutsXiCascades->get(0u, 6u));
    helper.setCutsSingleTrackCharmBaryonBachelor(cutsTrackCharmBaryonBachelor);
    helper.setMinPtCharmBaryonBachelor(minPtCuts->get(0u, 3u));
    helper.setNsigmaPiCutsForCharmBaryonBachelor(nSigmaPidCuts->get(0u, 4u), nSigmaPidCuts->get(1u, 4u));
    helper.setTpcPidCalibrationOption(setTPCCalib);

    hProcessedEvents = registry.add<TH1>("fProcessedEvents", "HF - event filtered;;counts", HistType::kTH1F, {{kNtriggersHF + 2, -0.5, kNtriggersHF + 1.5}});
    for (auto iBin = 0; iBin < kNtriggersHF + 2; ++iBin) {
      hProcessedEvents->GetXaxis()->SetBinLabel(iBin + 1, eventTitles[iBin].data());
    }

    if (activateQA) {
      hN2ProngCharmCand = registry.add<TH1>("fN2ProngCharmCand", "Number of 2-prong charm candidates per event;#it{N}_{candidates};counts", HistType::kTH1F, {{50, -0.5, 49.5}});
      hN3ProngCharmCand = registry.add<TH1>("fN3ProngCharmCand", "Number of 3-prong charm candidates per event;#it{N}_{candidates};counts", HistType::kTH1F, {{50, -0.5, 49.5}});
      for (int iCharmPart{0}; iCharmPart < kNCharmParticles; ++iCharmPart) {
        hCharmHighPt[iCharmPart] = registry.add<TH1>(Form("f%sHighPt", charmParticleNames[iCharmPart].data()), Form("#it{p}_{T} distribution of triggered high-#it{p}_{T} %s candidates;#it{p}_{T} (GeV/#it{c});counts", charmParticleNames[iCharmPart].data()), HistType::kTH1F, {ptAxis});
        hCharmProtonKstarDistr[iCharmPart] = registry.add<TH1>(Form("f%sProtonKstarDistr", charmParticleNames[iCharmPart].data()), Form("#it{k}* distribution of triggered p#minus%s pairs;#it{k}* (GeV/#it{c});counts", charmParticleNames[iCharmPart].data()), HistType::kTH1F, {kstarAxis});
        hMassVsPtC[iCharmPart] = registry.add<TH2>(Form("fMassVsPt%s", charmParticleNames[iCharmPart].data()), Form("#it{M} vs. #it{p}_{T} distribution of triggered %s candidates;#it{p}_{T} (GeV/#it{c});#it{M} (GeV/#it{c}^{2});counts", charmParticleNames[iCharmPart].data()), HistType::kTH2F, {ptAxis, massAxisC[iCharmPart]});
        if (applyML && activateQA > 1) {
          hBDTScoreBkg[iCharmPart] = registry.add<TH1>(Form("f%sBDTScoreBkgDistr", charmParticleNames[iCharmPart].data()), Form("BDT background score distribution for %s;BDT background score;counts", charmParticleNames[iCharmPart].data()), HistType::kTH1F, {bdtAxis});
          hBDTScorePrompt[iCharmPart] = registry.add<TH1>(Form("f%sBDTScorePromptDistr", charmParticleNames[iCharmPart].data()), Form("BDT prompt score distribution for %s;BDT prompt score;counts", charmParticleNames[iCharmPart].data()), HistType::kTH1F, {bdtAxis});
          hBDTScoreNonPrompt[iCharmPart] = registry.add<TH1>(Form("f%sBDTScoreNonPromptDistr", charmParticleNames[iCharmPart].data()), Form("BDT nonprompt score distribution for %s;BDT nonprompt score;counts", charmParticleNames[iCharmPart].data()), HistType::kTH1F, {bdtAxis});
        }
      }
      // charm resonances
      hMassVsPtC[kNCharmParticles] = registry.add<TH2>("fMassVsPtDStarPlus", "#Delta#it{M} vs. #it{p}_{T} distribution of triggered DStarPlus candidates;#it{p}_{T} (GeV/#it{c});#Delta#it{M} (GeV/#it{c}^{2});counts", HistType::kTH2F, {ptAxis, massAxisC[kNCharmParticles]});
      hMassVsPtC[kNCharmParticles + 1] = registry.add<TH2>("fMassVsPtDStarZero", "#Delta#it{M} vs. #it{p}_{T} distribution of triggered DStarZero candidates;#it{p}_{T} (GeV/#it{c});#Delta#it{M} (GeV/#it{c}^{2});counts", HistType::kTH2F, {ptAxis, massAxisC[kNCharmParticles + 1]});
      hMassVsPtC[kNCharmParticles + 2] = registry.add<TH2>("fMassVsPtDStarS", "#Delta#it{M} vs. #it{p}_{T} distribution of triggered DStarS candidates;#it{p}_{T} (GeV/#it{c});#Delta#it{M} (GeV/#it{c}^{2});counts", HistType::kTH2F, {ptAxis, massAxisC[kNCharmParticles + 2]});
      hMassVsPtC[kNCharmParticles + 3] = registry.add<TH2>("fMassVsPtDs1Plus", "#Delta#it{M} vs. #it{p}_{T} distribution of triggered Ds1Plus candidates;#it{p}_{T} (GeV/#it{c});#Delta#it{M} (GeV/#it{c}^{2});counts", HistType::kTH2F, {ptAxis, massAxisC[kNCharmParticles + 3]});
      hMassVsPtC[kNCharmParticles + 4] = registry.add<TH2>("fMassVsPtDs2StarPlus", "#Delta#it{M} vs. #it{p}_{T} distribution of triggered Ds2StarPlus candidates;#it{p}_{T} (GeV/#Delta#it{c});#it{M} (GeV/#it{c}^{2});counts", HistType::kTH2F, {ptAxis, massAxisC[kNCharmParticles + 4]});
      hMassVsPtC[kNCharmParticles + 5] = registry.add<TH2>("fMassVsPtXicStar", "#Delta#it{M} vs. #it{p}_{T} distribution of triggered XicStar candidates;#it{p}_{T} (GeV/#it{c});#Delta#it{M} (GeV/#it{c}^{2});counts", HistType::kTH2F, {ptAxis, massAxisC[kNCharmParticles + 5]});
      // charm baryons to LF cascades
      hMassVsPtC[kNCharmParticles + 6] = registry.add<TH2>("fMassVsPtCharmBaryonToXiPi", "#it{M} vs. #it{p}_{T} distribution of triggered #Xi+#pi candidates;#it{p}_{T} (GeV/#it{c});#it{M} (GeV/#it{c}^{2});counts", HistType::kTH2F, {ptAxis, massAxisC[kNCharmParticles + 6]});
      hMassVsPtC[kNCharmParticles + 7] = registry.add<TH2>("fMassVsPtCharmBaryonToXiKa", "#it{M} vs. #it{p}_{T} distribution of triggered #Xi+K candidates;#it{p}_{T} (GeV/#it{c});#it{M} (GeV/#it{c}^{2});counts", HistType::kTH2F, {ptAxis, massAxisC[kNCharmParticles + 7]});
      for (int iBeautyPart{0}; iBeautyPart < kNBeautyParticles; ++iBeautyPart) {
        hMassVsPtB[iBeautyPart] = registry.add<TH2>(Form("fMassVsPt%s", beautyParticleNames[iBeautyPart].data()), Form("#it{M} vs. #it{p}_{T} distribution of triggered %s candidates;#it{p}_{T} (GeV/#it{c});#it{M} (GeV/#it{c}^{2});counts", beautyParticleNames[iBeautyPart].data()), HistType::kTH2F, {ptAxis, massAxisB[iBeautyPart]});
      }
      for (int iV0{kPhoton}; iV0 < kNV0; ++iV0) {
        hArmPod[iV0] = registry.add<TH2>(Form("fArmPod%s", v0Names[iV0].data()), Form("Armenteros Podolanski plot for selected %s;#it{#alpha};#it{q}_{T} (GeV/#it{c})", v0Labels[iV0].data()), HistType::kTH2F, {alphaAxis, qtAxis});
      }
      hMassXi = registry.add<TH1>("fMassXi", "#it{M} distribution of #Xi candidates;#it{M} (GeV/#it{c}^{2});counts", HistType::kTH1F, {{100, 1.28f, 1.36f}});

      if (activateQA > 1) {
        hProtonTPCPID = registry.add<TH2>("fProtonTPCPID", "#it{N}_{#sigma}^{TPC} vs. #it{p} for selected protons;#it{p} (GeV/#it{c});#it{N}_{#sigma}^{TPC}", HistType::kTH2F, {pAxis, nSigmaAxis});
        hProtonTOFPID = registry.add<TH2>("fProtonTOFPID", "#it{N}_{#sigma}^{TOF} vs. #it{p} for selected protons;#it{p} (GeV/#it{c});#it{N}_{#sigma}^{TOF}", HistType::kTH2F, {pAxis, nSigmaAxis});
        hV0Selected = registry.add<TH2>("fV0Selected", "Selections for V0s;;counts", HistType::kTH2F, {{10, -0.5, 9.5}, {kNV0, -0.5, kNV0 - 0.5}});

        for (int iV0{kPhoton}; iV0 < kNV0; ++iV0) {
          hV0Selected->GetYaxis()->SetBinLabel(iV0 + 1, v0Labels[iV0].data());
        }
        hV0Selected->GetXaxis()->SetBinLabel(1, "analysed");
        hV0Selected->GetXaxis()->SetBinLabel(2, "rej. |#eta|");
        hV0Selected->GetXaxis()->SetBinLabel(3, "rej. radius");
        hV0Selected->GetXaxis()->SetBinLabel(4, "rej. cos(#theta_{P})");
        hV0Selected->GetXaxis()->SetBinLabel(5, "rej. AP / Mass");
        hV0Selected->GetXaxis()->SetBinLabel(6, "rej. DCA V0");
        hV0Selected->GetXaxis()->SetBinLabel(7, "rej. DCA V0 daughters");
        hV0Selected->GetXaxis()->SetBinLabel(8, "rej. pair cut");
        hV0Selected->GetXaxis()->SetBinLabel(9, "rej. PID");
        hV0Selected->GetXaxis()->SetBinLabel(10, "selected");
      }
    }

    ccdb->setURL(url.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    ccdbApi.init(url);
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>("GLO/Param/MatLUT"));

    thresholdBDTScores = {
      thresholdBDTScoreD0ToKPi,
      thresholdBDTScoreDPlusToPiKPi,
      thresholdBDTScoreDSToPiKK,
      thresholdBDTScoreLcToPiKP,
      thresholdBDTScoreXicToPiKP};

    onnxFiles = {
      onnxFileD0ToKPiConf,
      onnxFileDPlusToPiKPiConf,
      onnxFileDSToPiKKConf,
      onnxFileLcToPiKPConf,
      onnxFileXicToPiKPConf};

    // init ONNX runtime session
    if (applyML && (!loadModelsFromCCDB || timestampCCDB != 0)) {
      for (auto iCharmPart{0}; iCharmPart < kNCharmParticles; ++iCharmPart) {
        if (onnxFiles[iCharmPart] != "") {
          sessionML[iCharmPart].reset(helper.initONNXSession(onnxFiles[iCharmPart], charmParticleNames[iCharmPart], envML[iCharmPart], sessionOptions[iCharmPart], inputShapesML[iCharmPart], dataTypeML[iCharmPart], loadModelsFromCCDB, ccdbApi, mlModelPathCCDB.value, timestampCCDB));
        }
      }
    }
    // safety for optimisation tree
    if (applyOptimisation && !applyML) {
      LOG(fatal) << "Can't apply optimisation if ML is not applied.";
    }
  }

  using BigTracksMCPID = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::McTrackLabels>;
  using BigTracksPID = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr>;
  using CollsWithEvSel = soa::Join<aod::Collisions, aod::EvSels>;

  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;
  Preslice<aod::V0Datas> v0sPerCollision = aod::v0data::collisionId;
  Preslice<aod::Hf2Prongs> hf2ProngPerCollision = aod::track_association::collisionId;
  Preslice<aod::Hf3Prongs> hf3ProngPerCollision = aod::track_association::collisionId;
  Preslice<aod::CascDatas> cascPerCollision = aod::cascdata::collisionId;

  void process(CollsWithEvSel const& collisions,
               aod::BCsWithTimestamps const&,
               aod::V0Datas const& theV0s,
               aod::V0sLinked const& v0Links,
               aod::CascDatas const& cascades,
               aod::Hf2Prongs const& cand2Prongs,
               aod::Hf3Prongs const& cand3Prongs,
               aod::TrackAssoc const& trackIndices,
               BigTracksPID const& tracks)
  {
    for (const auto& collision : collisions) {

      bool keepEvent[kNtriggersHF]{false};
      if (applyEventSelection && (!collision.sel8() || std::fabs(collision.posZ()) > 11.f)) { // safety margin for Zvtx
        tags(keepEvent[kHighPt2P], keepEvent[kHighPt3P], keepEvent[kBeauty3P], keepEvent[kBeauty4P], keepEvent[kFemto2P], keepEvent[kFemto3P], keepEvent[kDoubleCharm2P], keepEvent[kDoubleCharm3P], keepEvent[kDoubleCharmMix], keepEvent[kV0Charm2P], keepEvent[kV0Charm3P], keepEvent[kCharmBarToXiBach]);
        continue;
      }

      auto thisCollId = collision.globalIndex();

      if (applyOptimisation) {
        optimisationTreeCollisions(thisCollId);
      }

      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();

      if (applyML && (loadModelsFromCCDB && timestampCCDB == 0) && !sessionML[kD0]) {
        for (auto iCharmPart{0}; iCharmPart < kNCharmParticles; ++iCharmPart) {
          if (onnxFiles[iCharmPart] != "") {
            sessionML[iCharmPart].reset(helper.initONNXSession(onnxFiles[iCharmPart], charmParticleNames[iCharmPart], envML[iCharmPart], sessionOptions[iCharmPart], inputShapesML[iCharmPart], dataTypeML[iCharmPart], loadModelsFromCCDB, ccdbApi, mlModelPathCCDB.value, bc.timestamp()));
          }
        }
      }

      // needed for track propagation
      if (currentRun != bc.runNumber()) {
        o2::parameters::GRPMagField* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", bc.timestamp());
        o2::base::Propagator::initFieldFromGRP(grpo);

        // needed for TPC PID postcalibrations
        if (setTPCCalib == 1) {
          helper.setTpcRecalibMaps(ccdb, bc, ccdbPathTPC);
        } else if (setTPCCalib > 1) {
          helper.setValuesBB(ccdbApi, bc, std::array{ccdbBBPion.value, ccdbBBAntiPion.value, ccdbBBKaon.value, ccdbBBAntiKaon.value, ccdbBBProton.value, ccdbBBAntiProton.value});
        }

        currentRun = bc.runNumber();
      }

      hProcessedEvents->Fill(0);

      std::vector<std::vector<int64_t>> indicesDau2Prong{};

      auto cand2ProngsThisColl = cand2Prongs.sliceBy(hf2ProngPerCollision, thisCollId);
      for (const auto& cand2Prong : cand2ProngsThisColl) {                                // start loop over 2 prongs
        if (!TESTBIT(cand2Prong.hfflag(), o2::aod::hf_cand_2prong::DecayType::D0ToPiK)) { // check if it's a D0
          continue;
        }

        auto trackPos = cand2Prong.prong0_as<BigTracksPID>(); // positive daughter
        auto trackNeg = cand2Prong.prong1_as<BigTracksPID>(); // negative daughter

        auto preselD0 = helper.isDzeroPreselected(trackPos, trackNeg);
        if (!preselD0) {
          continue;
        }

        auto trackParPos = getTrackPar(trackPos);
        auto trackParNeg = getTrackPar(trackNeg);
        o2::gpu::gpustd::array<float, 2> dcaPos{trackPos.dcaXY(), trackPos.dcaZ()};
        o2::gpu::gpustd::array<float, 2> dcaNeg{trackNeg.dcaXY(), trackNeg.dcaZ()};
        std::array<float, 3> pVecPos{trackPos.px(), trackPos.py(), trackPos.pz()};
        std::array<float, 3> pVecNeg{trackNeg.px(), trackNeg.py(), trackNeg.pz()};
        if (trackPos.collisionId() != thisCollId) {
          o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParPos, 2.f, noMatCorr, &dcaPos);
          getPxPyPz(trackParPos, pVecPos);
        }
        if (trackNeg.collisionId() != thisCollId) {
          o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParNeg, 2.f, noMatCorr, &dcaNeg);
          getPxPyPz(trackParNeg, pVecNeg);
        }

        bool isCharmTagged{true}, isBeautyTagged{true};

        // apply ML models
        int tagBDT = 0;
        float scoresToFill[3] = {-1., -1., -1.};
        if (applyML && onnxFiles[kD0] != "") {
          isCharmTagged = false;
          isBeautyTagged = false;

          // TODO: add more feature configurations
          std::vector<float> inputFeaturesD0{trackParPos.getPt(), dcaPos[0], dcaPos[1], trackParNeg.getPt(), dcaNeg[0], dcaNeg[1]};
          std::vector<double> inputFeaturesDoD0{trackParPos.getPt(), dcaPos[0], dcaPos[1], trackParNeg.getPt(), dcaNeg[0], dcaNeg[1]};

          if (dataTypeML[kD0] == 1) {
            auto scores = helper.predictONNX(inputFeaturesD0, sessionML[kD0], inputShapesML[kD0]);
            tagBDT = helper.isBDTSelected(scores, thresholdBDTScores[kD0]);
            for (int iScore{0}; iScore < 3; ++iScore) {
              scoresToFill[iScore] = scores[iScore];
            }
          } else if (dataTypeML[kD0] == 11) {
            auto scores = helper.predictONNX(inputFeaturesDoD0, sessionML[kD0], inputShapesML[kD0]);
            tagBDT = helper.isBDTSelected(scores, thresholdBDTScores[kD0]);
            for (int iScore{0}; iScore < 3; ++iScore) {
              scoresToFill[iScore] = scores[iScore];
            }
          } else {
            LOG(fatal) << "Error running model inference for D0: Unexpected input data type.";
          }

          if (applyML && activateQA > 1) {
            hBDTScoreBkg[kD0]->Fill(scoresToFill[0]);
            hBDTScorePrompt[kD0]->Fill(scoresToFill[1]);
            hBDTScoreNonPrompt[kD0]->Fill(scoresToFill[2]);
          }

          isCharmTagged = TESTBIT(tagBDT, RecoDecay::OriginType::Prompt);
          isBeautyTagged = TESTBIT(tagBDT, RecoDecay::OriginType::NonPrompt);
        }

        if (!isCharmTagged && !isBeautyTagged) {
          continue;
        }

        auto pVec2Prong = RecoDecay::pVec(pVecPos, pVecNeg);
        auto pt2Prong = RecoDecay::pt(pVec2Prong);

        if (applyOptimisation) {
          optimisationTreeCharm(thisCollId, pdg::Code::kD0, pt2Prong, scoresToFill[0], scoresToFill[1], scoresToFill[2]);
        }

        auto selD0 = helper.isSelectedD0InMassRange(pVecPos, pVecNeg, pt2Prong, preselD0, activateQA, hMassVsPtC[kD0]);

        if (pt2Prong >= ptThresholds->get(0u, 0u)) {
          keepEvent[kHighPt2P] = true;
          if (activateQA) {
            hCharmHighPt[kD0]->Fill(pt2Prong);
          }
        } // end high-pT selection

        if (isCharmTagged) {
          indicesDau2Prong.push_back(std::vector<int64_t>{trackPos.globalIndex(), trackNeg.globalIndex()});
        } // end multi-charm selection

        // compute masses already here, needed both for B0 --> D* (--> D0 Pi) Pi and Ds1 --> D* (--> D0 Pi) K0S
        auto massD0Cand = RecoDecay::m(std::array{pVecPos, pVecNeg}, std::array{massPi, massKa});
        auto massD0BarCand = RecoDecay::m(std::array{pVecPos, pVecNeg}, std::array{massKa, massPi});

        auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
        for (const auto& trackId : trackIdsThisCollision) { // start loop over tracks
          auto track = trackId.track_as<BigTracksPID>();

          if (track.globalIndex() == trackPos.globalIndex() || track.globalIndex() == trackNeg.globalIndex()) {
            continue;
          }

          auto trackParThird = getTrackPar(track);
          o2::gpu::gpustd::array<float, 2> dcaThird{track.dcaXY(), track.dcaZ()};
          std::array<float, 3> pVecThird = {track.px(), track.py(), track.pz()};
          if (track.collisionId() != thisCollId) {
            o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParThird, 2.f, noMatCorr, &dcaThird);
            getPxPyPz(trackParThird, pVecThird);
          }

          if (!keepEvent[kBeauty3P] && isBeautyTagged) {
            auto isTrackSelected = helper.isSelectedTrackForSoftPionOrBeauty(track, trackParThird, dcaThird, kBeauty3P);
            if (isTrackSelected && ((TESTBIT(selD0, 0) && track.sign() < 0) || (TESTBIT(selD0, 1) && track.sign() > 0))) {
              auto massCand = RecoDecay::m(std::array{pVec2Prong, pVecThird}, std::array{massD0, massPi});
              auto pVecBeauty3Prong = RecoDecay::pVec(pVec2Prong, pVecThird);
              auto ptCand = RecoDecay::pt(pVecBeauty3Prong);
              if (TESTBIT(isTrackSelected, kForBeauty) && std::fabs(massCand - massBPlus) <= deltaMassBeauty->get(0u, 0u)) {
                keepEvent[kBeauty3P] = true;
                // fill optimisation tree for D0
                if (applyOptimisation) {
                  optimisationTreeBeauty(thisCollId, pdg::Code::kD0, pt2Prong, scoresToFill[0], scoresToFill[1], scoresToFill[2], dcaThird[0]);
                }
                if (activateQA) {
                  hMassVsPtB[kBplus]->Fill(ptCand, massCand);
                }
              } else if (TESTBIT(isTrackSelected, kSoftPionForBeauty)) {
                std::array<float, 2> massDausD0{massPi, massKa};
                auto massD0dau = massD0Cand;
                if (track.sign() < 0) {
                  massDausD0[0] = massKa;
                  massDausD0[1] = massPi;
                  massD0dau = massD0BarCand;
                }
                auto massDstarCand = RecoDecay::m(std::array{pVecPos, pVecNeg, pVecThird}, std::array{massDausD0[0], massDausD0[1], massPi});
                auto massDiffDstar = massDstarCand - massD0dau;

                if (std::fabs(massDiffDstar - (massDStar - massD0)) <= maxDeltaMassCharmReso->get(0u, 0u)) { // additional check for B0->D*pi polarization studies
                  if (activateQA) {
                    hMassVsPtC[kNCharmParticles]->Fill(ptCand, massDiffDstar);
                  }
                  for (const auto& trackIdB : trackIdsThisCollision) { // start loop over tracks
                    auto trackB = trackIdB.track_as<BigTracksPID>();
                    if (track.globalIndex() == trackB.globalIndex()) {
                      continue;
                    }
                    auto trackParFourth = getTrackPar(trackB);
                    o2::gpu::gpustd::array<float, 2> dcaFourth{trackB.dcaXY(), trackB.dcaZ()};
                    std::array<float, 3> pVecFourth = {trackB.px(), trackB.py(), trackB.pz()};
                    if (trackB.collisionId() != thisCollId) {
                      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParFourth, 2.f, noMatCorr, &dcaFourth);
                      getPxPyPz(trackParFourth, pVecFourth);
                    }

                    auto isTrackFourthSelected = helper.isSelectedTrackForSoftPionOrBeauty(trackB, trackParFourth, dcaFourth, kBeauty3P);
                    if (track.sign() * trackB.sign() < 0 && TESTBIT(isTrackFourthSelected, kForBeauty)) {
                      auto massCandB0 = RecoDecay::m(std::array{pVecBeauty3Prong, pVecFourth}, std::array{massDStar, massPi});
                      if (std::fabs(massCandB0 - massB0) <= deltaMassBeauty->get(0u, 2u)) {
                        keepEvent[kBeauty3P] = true;
                        // fill optimisation tree for D0
                        if (applyOptimisation) {
                          optimisationTreeBeauty(thisCollId, 413, pt2Prong, scoresToFill[0], scoresToFill[1], scoresToFill[2], dcaFourth[0]); // pdgCode of D*(2010)+: 413
                        }
                        if (activateQA) {
                          auto pVecBeauty4Prong = RecoDecay::pVec(pVec2Prong, pVecThird, pVecFourth);
                          auto ptCandBeauty4Prong = RecoDecay::pt(pVecBeauty4Prong);
                          hMassVsPtB[kB0toDStar]->Fill(ptCandBeauty4Prong, massCandB0);
                        }
                      }
                    }
                  }
                }
              }
            }
          } // end beauty selection

          // 2-prong femto
          if (!keepEvent[kFemto2P] && enableFemtoChannels->get(0u, 0u) && isCharmTagged && track.collisionId() == thisCollId && (TESTBIT(selD0, 0) || TESTBIT(selD0, 1) || !requireCharmMassForFemto)) {
            bool isProton = helper.isSelectedProton4Femto(track, trackParThird, activateQA, hProtonTPCPID, hProtonTOFPID);
            if (isProton) {
              float relativeMomentum = helper.computeRelativeMomentum(pVecThird, pVec2Prong, massD0);
              if (applyOptimisation) {
                optimisationTreeFemto(thisCollId, pdg::Code::kD0, pt2Prong, scoresToFill[0], scoresToFill[1], scoresToFill[2], relativeMomentum, track.tpcNSigmaPr(), track.tofNSigmaPr());
              }
              if (relativeMomentum < femtoMaxRelativeMomentum) {
                keepEvent[kFemto2P] = true;
                if (activateQA) {
                  hCharmProtonKstarDistr[kD0]->Fill(relativeMomentum);
                }
              }
            }
          } // end femto selection

        } // end loop over tracks

        // 2-prong with Gamma (conversion photon)
        auto v0sThisCollision = theV0s.sliceBy(v0sPerCollision, thisCollId);
        for (const auto& v0 : v0sThisCollision) {
          if (!keepEvent[kV0Charm2P] && (isCharmTagged || isBeautyTagged) && (TESTBIT(selD0, 0) || TESTBIT(selD0, 1))) {
            auto posTrack = v0.posTrack_as<BigTracksPID>();
            auto negTrack = v0.negTrack_as<BigTracksPID>();
            auto selV0 = helper.isSelectedV0(v0, std::array{posTrack, negTrack}, collision, activateQA, hV0Selected, hArmPod);
            if (selV0) {
              // propagate to PV
              gpu::gpustd::array<float, 2> dcaInfo;
              std::array<float, 3> pVecV0 = {v0.px(), v0.py(), v0.pz()};
              auto trackParV0 = o2::track::TrackPar(std::array{v0.x(), v0.y(), v0.z()}, pVecV0, 0, true);
              o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParV0, 2.f, matCorr, &dcaInfo);
              getPxPyPz(trackParV0, pVecV0);
              if (TESTBIT(selV0, kPhoton)) {
                float massDStarCand{-1.}, massDStarBarCand{999.};
                float massDiffDstar{-1.}, massDiffDstarBar{999.};
                if (TESTBIT(selD0, 0)) {
                  massDStarCand = RecoDecay::m(std::array{pVecPos, pVecNeg, pVecV0}, std::array{massPi, massKa, massGamma});
                  massDiffDstar = massDStarCand - massD0Cand;
                }
                if (TESTBIT(selD0, 1)) {
                  massDStarBarCand = RecoDecay::m(std::array{pVecPos, pVecNeg, pVecV0}, std::array{massKa, massPi, massGamma});
                  massDiffDstarBar = massDStarBarCand - massD0BarCand;
                }
                bool isGoodDstar = (massDiffDstar < maxDeltaMassCharmReso->get(0u, 1u));
                bool isGoodDstarBar = (massDiffDstarBar < maxDeltaMassCharmReso->get(0u, 1u));

                if (isGoodDstar || isGoodDstarBar) {
                  keepEvent[kV0Charm2P] = true;
                  if (activateQA) {
                    auto pVecReso2Prong = RecoDecay::pVec(pVec2Prong, pVecV0);
                    auto ptCand = RecoDecay::pt(pVecReso2Prong);
                    if (isGoodDstar) {
                      hMassVsPtC[kNCharmParticles + 1]->Fill(ptCand, massDiffDstar);
                    }
                    if (isGoodDstarBar) {
                      hMassVsPtC[kNCharmParticles + 1]->Fill(ptCand, massDiffDstarBar);
                    }
                  }
                }
              }
              if (!keepEvent[kV0Charm2P] && TESTBIT(selV0, kK0S)) {

                // we first look for a D*+
                for (const auto& trackBachelorId : trackIdsThisCollision) { // start loop over tracks
                  auto trackBachelor = trackBachelorId.track_as<BigTracksPID>();
                  if (trackBachelor.globalIndex() == trackPos.globalIndex() || trackBachelor.globalIndex() == trackNeg.globalIndex()) {
                    continue;
                  }

                  auto trackParBachelor = getTrackPar(trackBachelor);
                  o2::gpu::gpustd::array<float, 2> dcaBachelor{trackBachelor.dcaXY(), trackBachelor.dcaZ()};
                  std::array<float, 3> pVecBachelor = {trackBachelor.px(), trackBachelor.py(), trackBachelor.pz()};
                  if (trackBachelor.collisionId() != thisCollId) {
                    o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParBachelor, 2.f, noMatCorr, &dcaBachelor);
                    getPxPyPz(trackParBachelor, pVecBachelor);
                  }

                  int isTrackSelected = helper.isSelectedTrackForSoftPionOrBeauty(trackBachelor, trackParBachelor, dcaBachelor, -1);
                  if (TESTBIT(isTrackSelected, kSoftPion) && ((TESTBIT(selD0, 0) && trackBachelor.sign() < 0) || (TESTBIT(selD0, 1) && trackBachelor.sign() > 0))) {
                    std::array<float, 2> massDausD0{massPi, massKa};
                    auto massD0dau = massD0Cand;
                    if (trackBachelor.sign() < 0) {
                      massDausD0[0] = massKa;
                      massDausD0[1] = massPi;
                      massD0dau = massD0BarCand;
                    }
                    auto massDStarCand = RecoDecay::m(std::array{pVecPos, pVecNeg, pVecBachelor}, std::array{massDausD0[0], massDausD0[1], massPi});
                    auto massDiffDstar = massDStarCand - massD0dau;
                    auto pVecDStarCand = RecoDecay::pVec(pVec2Prong, pVecBachelor);
                    auto ptDStarCand = RecoDecay::pt(pVecDStarCand);
                    if (std::fabs(massDiffDstar - (massDStar - massD0)) <= maxDeltaMassCharmReso->get(0u, 0u)) {
                      if (activateQA) {
                        hMassVsPtC[kNCharmParticles]->Fill(ptDStarCand, massDiffDstar);
                      }
                      auto massDStarK0S = RecoDecay::m(std::array{pVecPos, pVecNeg, pVecBachelor, pVecV0}, std::array{massDausD0[0], massDausD0[1], massPi, massK0S});
                      auto massDiffDsReso = massDStarK0S - massDStarCand;
                      if (massDiffDsReso < maxDeltaMassCharmReso->get(0u, 3u)) {
                        if (activateQA) {
                          auto pVecReso2Prong = RecoDecay::pVec(pVecDStarCand, pVecV0);
                          auto ptCand = RecoDecay::pt(pVecReso2Prong);
                          hMassVsPtC[kNCharmParticles + 3]->Fill(ptCand, massDiffDsReso);
                        }
                        keepEvent[kV0Charm2P] = true;
                      }
                    }
                  }
                }
              }
            }
          }
        } // end V0 selection

      } // end loop over 2-prong candidates

      std::vector<std::vector<int64_t>> indicesDau3Prong{};
      auto cand3ProngsThisColl = cand3Prongs.sliceBy(hf3ProngPerCollision, thisCollId);
      for (const auto& cand3Prong : cand3ProngsThisColl) { // start loop over 3 prongs
        std::array<int8_t, kNCharmParticles - 1> is3Prong = {
          TESTBIT(cand3Prong.hfflag(), o2::aod::hf_cand_3prong::DecayType::DplusToPiKPi),
          TESTBIT(cand3Prong.hfflag(), o2::aod::hf_cand_3prong::DecayType::DsToKKPi),
          TESTBIT(cand3Prong.hfflag(), o2::aod::hf_cand_3prong::DecayType::LcToPKPi),
          TESTBIT(cand3Prong.hfflag(), o2::aod::hf_cand_3prong::DecayType::XicToPKPi)};
        if (!std::accumulate(is3Prong.begin(), is3Prong.end(), 0)) { // check if it's a D+, Ds+, Lc+ or Xic+
          continue;
        }

        auto trackFirst = cand3Prong.prong0_as<BigTracksPID>();
        auto trackSecond = cand3Prong.prong1_as<BigTracksPID>();
        auto trackThird = cand3Prong.prong2_as<BigTracksPID>();

        auto trackParFirst = getTrackPar(trackFirst);
        auto trackParSecond = getTrackPar(trackSecond);
        auto trackParThird = getTrackPar(trackThird);
        o2::gpu::gpustd::array<float, 2> dcaFirst{trackFirst.dcaXY(), trackFirst.dcaZ()};
        o2::gpu::gpustd::array<float, 2> dcaSecond{trackSecond.dcaXY(), trackSecond.dcaZ()};
        o2::gpu::gpustd::array<float, 2> dcaThird{trackThird.dcaXY(), trackThird.dcaZ()};
        std::array<float, 3> pVecFirst = {trackFirst.px(), trackFirst.py(), trackFirst.pz()};
        std::array<float, 3> pVecSecond = {trackSecond.px(), trackSecond.py(), trackSecond.pz()};
        std::array<float, 3> pVecThird = {trackThird.px(), trackThird.py(), trackThird.pz()};
        if (trackFirst.collisionId() != thisCollId) {
          o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParFirst, 2.f, noMatCorr, &dcaFirst);
          getPxPyPz(trackParFirst, pVecFirst);
        }
        if (trackSecond.collisionId() != thisCollId) {
          o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParSecond, 2.f, noMatCorr, &dcaSecond);
          getPxPyPz(trackParSecond, pVecSecond);
        }
        if (trackThird.collisionId() != thisCollId) {
          o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParThird, 2.f, noMatCorr, &dcaThird);
          getPxPyPz(trackParThird, pVecThird);
        }

        if (is3Prong[0]) { // D+ preselections
          is3Prong[0] = helper.isDplusPreselected(trackSecond);
        }
        if (is3Prong[1]) { // Ds preselections
          is3Prong[1] = helper.isDsPreselected(pVecFirst, pVecThird, pVecSecond, trackSecond);
        }
        if (is3Prong[2] || is3Prong[3]) { // charm baryon preselections
          auto presel = helper.isCharmBaryonPreselected(trackFirst, trackThird, trackSecond);
          if (is3Prong[2]) {
            is3Prong[2] = presel;
          }
          if (is3Prong[3]) {
            is3Prong[3] = presel;
          }
        }

        std::array<int8_t, kNCharmParticles - 1> isCharmTagged = is3Prong;
        std::array<int8_t, kNCharmParticles - 1> isBeautyTagged = is3Prong;

        float scoresToFill[kNCharmParticles - 1][3];
        for (int i = 0; i < kNCharmParticles - 1; i++) {
          std::fill_n(scoresToFill[i], 3, -1);
        } // initialize BDT scores array outside ML loop
        // apply ML models
        if (applyML) {
          isCharmTagged = std::array<int8_t, kNCharmParticles - 1>{0};
          isBeautyTagged = std::array<int8_t, kNCharmParticles - 1>{0};

          // TODO: add more feature configurations
          std::vector<float> inputFeatures{trackParFirst.getPt(), dcaFirst[0], dcaFirst[1], trackParSecond.getPt(), dcaSecond[0], dcaSecond[1], trackParThird.getPt(), dcaThird[0], dcaThird[1]};
          std::vector<double> inputFeaturesD{trackParFirst.getPt(), dcaFirst[0], dcaFirst[1], trackParSecond.getPt(), dcaSecond[0], dcaSecond[1], trackParThird.getPt(), dcaThird[0], dcaThird[1]};
          for (auto iCharmPart{0}; iCharmPart < kNCharmParticles - 1; ++iCharmPart) {
            if (!is3Prong[iCharmPart] || onnxFiles[iCharmPart + 1] == "") {
              continue;
            }

            int tagBDT = 0;
            if (dataTypeML[iCharmPart + 1] == 1) {
              auto scores = helper.predictONNX(inputFeatures, sessionML[iCharmPart + 1], inputShapesML[iCharmPart + 1]);
              tagBDT = helper.isBDTSelected(scores, thresholdBDTScores[iCharmPart + 1]);
              for (int iScore{0}; iScore < 3; ++iScore) {
                scoresToFill[iCharmPart][iScore] = scores[iScore];
              }
            } else if (dataTypeML[iCharmPart + 1] == 11) {
              auto scores = helper.predictONNX(inputFeaturesD, sessionML[iCharmPart + 1], inputShapesML[iCharmPart + 1]);
              tagBDT = helper.isBDTSelected(scores, thresholdBDTScores[iCharmPart + 1]);
              for (int iScore{0}; iScore < 3; ++iScore) {
                scoresToFill[iCharmPart][iScore] = scores[iScore];
              }
            } else {
              LOG(error) << "Error running model inference for " << charmParticleNames[iCharmPart + 1].data() << ": Unexpected input data type.";
            }

            isCharmTagged[iCharmPart] = TESTBIT(tagBDT, RecoDecay::OriginType::Prompt);
            isBeautyTagged[iCharmPart] = TESTBIT(tagBDT, RecoDecay::OriginType::NonPrompt);

            if (activateQA > 1) {
              hBDTScoreBkg[iCharmPart + 1]->Fill(scoresToFill[iCharmPart][0]);
              hBDTScorePrompt[iCharmPart + 1]->Fill(scoresToFill[iCharmPart][1]);
              hBDTScoreNonPrompt[iCharmPart + 1]->Fill(scoresToFill[iCharmPart][2]);
            }
          }
        }

        if (!std::accumulate(isCharmTagged.begin(), isCharmTagged.end(), 0) && !std::accumulate(isBeautyTagged.begin(), isBeautyTagged.end(), 0)) {
          continue;
        }

        if (std::accumulate(isCharmTagged.begin(), isCharmTagged.end(), 0)) {
          indicesDau3Prong.push_back(std::vector<int64_t>{trackFirst.globalIndex(), trackSecond.globalIndex(), trackThird.globalIndex()});
        } // end multiple 3-prong selection

        auto pVec3Prong = RecoDecay::pVec(pVecFirst, pVecSecond, pVecThird);
        auto pt3Prong = RecoDecay::pt(pVec3Prong);
        float sign3Prong = -1 * trackFirst.sign() * trackSecond.sign() * trackThird.sign();

        std::array<int8_t, kNCharmParticles - 1> is3ProngInMass{0};
        if (is3Prong[0]) {
          is3ProngInMass[0] = helper.isSelectedDplusInMassRange(pVecFirst, pVecThird, pVecSecond, pt3Prong, activateQA, hMassVsPtC[kDplus]);
          if (applyOptimisation) {
            optimisationTreeCharm(thisCollId, pdg::Code::kDPlus, pt3Prong, scoresToFill[0][0], scoresToFill[0][1], scoresToFill[0][2]);
          }
        }
        if (is3Prong[1]) {
          is3ProngInMass[1] = helper.isSelectedDsInMassRange(pVecFirst, pVecThird, pVecSecond, pt3Prong, is3Prong[1], activateQA, hMassVsPtC[kDs]);
          if (applyOptimisation) {
            optimisationTreeCharm(thisCollId, pdg::Code::kDS, pt3Prong, scoresToFill[1][0], scoresToFill[1][1], scoresToFill[1][2]);
          }
        }
        if (is3Prong[2]) {
          is3ProngInMass[2] = helper.isSelectedLcInMassRange(pVecFirst, pVecThird, pVecSecond, pt3Prong, is3Prong[2], activateQA, hMassVsPtC[kLc]);
          if (applyOptimisation) {
            optimisationTreeCharm(thisCollId, pdg::Code::kLambdaCPlus, pt3Prong, scoresToFill[2][0], scoresToFill[2][1], scoresToFill[2][2]);
          }
        }
        if (is3Prong[3]) {
          is3ProngInMass[3] = helper.isSelectedXicInMassRange(pVecFirst, pVecThird, pVecSecond, pt3Prong, is3Prong[3], activateQA, hMassVsPtC[kXic]);
          if (applyOptimisation) {
            optimisationTreeCharm(thisCollId, pdg::Code::kXiCPlus, pt3Prong, scoresToFill[3][0], scoresToFill[3][1], scoresToFill[3][2]);
          }
        }

        if (pt3Prong >= ptThresholds->get(0u, 1u)) {
          keepEvent[kHighPt3P] = true;
          if (activateQA) {
            for (auto iCharmPart{1}; iCharmPart < kNCharmParticles; ++iCharmPart) {
              if (is3Prong[iCharmPart - 1] && (isCharmTagged[iCharmPart - 1] || isBeautyTagged[iCharmPart - 1])) {
                hCharmHighPt[iCharmPart]->Fill(pt3Prong);
              }
            }
          }
        } // end high-pT selection

        auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);

        for (const auto& trackId : trackIdsThisCollision) { // start loop over track indices as associated to this collision in HF code
          auto track = trackId.track_as<BigTracksPID>();
          if (track.globalIndex() == trackFirst.globalIndex() || track.globalIndex() == trackSecond.globalIndex() || track.globalIndex() == trackThird.globalIndex()) {
            continue;
          }

          auto trackParFourth = getTrackPar(track);
          o2::gpu::gpustd::array<float, 2> dcaFourth{track.dcaXY(), track.dcaZ()};
          std::array<float, 3> pVecFourth = {track.px(), track.py(), track.pz()};
          if (track.collisionId() != thisCollId) {
            o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParFourth, 2.f, noMatCorr, &dcaFourth);
            getPxPyPz(trackParFourth, pVecFourth);
          }

          int charmParticleID[kNBeautyParticles - 2] = {pdg::Code::kDPlus, pdg::Code::kDS, pdg::Code::kLambdaCPlus, pdg::Code::kXiCPlus};

          float massCharmHypos[kNBeautyParticles - 2] = {massDPlus, massDs, massLc, massXic};
          float massBeautyHypos[kNBeautyParticles - 2] = {massB0, massBs, massLb, massXib};
          float deltaMassHypos[kNBeautyParticles - 2] = {deltaMassBeauty->get(0u, 1u), deltaMassBeauty->get(0u, 3u), deltaMassBeauty->get(0u, 4u), deltaMassBeauty->get(0u, 5u)};
          auto isTrackSelected = helper.isSelectedTrackForSoftPionOrBeauty(track, trackParFourth, dcaFourth, kBeauty4P);
          if (track.sign() * sign3Prong < 0 && TESTBIT(isTrackSelected, kForBeauty)) {
            for (int iHypo{0}; iHypo < kNBeautyParticles - 2 && !keepEvent[kBeauty4P]; ++iHypo) {
              if (isBeautyTagged[iHypo] && (TESTBIT(is3ProngInMass[iHypo], 0) || TESTBIT(is3ProngInMass[iHypo], 1))) {
                auto massCandB = RecoDecay::m(std::array{pVec3Prong, pVecFourth}, std::array{massCharmHypos[iHypo], massPi});
                if (std::fabs(massCandB - massBeautyHypos[iHypo]) <= deltaMassHypos[iHypo]) {
                  keepEvent[kBeauty4P] = true;
                  if (applyOptimisation) {
                    optimisationTreeBeauty(thisCollId, charmParticleID[iHypo], pt3Prong, scoresToFill[iHypo][0], scoresToFill[iHypo][1], scoresToFill[iHypo][2], dcaFourth[0]);
                  }
                  if (activateQA) {
                    auto pVecBeauty4Prong = RecoDecay::pVec(pVec3Prong, pVecFourth);
                    auto ptCandBeauty4Prong = RecoDecay::pt(pVecBeauty4Prong);
                    hMassVsPtB[iHypo + 2]->Fill(ptCandBeauty4Prong, massCandB);
                  }
                }
              }
            }
          } // end beauty selection

          // 3-prong femto
          bool isProton = helper.isSelectedProton4Femto(track, trackParFourth, activateQA, hProtonTPCPID, hProtonTOFPID);
          if (isProton && track.collisionId() == thisCollId) {
            for (int iHypo{0}; iHypo < kNCharmParticles - 1 && !keepEvent[kFemto3P]; ++iHypo) {
              if (isCharmTagged[iHypo] && enableFemtoChannels->get(0u, iHypo + 1) && (TESTBIT(is3ProngInMass[iHypo], 0) || TESTBIT(is3ProngInMass[iHypo], 1) || !requireCharmMassForFemto)) {
                float relativeMomentum = helper.computeRelativeMomentum(pVecFourth, pVec3Prong, massCharmHypos[iHypo]);
                if (applyOptimisation) {
                  optimisationTreeFemto(thisCollId, charmParticleID[iHypo], pt3Prong, scoresToFill[iHypo][0], scoresToFill[iHypo][1], scoresToFill[iHypo][2], relativeMomentum, track.tpcNSigmaPr(), track.tofNSigmaPr());
                }
                if (relativeMomentum < femtoMaxRelativeMomentum) {
                  keepEvent[kFemto3P] = true;
                  if (activateQA) {
                    hCharmProtonKstarDistr[iHypo + 1]->Fill(relativeMomentum);
                  }
                }
              }
            }
          } // end femto selection

        } // end loop over tracks

        // 3-prong with V0 (Ds gamma, D+ K0S, D+ Lambda)
        auto v0sThisCollision = theV0s.sliceBy(v0sPerCollision, thisCollId);
        bool isGoodDsToKKPi = (isCharmTagged[kDs - 1] || isBeautyTagged[kDs - 1]) && TESTBIT(is3ProngInMass[kDs - 1], 0);
        bool isGoodDsToPiKK = (isCharmTagged[kDs - 1] || isBeautyTagged[kDs - 1]) && TESTBIT(is3ProngInMass[kDs - 1], 1);
        bool isGoodDPlus = (isCharmTagged[kDplus - 1] || isBeautyTagged[kDplus - 1]) && is3ProngInMass[kDplus - 1];
        auto massDPlusCand = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massPi, massKa, massPi});
        auto massDsKKPi = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massKa, massKa, massPi});
        auto massDsPiKK = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massPi, massKa, massKa});
        for (const auto& v0 : v0sThisCollision) {
          if (!keepEvent[kV0Charm3P] && (isGoodDsToKKPi || isGoodDsToPiKK || isGoodDPlus)) {
            auto posTrack = v0.posTrack_as<BigTracksPID>();
            auto negTrack = v0.negTrack_as<BigTracksPID>();
            auto selV0 = helper.isSelectedV0(v0, std::array{posTrack, negTrack}, collision, activateQA, hV0Selected, hArmPod);
            if (selV0 > 0) {
              // propagate to PV
              gpu::gpustd::array<float, 2> dcaInfo;
              std::array<float, 3> pVecV0 = {v0.px(), v0.py(), v0.pz()};
              auto trackParV0 = o2::track::TrackPar(std::array{v0.x(), v0.y(), v0.z()}, pVecV0, 0, true);
              o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParV0, 2.f, matCorr, &dcaInfo);
              getPxPyPz(trackParV0, pVecV0);

              if (!keepEvent[kV0Charm3P] && (isGoodDsToKKPi || isGoodDsToPiKK) && TESTBIT(selV0, kPhoton)) {
                float massDsStarToKKPiCand{-1.}, massDsStarToPiKKCand{999.};
                float massDiffDsStarToKKPi{-1.}, massDiffDsStarToPiKK{999.};
                if (isGoodDsToKKPi) {
                  massDsStarToKKPiCand = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird, pVecV0}, std::array{massKa, massKa, massPi, massGamma});
                  massDiffDsStarToKKPi = massDsStarToKKPiCand - massDsKKPi;
                }
                if (isGoodDsToPiKK) {
                  massDsStarToPiKKCand = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird, pVecV0}, std::array{massPi, massKa, massKa, massGamma});
                  massDiffDsStarToPiKK = massDsStarToPiKKCand - massDsPiKK;
                }
                bool isGoodDsStarToKKPi = (massDiffDsStarToKKPi < maxDeltaMassCharmReso->get(0u, 1u));
                bool isGoodDsStarToPiKK = (massDiffDsStarToPiKK < maxDeltaMassCharmReso->get(0u, 1u));

                if (isGoodDsStarToKKPi || isGoodDsStarToPiKK) {
                  keepEvent[kV0Charm3P] = true;
                  if (activateQA) {
                    auto pVecReso3Prong = RecoDecay::pVec(pVec3Prong, pVecV0);
                    auto ptCand = RecoDecay::pt(pVecReso3Prong);
                    if (isGoodDsStarToKKPi) {
                      hMassVsPtC[kNCharmParticles + 2]->Fill(ptCand, massDiffDsStarToKKPi);
                    }
                    if (isGoodDsStarToPiKK) {
                      hMassVsPtC[kNCharmParticles + 2]->Fill(ptCand, massDiffDsStarToPiKK);
                    }
                  }
                }
              }
              if (!keepEvent[kV0Charm3P] && isGoodDPlus) {
                if (TESTBIT(selV0, kK0S)) { // Ds2*
                  auto massDsStarCand = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird, pVecV0}, std::array{massPi, massKa, massPi, massK0S});
                  auto massDiffDsStar = massDsStarCand - massDPlusCand;
                  if (massDiffDsStar < maxDeltaMassCharmReso->get(0u, 4u)) {
                    if (activateQA) {
                      auto pVecReso3Prong = RecoDecay::pVec(pVec3Prong, pVecV0);
                      auto ptCand = RecoDecay::pt(pVecReso3Prong);
                      hMassVsPtC[kNCharmParticles + 4]->Fill(ptCand, massDiffDsStar);
                    }
                    keepEvent[kV0Charm3P] = true;
                  }
                }
                if ((TESTBIT(selV0, kLambda) && sign3Prong > 0) || (TESTBIT(selV0, kAntiLambda) && sign3Prong < 0)) { // Xic(3055) and Xic(3080)
                  auto massXicStarCand = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird, pVecV0}, std::array{massPi, massKa, massPi, massLambda});
                  auto massDiffXicStar = massXicStarCand - massDPlusCand;
                  if (massDiffXicStar < maxDeltaMassCharmReso->get(0u, 5u)) {
                    if (activateQA) {
                      auto pVecReso3Prong = RecoDecay::pVec(pVec3Prong, pVecV0);
                      auto ptCand = RecoDecay::pt(pVecReso3Prong);
                      hMassVsPtC[kNCharmParticles + 5]->Fill(ptCand, massDiffXicStar);
                    }
                    keepEvent[kV0Charm3P] = true;
                  }
                }
              }
            }
          }
        } // end gamma selection

      } // end loop over 3-prong candidates

      if (!keepEvent[kCharmBarToXiBach]) {
        auto cascThisColl = cascades.sliceBy(cascPerCollision, thisCollId);
        for (const auto& casc : cascThisColl) {
          if (!casc.v0_as<aod::V0sLinked>().has_v0Data()) { // check that V0 data are stored
            continue;
          }
          auto bachelorCasc = casc.bachelor_as<BigTracksPID>();
          auto v0 = casc.v0_as<aod::V0sLinked>();
          auto v0Element = v0.v0Data_as<aod::V0Datas>();
          auto v0DauPos = v0Element.posTrack_as<BigTracksPID>();
          auto v0DauNeg = v0Element.negTrack_as<BigTracksPID>();
          if (!helper.isSelectedCascade(casc, v0Element, std::array{bachelorCasc, v0DauPos, v0DauNeg}, collision)) {
            continue;
          }
          if (activateQA) {
            hMassXi->Fill(casc.mXi());
          }

          auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
          for (const auto& trackId : trackIdsThisCollision) { // start loop over tracks
            auto track = trackId.track_as<BigTracksPID>();

            // ask for opposite sign daughters (omegac daughters)
            if (track.sign() * bachelorCasc.sign() >= 0) {
              continue;
            }

            // check if track is one of the Xi daughters
            if (track.globalIndex() == bachelorCasc.globalIndex() || track.globalIndex() == v0DauPos.globalIndex() || track.globalIndex() == v0DauNeg.globalIndex()) {
              continue;
            }

            // propagate to PV
            gpu::gpustd::array<float, 2> dcaInfo;
            std::array<float, 3> pVecCascade = {casc.px(), casc.py(), casc.pz()};
            auto trackParCasc = o2::track::TrackPar(std::array{casc.x(), casc.y(), casc.z()}, pVecCascade, bachelorCasc.sign(), true);
            o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParCasc, 2.f, matCorr, &dcaInfo);
            getPxPyPz(trackParCasc, pVecCascade);

            auto trackParBachelor = getTrackPar(track);
            std::array<float, 3> pVecBachelor = {track.px(), track.py(), track.pz()};
            if (track.collisionId() != thisCollId) {
              o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParBachelor, 2.f, noMatCorr, &dcaInfo);
              getPxPyPz(trackParBachelor, pVecBachelor);
            }

            auto isSelBachelor = helper.isSelectedBachelorForCharmBaryon(track, dcaInfo);
            if (isSelBachelor == kRejected) {
              continue;
            }

            auto ptCharmBaryon = RecoDecay::pt(RecoDecay::pVec(pVecCascade, pVecBachelor));

            if (!keepEvent[kCharmBarToXiBach] && TESTBIT(isSelBachelor, kPionForCharmBaryon)) {
              auto massXiPi = RecoDecay::m(std::array{pVecCascade, pVecBachelor}, std::array{massXi, massPi});
              if (ptCharmBaryon > cutsXiBachelor->get(0u, 0u) && massXiPi >= cutsXiBachelor->get(0u, 2u) && massXiPi <= 2.8f) {
                keepEvent[kCharmBarToXiBach] = true;
                if (activateQA) {
                  hMassVsPtC[kNCharmParticles + 6]->Fill(ptCharmBaryon, massXiPi);
                }
              }
            }
            if (!keepEvent[kCharmBarToXiBach] && TESTBIT(isSelBachelor, kKaonForCharmBaryon)) {
              auto massXiKa = RecoDecay::m(std::array{pVecCascade, pVecBachelor}, std::array{massXi, massKa});
              if (ptCharmBaryon > cutsXiBachelor->get(0u, 1u) && massXiKa >= cutsXiBachelor->get(0u, 3u) && massXiKa <= 2.8f) {
                keepEvent[kCharmBarToXiBach] = true;
                if (activateQA) {
                  hMassVsPtC[kNCharmParticles + 7]->Fill(ptCharmBaryon, massXiKa);
                }
              }
            }
          }
        }
      }

      auto n2Prongs = helper.computeNumberOfCandidates(indicesDau2Prong);
      auto n3Prongs = helper.computeNumberOfCandidates(indicesDau3Prong);
      indicesDau2Prong.insert(indicesDau2Prong.end(), indicesDau3Prong.begin(), indicesDau3Prong.end());
      auto n23Prongs = helper.computeNumberOfCandidates(indicesDau2Prong);

      if (activateQA) {
        hN2ProngCharmCand->Fill(n2Prongs);
        hN3ProngCharmCand->Fill(n3Prongs);
      }

      if (n2Prongs > 1 && enableDoubleCharmChannels->get(0u, 0u)) {
        keepEvent[kDoubleCharm2P] = true;
      }
      if (n3Prongs > 1 && enableDoubleCharmChannels->get(0u, 1u)) {
        keepEvent[kDoubleCharm3P] = true;
      }
      if (n23Prongs > 1 && enableDoubleCharmChannels->get(0u, 2u)) {
        keepEvent[kDoubleCharmMix] = true;
      }

      tags(keepEvent[kHighPt2P], keepEvent[kHighPt3P], keepEvent[kBeauty3P], keepEvent[kBeauty4P], keepEvent[kFemto2P], keepEvent[kFemto3P], keepEvent[kDoubleCharm2P], keepEvent[kDoubleCharm3P], keepEvent[kDoubleCharmMix], keepEvent[kV0Charm2P], keepEvent[kV0Charm3P], keepEvent[kCharmBarToXiBach]);

      if (!std::accumulate(keepEvent, keepEvent + kNtriggersHF, 0)) {
        hProcessedEvents->Fill(1);
      } else {
        for (int iTrigger{0}; iTrigger < kNtriggersHF; ++iTrigger) {
          if (keepEvent[iTrigger]) {
            hProcessedEvents->Fill(iTrigger + 2);
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{

  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<HfFilter>(cfg));

  return workflow;
}
