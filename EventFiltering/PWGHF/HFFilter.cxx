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

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include <onnxruntime/core/session/experimental_onnxruntime_cxx_api.h>

#include "EventFiltering/filterTables.h"
#include "HFFilterHelpers.h"

#include "Common/DataModel/CollisionAssociation.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hffilters;
using namespace hf_cuts_single_track;
using namespace hf_cuts_bdt_multiclass;

struct HfFilter { // Main struct for HF triggers

  Produces<aod::HfFilters> tags;
  Produces<aod::HFOptimisationTreeBeauty> optimisationTreeBeauty;
  Produces<aod::HFOptimisationTreeCharm> optimisationTreeCharm;
  Produces<aod::HFOptimisationTreeFemto> optimisationTreeFemto;
  Produces<aod::HFOptimisationTreeCollisions> optimisationTreeCollisions;

  Configurable<int> activateQA{"activateQA", 0, "flag to enable QA histos (0 no QA, 1 basic QA, 2 extended QA, 3 very extended QA)"};

  // parameters for all triggers
  Configurable<float> nsigmaTPCProtonLc{"nsigmaTPCProtonLc", 3., "Maximum value for TPC PID proton Nsigma for Lc candidates"};
  Configurable<float> nsigmaTOFProtonLc{"nsigmaTOFProtonLc", 3., "Maximum value for TOF PID proton Nsigma for Lc candidates"};
  Configurable<float> nsigmaTPCPionKaonDzero{"nsigmaTPCPionKaonDzero", 4., "Maximum value for TPC PID pion/kaon Nsigma for D0 candidates"};
  Configurable<float> nsigmaTOFPionKaonDzero{"nsigmaTOFPionKaonDzero", 4., "Maximum value for TOF PID pion/kaon Nsigma for D0 candidates"};
  Configurable<float> nsigmaTPCKaon3Prong{"nsigmaTPCKaon3Prong", 4., "Maximum value for TPC PID kaon Nsigma for all 3-prongs candidates"};
  Configurable<float> nsigmaTOFKaon3Prong{"nsigmaTOFKaon3Prong", 4., "Maximum value for TOF PID kaon Nsigma for all 3-prongs candidates"};

  // parameters for high-pT triggers
  Configurable<float> pTThreshold2Prong{"pTThreshold2Prong", 8., "pT treshold for high pT 2-prong candidates for kHighPt triggers in GeV/c"};
  Configurable<float> pTThreshold3Prong{"pTThreshold3Prong", 8., "pT treshold for high pT 3-prong candidates for kHighPt triggers in GeV/c"};

  // parameters for beauty triggers
  Configurable<float> deltaMassBPlus{"deltaMassBPlus", 0.3, "invariant-mass delta with respect to the B+ mass"};
  Configurable<float> deltaMassB0{"deltaMassB0", 0.3, "invariant-mass delta with respect to the B0 mass"};
  Configurable<float> deltaMassBs{"deltaMassBs", 0.3, "invariant-mass delta with respect to the Bs mass"};
  Configurable<float> deltaMassCharmHadronForBeauty{"deltaMassCharmHadronForBeauty", 0.04, "invariant-mass delta for charm"};
  Configurable<float> deltaMassLb{"deltaMassLb", 0.3, "invariant-mass delta with respect to the Lb mass"};
  Configurable<float> deltaMassXib{"deltaMassXib", 0.3, "invariant-mass delta with respect to the Lb mass"};
  Configurable<float> deltaMassDStar{"deltaMassDStar", 0.04, "invariant-mass delta with respect to the D* mass for B0 -> D*pi"};
  Configurable<float> pTMinBeautyBachelor{"pTMinBeautyBachelor", 0.5, "minimum pT for bachelor pion track used to build b-hadron candidates"};
  Configurable<float> pTMinSoftPion{"pTMinSoftPion", 0.1, "minimum pT for soft pion track used to build D* mesons in the b-hadron decay chain"};
  Configurable<std::vector<double>> pTBinsTrack{"pTBinsTrack", std::vector<double>{hf_cuts_single_track::vecBinsPtTrack}, "track pT bin limits for DCAXY pT-dependent cut"};
  Configurable<LabeledArray<double>> cutsTrackBeauty3Prong{"cutsTrackBeauty3Prong", {hf_cuts_single_track::cutsTrack[0], hf_cuts_single_track::nBinsPtTrack, hf_cuts_single_track::nCutVarsTrack, hf_cuts_single_track::labelsPtTrack, hf_cuts_single_track::labelsCutVarTrack}, "Single-track selections per pT bin for 3-prong beauty candidates"};
  Configurable<LabeledArray<double>> cutsTrackBeauty4Prong{"cutsTrackBeauty4Prong", {hf_cuts_single_track::cutsTrack[0], hf_cuts_single_track::nBinsPtTrack, hf_cuts_single_track::nCutVarsTrack, hf_cuts_single_track::labelsPtTrack, hf_cuts_single_track::labelsCutVarTrack}, "Single-track selections per pT bin for 4-prong beauty candidates"};
  std::array<LabeledArray<double>, 2> cutsSingleTrackBeauty;

  // parameters for femto triggers
  Configurable<float> femtoMaxRelativeMomentum{"femtoMaxRelativeMomentum", 2., "Maximal allowed value for relative momentum between charm-proton pairs in GeV/c"};
  Configurable<float> femtoMinProtonPt{"femtoMinProtonPt", 0.5, "Minimal required transverse momentum for proton in GeV/c"};
  Configurable<bool> femtoProtonOnlyTOF{"femtoProtonOnlyTOF", false, "Use only TOF information for proton identification if true"};
  Configurable<float> femtoMaxNsigmaProton{"femtoMaxNsigmaProton", 3., "Maximum value for PID proton Nsigma for femto triggers"};

  // parameters for photon triggers
  Configurable<float> photonMinCosPA{"photonMinCosPA", 0.85, "Minimal required cosine of pointing angle for photons"};

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
  Configurable<long> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB. Exceptions: > 0 for the specific timestamp, 0 gets the run dependent timestamp"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<string> ccdbPathTPC{"ccdbPathTPC", "Users/i/iarsene/Calib/TPCpostCalib", "base path to the ccdb object"};

  // TPC PID calibrations
  Configurable<int> setTPCCalib{"setTPCCalib", 0, "0 is not use re-calibrations, 1 is compute TPC post-calibrated n-sigmas, 2 is using TPC Spline"};
  Configurable<std::string> ccdbBBProton{"ccdbBBProton", "Users/l/lserksny/PIDProton", "Path to the CCDB ocject for proton BB param"};
  Configurable<std::string> ccdbBBAntiProton{"ccdbBBAntiProton", "Users/l/lserksny/PIDAntiProton", "Path to the CCDB ocject for antiproton BB param"};
  Configurable<std::string> ccdbBBPion{"ccdbBBPion", "Users/l/lserksny/PIDPion", "Path to the CCDB ocject for Pion BB param"};
  Configurable<std::string> ccdbBBAntiPion{"ccdbBBAntiPion", "Users/l/lserksny/PIDAntiPion", "Path to the CCDB ocject for antiPion BB param"};
  Configurable<std::string> ccdbBBKaon{"ccdbBBKaon", "Users/l/lserksny/PIDPion", "Path to the CCDB ocject for Kaon BB param"};
  Configurable<std::string> ccdbBBAntiKaon{"ccdbBBAntiKaon", "Users/l/lserksny/PIDAntiPion", "Path to the CCDB ocject for antiKaon BB param"};

  // parameter for Optimisation Tree
  Configurable<bool> applyOptimisation{"applyOptimisation", false, "Flag to enable or disable optimisation"};
  int currentRun = 0; // needed to detect if the run changed and trigger update of calibrations etc.

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
  std::array<std::shared_ptr<TH2>, kNCharmParticles + 3> hMassVsPtC{}; // +3 for resonances (D*+, D*0, Ds*+)
  std::array<std::shared_ptr<TH2>, kNCharmParticles> hMassVsPhiC{};
  std::shared_ptr<TH2> hProtonTPCPID, hProtonTOFPID;
  std::array<std::shared_ptr<TH1>, kNCharmParticles> hBDTScoreBkg{};
  std::array<std::shared_ptr<TH1>, kNCharmParticles> hBDTScorePrompt{};
  std::array<std::shared_ptr<TH1>, kNCharmParticles> hBDTScoreNonPrompt{};
  std::shared_ptr<TH1> hGammaSelected, hGammaEtaBefore, hGammaEtaAfter;
  std::shared_ptr<TH2> hGammaArmPodBefore, hGammaArmPodAfter;

  // Histograms of TPC calibration for pion and proton
  std::array<TH3F*, 2> hMapPion = {nullptr, nullptr};
  std::array<TH3F*, 2> hMapProton = {nullptr, nullptr};
  std::array<std::vector<double>, 2> hBBProton{};
  std::array<std::vector<double>, 2> hBBPion{};
  std::array<std::vector<double>, 2> hBBKaon{};
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
  o2::base::Propagator::MatCorrType noMatCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  void init(o2::framework::InitContext&)
  {
    cutsSingleTrackBeauty = {cutsTrackBeauty3Prong, cutsTrackBeauty4Prong};

    hProcessedEvents = registry.add<TH1>("fProcessedEvents", "HF - event filtered;;counts", HistType::kTH1F, {{kNtriggersHF + 2, -0.5, kNtriggersHF + 1.5}});
    std::array<std::string, kNtriggersHF + 2> eventTitles = {"all", "rejected", "w/ high-#it{p}_{T} 2p charm", "w/ high-#it{p}_{T} 3p charm", "w/ 3p beauty", "w/ 4p beauty", "w/ 2p femto", "w/ 3p femto", "w/ 2p double charm", "w/ 3p double charm", "w/ 2p and 3p double charm", "w/ 2p soft gamma", "w/ 3p soft gamma"};
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
        if (activateQA > 2) {
          hMassVsPhiC[iCharmPart] = registry.add<TH2>(Form("fMassVsPhi%s", charmParticleNames[iCharmPart].data()), Form("#it{M} vs. #varphi distribution of triggered %s candidates;#varphi;#it{M} (GeV/#it{c}^{2});counts", charmParticleNames[iCharmPart].data()), HistType::kTH2F, {phiAxis, massAxisC[iCharmPart]});
        }
      }
      hMassVsPtC[kNCharmParticles] = registry.add<TH2>("fMassVsPtDStarPlus", "#it{M} vs. #it{p}_{T} distribution of triggered DStarPlus candidates;#it{p}_{T} (GeV/#it{c});#it{M} (GeV/#it{c}^{2});counts", HistType::kTH2F, {ptAxis, massAxisC[kNCharmParticles]});
      hMassVsPtC[kNCharmParticles + 1] = registry.add<TH2>("fMassVsPtDStarZero", "#it{M} vs. #it{p}_{T} distribution of triggered DStarZero candidates;#it{p}_{T} (GeV/#it{c});#it{M} (GeV/#it{c}^{2});counts", HistType::kTH2F, {ptAxis, massAxisC[kNCharmParticles + 1]});
      hMassVsPtC[kNCharmParticles + 2] = registry.add<TH2>("fMassVsPtDStarS", "#it{M} vs. #it{p}_{T} distribution of triggered DStarS candidates;#it{p}_{T} (GeV/#it{c});#it{M} (GeV/#it{c}^{2});counts", HistType::kTH2F, {ptAxis, massAxisC[kNCharmParticles + 2]});
      for (int iBeautyPart{0}; iBeautyPart < kNBeautyParticles; ++iBeautyPart) {
        hMassVsPtB[iBeautyPart] = registry.add<TH2>(Form("fMassVsPt%s", beautyParticleNames[iBeautyPart].data()), Form("#it{M} vs. #it{p}_{T} distribution of triggered %s candidates;#it{p}_{T} (GeV/#it{c});#it{M} (GeV/#it{c}^{2});counts", beautyParticleNames[iBeautyPart].data()), HistType::kTH2F, {ptAxis, massAxisB[iBeautyPart]});
      }
      if (activateQA > 1) {
        hProtonTPCPID = registry.add<TH2>("fProtonTPCPID", "#it{N}_{#sigma}^{TPC} vs. #it{p} for selected protons;#it{p} (GeV/#it{c});#it{N}_{#sigma}^{TPC}", HistType::kTH2F, {pAxis, nSigmaAxis});
        hProtonTOFPID = registry.add<TH2>("fProtonTOFPID", "#it{N}_{#sigma}^{TOF} vs. #it{p} for selected protons;#it{p} (GeV/#it{c});#it{N}_{#sigma}^{TOF}", HistType::kTH2F, {pAxis, nSigmaAxis});
        hGammaSelected = registry.add<TH1>("fGammaSelected", "Selections for converted gamma;;counts", HistType::kTH1F, {{7, -0.5, 6.5}});
        hGammaEtaBefore = registry.add<TH1>("fGammaEtaBefore", "#eta of converted gamma before selections;;counts", HistType::kTH1F, {etaAxis});
        hGammaEtaAfter = registry.add<TH1>("hGammaEtaAfter", "#eta of converted gamma after selections;;counts", HistType::kTH1F, {etaAxis});
        hGammaArmPodBefore = registry.add<TH2>("fGammaAPbefore", "Armenteros Podolanski plot for converted gamma, before selections;#it{#alpha};#it{q}_{T} (GeV/#it{c})", HistType::kTH2F, {alphaAxis, qtAxis});
        hGammaArmPodAfter = registry.add<TH2>("fGammaAPafter", "Armenteros Podolanski plot for converted gamma, after selections;#it{#alpha};#it{q}_{T} (GeV/#it{c})", HistType::kTH2F, {alphaAxis, qtAxis});
      }
    }

    ccdb->setURL(url.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    ccdbApi.init(url);

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
          sessionML[iCharmPart].reset(InitONNXSession(onnxFiles[iCharmPart], charmParticleNames[iCharmPart], envML[iCharmPart], sessionOptions[iCharmPart], inputShapesML[iCharmPart], dataTypeML[iCharmPart], loadModelsFromCCDB, ccdbApi, mlModelPathCCDB.value, timestampCCDB));
        }
      }
    }
    // safety for optimisation tree
    if (applyOptimisation && !applyML) {
      LOG(fatal) << "Can't apply optimisation if ML is not applied.";
    }
  }

  using BigTracksMCPID = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::McTrackLabels>;

  Filter trackFilter = requireGlobalTrackWoDCAInFilter();
  using BigTracksPID = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr>>;

  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;
  Preslice<aod::V0Datas> v0sPerCollision = aod::v0data::collisionId;
  Preslice<aod::Hf2Prongs> hf2ProngPerCollision = aod::track_association::collisionId;
  Preslice<aod::Hf3Prongs> hf3ProngPerCollision = aod::track_association::collisionId;

  void process(aod::Collisions const& collisions,
               aod::BCsWithTimestamps const&,
               aod::V0Datas const& theV0s,
               aod::Hf2Prongs const& cand2Prongs,
               aod::Hf3Prongs const& cand3Prongs,
               aod::TrackAssoc const& trackIndices,
               BigTracksPID const& tracks)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();

      if (applyOptimisation) {
        optimisationTreeCollisions(thisCollId);
      }

      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();

      if (applyML && (loadModelsFromCCDB && timestampCCDB == 0) && !sessionML[kD0]) {
        for (auto iCharmPart{0}; iCharmPart < kNCharmParticles; ++iCharmPart) {
          if (onnxFiles[iCharmPart] != "") {
            sessionML[iCharmPart].reset(InitONNXSession(onnxFiles[iCharmPart], charmParticleNames[iCharmPart], envML[iCharmPart], sessionOptions[iCharmPart], inputShapesML[iCharmPart], dataTypeML[iCharmPart], loadModelsFromCCDB, ccdbApi, mlModelPathCCDB.value, bc.timestamp()));
          }
        }
      }

      // needed for track propagation
      if (currentRun != bc.runNumber()) {
        o2::parameters::GRPMagField* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(ccdbPathGrpMag, bc.timestamp());
        o2::base::Propagator::initFieldFromGRP(grpo);

        // needed for TPC PID postcalibrations
        if (setTPCCalib == 1) {
          auto calibList = ccdb->getForTimeStamp<TList>(ccdbPathTPC.value, bc.timestamp());
          if (!calibList) {
            LOG(fatal) << "Can not find the TPC Post Calibration object!";
          }

          hMapPion[0] = (TH3F*)calibList->FindObject("mean_map_pion");
          hMapPion[1] = (TH3F*)calibList->FindObject("sigma_map_pion");
          hMapProton[0] = (TH3F*)calibList->FindObject("mean_map_proton");
          hMapProton[1] = (TH3F*)calibList->FindObject("sigma_map_proton");

          if (!hMapPion[0] || !hMapPion[1] || !hMapProton[0] || !hMapProton[1]) {
            LOG(fatal) << "Can not find histograms!";
          }
        } else if (setTPCCalib > 1) {

          hBBProton[0] = setValuesBB(ccdbApi, bc, ccdbBBProton);
          hBBProton[1] = setValuesBB(ccdbApi, bc, ccdbBBAntiProton);
          hBBPion[0] = setValuesBB(ccdbApi, bc, ccdbBBPion);
          hBBPion[1] = setValuesBB(ccdbApi, bc, ccdbBBAntiPion);
          hBBKaon[0] = setValuesBB(ccdbApi, bc, ccdbBBKaon);
          hBBKaon[1] = setValuesBB(ccdbApi, bc, ccdbBBAntiKaon);
        }

        currentRun = bc.runNumber();
      }

      hProcessedEvents->Fill(0);

      // collision process loop
      bool keepEvent[kNtriggersHF]{false};
      //

      std::vector<std::vector<long>> indicesDau2Prong{};

      auto cand2ProngsThisColl = cand2Prongs.sliceBy(hf2ProngPerCollision, thisCollId);
      for (const auto& cand2Prong : cand2ProngsThisColl) {                                // start loop over 2 prongs
        if (!TESTBIT(cand2Prong.hfflag(), o2::aod::hf_cand_2prong::DecayType::D0ToPiK)) { // check if it's a D0
          continue;
        }

        auto trackPos = cand2Prong.prong0_as<BigTracksPID>(); // positive daughter
        auto trackNeg = cand2Prong.prong1_as<BigTracksPID>(); // negative daughter

        auto preselD0 = isDzeroPreselected(trackPos, trackNeg, nsigmaTPCPionKaonDzero, nsigmaTOFPionKaonDzero, setTPCCalib, hMapPion, hBBPion, hBBKaon);
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
            auto scores = PredictONNX(inputFeaturesD0, sessionML[kD0], inputShapesML[kD0]);
            tagBDT = isBDTSelected(scores, thresholdBDTScores[kD0]);
            for (int iScore{0}; iScore < 3; ++iScore) {
              scoresToFill[iScore] = scores[iScore];
            }
          } else if (dataTypeML[kD0] == 11) {
            auto scores = PredictONNX(inputFeaturesDoD0, sessionML[kD0], inputShapesML[kD0]);
            tagBDT = isBDTSelected(scores, thresholdBDTScores[kD0]);
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
        auto phi2Prong = RecoDecay::phi(pVec2Prong);

        if (applyOptimisation) {
          optimisationTreeCharm(thisCollId, pdg::Code::kD0, pt2Prong, scoresToFill[0], scoresToFill[1], scoresToFill[2]);
        }

        auto selD0 = isSelectedD0InMassRange(pVecPos, pVecNeg, pt2Prong, phi2Prong, preselD0, deltaMassCharmHadronForBeauty, activateQA, hMassVsPtC[kD0], hMassVsPhiC[kD0]);

        if (pt2Prong >= pTThreshold2Prong) {
          keepEvent[kHighPt2P] = true;
          if (activateQA) {
            hCharmHighPt[kD0]->Fill(pt2Prong);
          }
        } // end high-pT selection

        if (isCharmTagged) {
          indicesDau2Prong.push_back(std::vector<long>{trackPos.globalIndex(), trackNeg.globalIndex()});
        } // end multi-charm selection

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
            int isTrackSelected = isSelectedTrackForBeauty(trackParThird, dcaThird, pTMinSoftPion, pTMinBeautyBachelor, pTBinsTrack, cutsSingleTrackBeauty[kBeauty3P - 2]);
            if (isTrackSelected && ((TESTBIT(selD0, 0) && track.sign() < 0) || (TESTBIT(selD0, 1) && track.sign() > 0))) {
              auto massCand = RecoDecay::m(std::array{pVec2Prong, pVecThird}, std::array{massD0, massPi});
              auto pVecBeauty3Prong = RecoDecay::pVec(pVec2Prong, pVecThird);
              auto ptCand = RecoDecay::pt(pVecBeauty3Prong);
              if (isTrackSelected == kRegular && std::abs(massCand - massBPlus) <= deltaMassBPlus) {
                keepEvent[kBeauty3P] = true;
                // fill optimisation tree for D0
                if (applyOptimisation) {
                  optimisationTreeBeauty(thisCollId, pdg::Code::kD0, pt2Prong, scoresToFill[0], scoresToFill[1], scoresToFill[2], dcaThird[0]);
                }
                if (activateQA) {
                  hMassVsPtB[kBplus]->Fill(ptCand, massCand);
                }
              } else if (std::abs(massCand - massDStar) <= deltaMassDStar) { // additional check for B0->D*pi polarization studies
                if (activateQA) {
                  hMassVsPtC[kNCharmParticles]->Fill(ptCand, massCand);
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

                  if (track.sign() * trackB.sign() < 0 && isSelectedTrackForBeauty(trackParFourth, dcaFourth, pTMinSoftPion, pTMinBeautyBachelor, pTBinsTrack, cutsSingleTrackBeauty[kBeauty3P - 2]) == kRegular) {
                    auto massCandB0 = RecoDecay::m(std::array{pVec2Prong, pVecThird, pVecFourth}, std::array{massD0, massPi, massPi});
                    if (std::abs(massCandB0 - massB0) <= deltaMassB0) {
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
          } // end beauty selection

          // 2-prong femto
          if (!keepEvent[kFemto2P] && isCharmTagged && track.collisionId() == thisCollId) {
            bool isProton = isSelectedProton4Femto(track, trackParThird, femtoMinProtonPt, femtoMaxNsigmaProton, femtoProtonOnlyTOF, setTPCCalib, hMapProton, hBBProton, activateQA, hProtonTPCPID, hProtonTOFPID);
            if (isProton) {
              float relativeMomentum = computeRelativeMomentum(pVecThird, pVec2Prong, massD0);
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
        for (auto& gamma : v0sThisCollision) {
          if (!keepEvent[kGammaCharm2P] && (isCharmTagged || isBeautyTagged) && (TESTBIT(selD0, 0) || (TESTBIT(selD0, 1)))) {
            float V0CosinePA = gamma.v0cosPA(collision.posX(), collision.posY(), collision.posZ());
            bool isGamma = isSelectedGamma(gamma, photonMinCosPA, V0CosinePA, activateQA, hGammaSelected, hGammaEtaBefore, hGammaEtaAfter, hGammaArmPodBefore, hGammaArmPodAfter);
            if (isGamma) {
              std::array<float, 3> gammaVec = {gamma.px(), gamma.py(), gamma.pz()};
              auto massGammaCharm = RecoDecay::m(std::array{pVec2Prong, gammaVec}, std::array{massD0, massGamma});
              if (massGammaCharm < 3.0) { // remove candidates with invariant mass above 3 GeV
                if (activateQA) {
                  auto pVecReso2Prong = RecoDecay::pVec(pVec2Prong, gammaVec);
                  auto ptCand = RecoDecay::pt(pVecReso2Prong);
                  hMassVsPtC[kNCharmParticles + 1]->Fill(ptCand, massGammaCharm);
                }
                keepEvent[kGammaCharm2P] = true;
              }
            }
          }
        } // end gamma selection

      } // end loop over 2-prong candidates

      std::vector<std::vector<long>> indicesDau3Prong{};
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
          is3Prong[0] = isDplusPreselected(trackSecond, nsigmaTPCKaon3Prong, nsigmaTOFKaon3Prong, setTPCCalib, hMapPion, hBBKaon);
        }
        if (is3Prong[1]) { // Ds preselections
          is3Prong[1] = isDsPreselected(pVecFirst, pVecThird, pVecSecond, trackSecond, nsigmaTPCKaon3Prong, nsigmaTOFKaon3Prong, setTPCCalib, hMapPion, hBBKaon);
        }
        if (is3Prong[2] || is3Prong[3]) { // charm baryon preselections
          auto presel = isCharmBaryonPreselected(trackFirst, trackThird, trackSecond, nsigmaTPCProtonLc, nsigmaTOFProtonLc, nsigmaTPCKaon3Prong, nsigmaTOFKaon3Prong, setTPCCalib, hMapProton, hBBProton, hMapPion, hBBKaon);
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
              auto scores = PredictONNX(inputFeatures, sessionML[iCharmPart + 1], inputShapesML[iCharmPart + 1]);
              tagBDT = isBDTSelected(scores, thresholdBDTScores[iCharmPart + 1]);
              for (int iScore{0}; iScore < 3; ++iScore) {
                scoresToFill[iCharmPart][iScore] = scores[iScore];
              }
            } else if (dataTypeML[iCharmPart + 1] == 11) {
              auto scores = PredictONNX(inputFeaturesD, sessionML[iCharmPart + 1], inputShapesML[iCharmPart + 1]);
              tagBDT = isBDTSelected(scores, thresholdBDTScores[iCharmPart + 1]);
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
          indicesDau3Prong.push_back(std::vector<long>{trackFirst.globalIndex(), trackSecond.globalIndex(), trackThird.globalIndex()});
        } // end multiple 3-prong selection

        auto pVec3Prong = RecoDecay::pVec(pVecFirst, pVecSecond, pVecThird);
        auto pt3Prong = RecoDecay::pt(pVec3Prong);
        auto phi3Prong = RecoDecay::phi(pVec3Prong);
        float sign3Prong = -1 * trackFirst.sign() * trackSecond.sign() * trackThird.sign();

        std::array<int8_t, kNCharmParticles - 1> is3ProngInMass{0};
        if (is3Prong[0]) {
          is3ProngInMass[0] = isSelectedDplusInMassRange(pVecFirst, pVecThird, pVecSecond, pt3Prong, phi3Prong, deltaMassCharmHadronForBeauty, activateQA, hMassVsPtC[kDplus], hMassVsPhiC[kDplus]);
          if (applyOptimisation) {
            optimisationTreeCharm(thisCollId, pdg::Code::kDPlus, pt3Prong, scoresToFill[0][0], scoresToFill[0][1], scoresToFill[0][2]);
          }
        }
        if (is3Prong[1]) {
          is3ProngInMass[1] = isSelectedDsInMassRange(pVecFirst, pVecThird, pVecSecond, pt3Prong, phi3Prong, is3Prong[1], deltaMassCharmHadronForBeauty, activateQA, hMassVsPtC[kDs], hMassVsPhiC[kDs]);
          if (applyOptimisation) {
            optimisationTreeCharm(thisCollId, pdg::Code::kDS, pt3Prong, scoresToFill[1][0], scoresToFill[1][1], scoresToFill[1][2]);
          }
        }
        if (is3Prong[2]) {
          is3ProngInMass[2] = isSelectedLcInMassRange(pVecFirst, pVecThird, pVecSecond, pt3Prong, phi3Prong, is3Prong[2], deltaMassCharmHadronForBeauty, activateQA, hMassVsPtC[kLc], hMassVsPhiC[kLc]);
          if (applyOptimisation) {
            optimisationTreeCharm(thisCollId, pdg::Code::kLambdaCPlus, pt3Prong, scoresToFill[2][0], scoresToFill[2][1], scoresToFill[2][2]);
          }
        }
        if (is3Prong[3]) {
          is3ProngInMass[3] = isSelectedXicInMassRange(pVecFirst, pVecThird, pVecSecond, pt3Prong, phi3Prong, is3Prong[3], deltaMassCharmHadronForBeauty, activateQA, hMassVsPtC[kXic], hMassVsPhiC[kXic]);
          if (applyOptimisation) {
            optimisationTreeCharm(thisCollId, pdg::Code::kXiCPlus, pt3Prong, scoresToFill[3][0], scoresToFill[3][1], scoresToFill[3][2]);
          }
        }

        if (pt3Prong >= pTThreshold3Prong) {
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
          float deltaMassHypos[kNBeautyParticles - 2] = {deltaMassB0, deltaMassBs, deltaMassLb, deltaMassXib};
          if (track.sign() * sign3Prong < 0 && isSelectedTrackForBeauty(trackParFourth, dcaFourth, pTMinBeautyBachelor, pTMinBeautyBachelor, pTBinsTrack, cutsSingleTrackBeauty[kBeauty4P - 2]) == kRegular) {
            for (int iHypo{0}; iHypo < kNBeautyParticles - 2 && !keepEvent[kBeauty4P]; ++iHypo) {
              if (isBeautyTagged[iHypo] && (TESTBIT(is3ProngInMass[iHypo], 0) || TESTBIT(is3ProngInMass[iHypo], 1))) {
                auto massCandB = RecoDecay::m(std::array{pVec3Prong, pVecFourth}, std::array{massCharmHypos[iHypo], massPi});
                if (std::abs(massCandB - massBeautyHypos[iHypo]) <= deltaMassHypos[iHypo]) {
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
          bool isProton = isSelectedProton4Femto(track, trackParFourth, femtoMinProtonPt, femtoMaxNsigmaProton, femtoProtonOnlyTOF, setTPCCalib, hMapProton, hBBProton, activateQA, hProtonTPCPID, hProtonTOFPID);
          if (isProton && track.collisionId() == thisCollId) {
            for (int iHypo{0}; iHypo < kNCharmParticles - 1 && !keepEvent[kFemto3P]; ++iHypo) {
              if (isCharmTagged[iHypo]) {
                float relativeMomentum = computeRelativeMomentum(pVecFourth, pVec3Prong, massCharmHypos[iHypo]);
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

        // 3-prong with Gamma (conversion photon)
        auto v0sThisCollision = theV0s.sliceBy(v0sPerCollision, thisCollId);
        for (auto& gamma : v0sThisCollision) {
          if (!keepEvent[kGammaCharm3P] && (isCharmTagged[kDs - 1] || isBeautyTagged[kDs - 1]) && (TESTBIT(is3ProngInMass[kDs - 1], 0) || TESTBIT(is3ProngInMass[kDs - 1], 1))) {
            float V0CosinePA = gamma.v0cosPA(collision.posX(), collision.posY(), collision.posZ());
            bool isGamma = isSelectedGamma(gamma, photonMinCosPA, V0CosinePA, activateQA, hGammaSelected, hGammaEtaBefore, hGammaEtaAfter, hGammaArmPodBefore, hGammaArmPodAfter);
            if (isGamma) {
              std::array<float, 3> gammaVec = {gamma.px(), gamma.py(), gamma.pz()};
              auto massGammaCharm = RecoDecay::m(std::array{pVec3Prong, gammaVec}, std::array{massDs, massGamma});
              if (massGammaCharm < 3.) { // remove candidates with invariant mass above some value (TODO: set me as a configurable)
                if (activateQA) {
                  auto pVecReso3Prong = RecoDecay::pVec(pVec3Prong, gammaVec);
                  auto ptCand = RecoDecay::pt(pVecReso3Prong);
                  hMassVsPtC[kNCharmParticles + 2]->Fill(ptCand, massGammaCharm);
                }
                keepEvent[kGammaCharm3P] = true;
              }
            }
          }
        } // end gamma selection

      } // end loop over 3-prong candidates

      auto n2Prongs = computeNumberOfCandidates(indicesDau2Prong);
      auto n3Prongs = computeNumberOfCandidates(indicesDau3Prong);
      indicesDau2Prong.insert(indicesDau2Prong.end(), indicesDau3Prong.begin(), indicesDau3Prong.end());
      auto n23Prongs = computeNumberOfCandidates(indicesDau2Prong);

      if (activateQA) {
        hN2ProngCharmCand->Fill(n2Prongs);
        hN3ProngCharmCand->Fill(n3Prongs);
      }

      if (n2Prongs > 1) {
        keepEvent[kDoubleCharm2P] = true;
      }
      if (n3Prongs > 1) {
        keepEvent[kDoubleCharm3P] = true;
      }
      if (n23Prongs > 1) {
        keepEvent[kDoubleCharmMix] = true;
      }

      tags(keepEvent[kHighPt2P], keepEvent[kHighPt3P], keepEvent[kBeauty3P], keepEvent[kBeauty4P], keepEvent[kFemto2P], keepEvent[kFemto3P], keepEvent[kDoubleCharm2P], keepEvent[kDoubleCharm3P], keepEvent[kDoubleCharmMix], keepEvent[kGammaCharm2P], keepEvent[kGammaCharm3P]);

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
