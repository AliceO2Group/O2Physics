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

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hffilters;
using namespace hf_cuts_single_track;
using namespace hf_cuts_bdt_multiclass;

struct AddCollisionId {

  Produces<o2::aod::Colls2Prong> colls2Prong;
  Produces<o2::aod::Colls3Prong> colls3Prong;

  void process(aod::Hf2Prongs const& cand2Prongs,
               aod::Hf3Prongs const& cand3Prongs,
               aod::Tracks const&)
  {
    for (const auto& cand2Prong : cand2Prongs) {
      colls2Prong(cand2Prong.prong0_as<aod::Tracks>().collisionId());
    }
    for (const auto& cand3Prong : cand3Prongs) {
      colls3Prong(cand3Prong.prong0_as<aod::Tracks>().collisionId());
    }
  }
};

struct HfFilter { // Main struct for HF triggers

  Produces<aod::HfFilters> tags;
  Produces<aod::HFTrigTrain2P> train2P;
  Produces<aod::HFTrigTrain3P> train3P;
  Produces<aod::HFOptimisationTreeBeauty> optimisationTreeBeauty;
  Produces<aod::HFOptimisationTreeCharm> optimisationTreeCharm;
  Produces<aod::HFOptimisationTreeFemto> optimisationTreeFemto;
  Produces<aod::HFOptimisationTreeCollisions> optimisationTreeCollisions;

  Configurable<bool> activateQA{"activateQA", false, "flag to enable QA histos"};

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
  Configurable<bool> femtoProtonOnlyTOF{"femtoProtonOnlyTOF", true, "Use only TOF information for proton identification if true"};
  Configurable<float> femtoMaxNsigmaProton{"femtoMaxNsigmaProton", 3., "Maximum value for PID proton Nsigma for femto triggers"};

  // parameters for production of training samples
  Configurable<bool> fillSignal{"fillSignal", true, "Flag to fill derived tables with signal for ML trainings"};
  Configurable<bool> fillBackground{"fillBackground", true, "Flag to fill derived tables with background for ML trainings"};
  Configurable<double> donwSampleBkgFactor{"donwSampleBkgFactor", 1., "Fraction of background candidates to keep for ML trainings"};

  // parameters for ML application with ONNX
  Configurable<bool> applyML{"applyML", false, "Flag to enable or disable ML application"};
  Configurable<std::vector<double>> pTBinsBDT{"pTBinsBDT", std::vector<double>{hf_cuts_bdt_multiclass::vecBinsPt}, "track pT bin limits for BDT cut"};

  Configurable<std::string> onnxFileD0ToKPiConf{"onnxFileD0ToKPiConf", "XGBoostModel.onnx", "ONNX file for ML model for D0 candidates"};
  Configurable<LabeledArray<double>> thresholdBDTScoreD0ToKPi{"thresholdBDTScoreD0ToKPi", {hf_cuts_bdt_multiclass::cuts[0], hf_cuts_bdt_multiclass::nBinsPt, hf_cuts_bdt_multiclass::nCutBdtScores, hf_cuts_bdt_multiclass::labelsPt, hf_cuts_bdt_multiclass::labelsCutBdt}, "Threshold values for BDT output scores of D0 candidates"};

  Configurable<std::string> onnxFileDPlusToPiKPiConf{"onnxFileDPlusToPiKPiConf", "", "ONNX file for ML model for D+ candidates"};
  Configurable<LabeledArray<double>> thresholdBDFScoreDPlusToPiKPi{"thresholdBDFScoreDPlusToPiKPi", {hf_cuts_bdt_multiclass::cuts[0], hf_cuts_bdt_multiclass::nBinsPt, hf_cuts_bdt_multiclass::nCutBdtScores, hf_cuts_bdt_multiclass::labelsPt, hf_cuts_bdt_multiclass::labelsCutBdt}, "Threshold values for BDT output scores of D+ candidates"};

  Configurable<std::string> onnxFileDSToPiKKConf{"onnxFileDSToPiKKConf", "", "ONNX file for ML model for Ds+ candidates"};
  Configurable<LabeledArray<double>> thresholdBDFScoreDSToPiKK{"thresholdBDFScoreDSToPiKK", {hf_cuts_bdt_multiclass::cuts[0], hf_cuts_bdt_multiclass::nBinsPt, hf_cuts_bdt_multiclass::nCutBdtScores, hf_cuts_bdt_multiclass::labelsPt, hf_cuts_bdt_multiclass::labelsCutBdt}, "Threshold values for BDT output scores of Ds+ candidates"};

  Configurable<std::string> onnxFileLcToPiKPConf{"onnxFileLcToPiKPConf", "", "ONNX file for ML model for Lc+ candidates"};
  Configurable<LabeledArray<double>> thresholdBDFScoreLcToPiKP{"thresholdBDFScoreLcToPiKP", {hf_cuts_bdt_multiclass::cuts[0], hf_cuts_bdt_multiclass::nBinsPt, hf_cuts_bdt_multiclass::nCutBdtScores, hf_cuts_bdt_multiclass::labelsPt, hf_cuts_bdt_multiclass::labelsCutBdt}, "Threshold values for BDT output scores of Lc+ candidates"};

  Configurable<std::string> onnxFileXicToPiKPConf{"onnxFileXicToPiKPConf", "", "ONNX file for ML model for Xic+ candidates"};
  Configurable<LabeledArray<double>> thresholdBDFScoreXicToPiKP{"thresholdBDFScoreXicToPiKP", {hf_cuts_bdt_multiclass::cuts[0], hf_cuts_bdt_multiclass::nBinsPt, hf_cuts_bdt_multiclass::nCutBdtScores, hf_cuts_bdt_multiclass::labelsPt, hf_cuts_bdt_multiclass::labelsCutBdt}, "Threshold values for BDT output scores of Xic+ candidates"};

  // CCDB configuration
  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> mlModelPathCCDB{"mlModelPathCCDB", "Analysis/PWGHF/ML/HFTrigger/", "Path on CCDB"};
  Configurable<long> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB. Exceptions: > 0 for the specific timestamp, 0 gets the run dependent timestamp"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};

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
  std::array<std::shared_ptr<TH2>, kNCharmParticles + 3> hMassVsPtC{}; // +3 for resonances (D*+, D*0, Ds*+)
  std::shared_ptr<TH2> hProtonTPCPID, hProtonTOFPID;
  std::array<std::shared_ptr<TH1>, kNCharmParticles> hBDTScoreBkg{};
  std::array<std::shared_ptr<TH1>, kNCharmParticles> hBDTScorePrompt{};
  std::array<std::shared_ptr<TH1>, kNCharmParticles> hBDTScoreNonPrompt{};
  std::shared_ptr<TH1> hGammaSelected;
  std::shared_ptr<TH2> hGammaAPbefore, hGammaAPafter;

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
        hCharmHighPt[iCharmPart] = registry.add<TH1>(Form("f%sHighPt", charmParticleNames[iCharmPart].data()), Form("#it{p}_{T} distribution of triggered high-#it{p}_{T} %s candidates;#it{p}_{T} (GeV/#it{c});counts", charmParticleNames[iCharmPart].data()), HistType::kTH1F, {{100, 0., 50.}});
        hCharmProtonKstarDistr[iCharmPart] = registry.add<TH1>(Form("f%sProtonKstarDistr", charmParticleNames[iCharmPart].data()), Form("#it{k}* distribution of triggered p#minus%s pairs;#it{k}* (GeV/#it{c});counts", charmParticleNames[iCharmPart].data()), HistType::kTH1F, {{100, 0., 1.}});
        hMassVsPtC[iCharmPart] = registry.add<TH2>(Form("fMassVsPt%s", charmParticleNames[iCharmPart].data()), Form("#it{M} vs. #it{p}_{T} distribution of triggered %s candidates;#it{p}_{T} (GeV/#it{c});#it{M} (GeV/#it{c}^{2});counts", charmParticleNames[iCharmPart].data()), HistType::kTH2F, {{100, 0., 50.}, {300, 1.60, 2.60}});
        if (applyML) {
          hBDTScoreBkg[iCharmPart] = registry.add<TH1>(Form("f%sBDTScoreBkgDistr", charmParticleNames[iCharmPart].data()), Form("BDT background score distribution for %s;BDT background score;counts", charmParticleNames[iCharmPart].data()), HistType::kTH1F, {{100, 0., 1.}});
          hBDTScorePrompt[iCharmPart] = registry.add<TH1>(Form("f%sBDTScorePromptDistr", charmParticleNames[iCharmPart].data()), Form("BDT prompt score distribution for %s;BDT prompt score;counts", charmParticleNames[iCharmPart].data()), HistType::kTH1F, {{100, 0., 1.}});
          hBDTScoreNonPrompt[iCharmPart] = registry.add<TH1>(Form("f%sBDTScoreNonPromptDistr", charmParticleNames[iCharmPart].data()), Form("BDT nonprompt score distribution for %s;BDT nonprompt score;counts", charmParticleNames[iCharmPart].data()), HistType::kTH1F, {{100, 0., 1.}});
        }
      }
      hMassVsPtC[kNCharmParticles] = registry.add<TH2>("fMassVsPtDStarPlus", "#it{M} vs. #it{p}_{T} distribution of triggered DStarPlus candidates;#it{p}_{T} (GeV/#it{c});#it{M} (GeV/#it{c}^{2});counts", HistType::kTH2F, {{100, 0., 50.}, {300, 1.60, 2.60}});
      hMassVsPtC[kNCharmParticles + 1] = registry.add<TH2>("fMassVsPtDStarZero", "#it{M} vs. #it{p}_{T} distribution of triggered DStarZero candidates;#it{p}_{T} (GeV/#it{c});#it{M} (GeV/#it{c}^{2});counts", HistType::kTH2F, {{100, 0., 50.}, {330, 1.90, 3.00}});
      hMassVsPtC[kNCharmParticles + 2] = registry.add<TH2>("fMassVsPtDStarS", "#it{M} vs. #it{p}_{T} distribution of triggered DStarS candidates;#it{p}_{T} (GeV/#it{c});#it{M} (GeV/#it{c}^{2});counts", HistType::kTH2F, {{100, 0., 50.}, {330, 1.90, 3.00}});
      for (int iBeautyPart{0}; iBeautyPart < kNBeautyParticles; ++iBeautyPart) {
        hMassVsPtB[iBeautyPart] = registry.add<TH2>(Form("fMassVsPt%s", beautyParticleNames[iBeautyPart].data()), Form("#it{M} vs. #it{p}_{T} distribution of triggered %s candidates;#it{p}_{T} (GeV/#it{c});#it{M} (GeV/#it{c}^{2});counts", beautyParticleNames[iBeautyPart].data()), HistType::kTH2F, {{100, 0., 50.}, {220, 4.9, 6.0}});
      }
      hProtonTPCPID = registry.add<TH2>("fProtonTPCPID", "#it{N}_{#sigma}^{TPC} vs. #it{p} for selected protons;#it{p} (GeV/#it{c});#it{N}_{#sigma}^{TPC}", HistType::kTH2F, {{100, 0., 10.}, {200, -10., 10.}});
      hProtonTOFPID = registry.add<TH2>("fProtonTOFPID", "#it{N}_{#sigma}^{TOF} vs. #it{p} for selected protons;#it{p} (GeV/#it{c});#it{N}_{#sigma}^{TOF}", HistType::kTH2F, {{100, 0., 10.}, {200, -10., 10.}});
      hGammaSelected = registry.add<TH1>("fGammaSelected", "Selections for converted gamma;;counts", HistType::kTH1F, {{7, -0.5, 6.5}});
      hGammaAPbefore = registry.add<TH2>("fGammaAPbefore", "Armenteros Podolanski plot for converted gamma, before selections;#it{#alpha};#it{q}_{T} (GeV/#it{c})", HistType::kTH2F, {{100, -1., 1.}, {100, 0., .25}});
      hGammaAPafter = registry.add<TH2>("fGammaAPafter", "Armenteros Podolanski plot for converted gamma, after selections;#it{#alpha};#it{q}_{T} (GeV/#it{c})", HistType::kTH2F, {{100, -1., 1.}, {100, 0., .25}});
    }

    ccdb->setURL(url.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    ccdbApi.init(url);

    thresholdBDTScores = {
      thresholdBDTScoreD0ToKPi,
      thresholdBDFScoreDPlusToPiKPi,
      thresholdBDFScoreDSToPiKK,
      thresholdBDFScoreLcToPiKP,
      thresholdBDFScoreXicToPiKP};

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
  /// Single-track cuts for bachelor track of beauty candidates
  /// \param track is a track
  /// \param candType candidate type (3-prong or 4-prong beauty candidate)
  /// \return 0 if track is rejected, 1 if track is soft pion, 2 if it is regular beauty
  template <typename T>
  int isSelectedTrackForBeauty(const T& track, const int candType)
  {
    auto pT = track.pt();
    auto pTBinTrack = findBin(pTBinsTrack, pT);
    if (pTBinTrack == -1) {
      return kRejected;
    }

    if (pT < pTMinSoftPion) { // soft pion should be more stringent than usual tracks
      return kRejected;
    }

    if (std::abs(track.eta()) > 0.8) {
      return kRejected;
    }

    if (std::abs(track.dcaZ()) > 2.f) {
      return kRejected;
    }

    if (std::abs(track.dcaXY()) < cutsSingleTrackBeauty[candType - 2].get(pTBinTrack, "min_dcaxytoprimary")) {
      return kRejected; // minimum DCAxy
    }
    if (std::abs(track.dcaXY()) > cutsSingleTrackBeauty[candType - 2].get(pTBinTrack, "max_dcaxytoprimary")) {
      return kRejected; // maximum DCAxy
    }

    // below only regular beauty tracks, not required for soft pions
    if (pT < pTMinBeautyBachelor) {
      return kSoftPion;
    }

    return kRegular;
  }

  /// Basic selection of gamma candidates
  /// \return true if gamma passes all cuts
  template <typename T>
  bool isSelectedGamma(const T& gamma, float GammaCosinePA)
  {
    if (activateQA) {
      hGammaSelected->Fill(0);
      hGammaAPbefore->Fill(gamma.alpha(), gamma.qtarm());
    }
    if (std::abs(gamma.eta()) > 0.8) {
      if (activateQA)
        hGammaSelected->Fill(1);
      return false;
    }

    if (gamma.v0radius() < 0. || gamma.v0radius() > 180.) {
      if (activateQA)
        hGammaSelected->Fill(2);
      return false;
    }

    if ((std::pow(gamma.alpha() / 0.95, 2) + std::pow(gamma.qtarm() / 0.05, 2)) >= 1) {
      if (activateQA)
        hGammaSelected->Fill(3);
      return false;
    }

    if (std::abs(gamma.psipair()) > 0.1) {
      if (activateQA)
        hGammaSelected->Fill(4);
      return false;
    }

    if (GammaCosinePA < 0.85) {
      if (activateQA)
        hGammaSelected->Fill(5);
      return false;
    }

    if (activateQA) {
      hGammaSelected->Fill(6);
      hGammaAPafter->Fill(gamma.alpha(), gamma.qtarm());
    }
    return true;
  }

  /// Basic selection of proton candidates
  /// \param track is a track
  /// \return true if track passes all cuts
  template <typename T>
  bool isSelectedProton4Femto(const T& track)
  {
    if (track.pt() < femtoMinProtonPt) {
      return false;
    }

    if (std::abs(track.eta()) > 0.8) {
      return false;
    }

    if (track.isGlobalTrack() != (uint8_t) true) {
      return false; // use only global tracks
    }

    float NSigmaTPC = track.tpcNSigmaPr();
    float NSigmaTOF = track.tofNSigmaPr();
    float NSigma;

    if (femtoProtonOnlyTOF) {
      NSigma = abs(NSigmaTOF);
    } else {
      NSigma = sqrt(NSigmaTPC * NSigmaTPC + NSigmaTOF * NSigmaTOF);
    }

    if (NSigma > femtoMaxNsigmaProton) {
      return false;
    }

    if (activateQA) {
      hProtonTPCPID->Fill(track.p(), NSigmaTPC);
      hProtonTOFPID->Fill(track.p(), NSigmaTOF);
    }

    return true;
  }

  /// Basic selection of proton candidates for Lc
  /// \param track is a track
  /// \return true if track passes all cuts
  template <typename T>
  bool isSelectedProton4CharmBaryons(const T& track)
  {
    float NSigmaTPC = track.tpcNSigmaPr();
    float NSigmaTOF = track.tofNSigmaPr();

    if (std::abs(NSigmaTPC) > nsigmaTPCProtonLc) {
      return false;
    }
    if (track.hasTOF() && std::abs(NSigmaTOF) > nsigmaTOFProtonLc) {
      return false;
    }

    return true;
  }

  /// Basic selection of kaon candidates for charm candidates
  /// \param track is a track
  /// \return true if track passes all cuts
  template <typename T>
  bool isSelectedKaon4Charm3Prong(const T& track)
  {
    float NSigmaTPC = track.tpcNSigmaKa();
    float NSigmaTOF = track.tofNSigmaKa();

    if (std::abs(NSigmaTPC) > nsigmaTPCKaon3Prong) {
      return false;
    }
    if (track.hasTOF() && std::abs(NSigmaTOF) > nsigmaTOFKaon3Prong) {
      return false;
    }

    return true;
  }

  /// Basic additional selection of D0 candidates
  /// \param trackPos is the positive track
  /// \param trackNeg is the negative track
  /// \return BIT(0) for D0, BIT(1) for D0bar
  template <typename T>
  int8_t isDzeroPreselected(const T& trackPos, const T& trackNeg)
  {
    int8_t retValue = 0;

    float NSigmaPiTPCPos = trackPos.tpcNSigmaPi();
    float NSigmaPiTOFPos = trackPos.tofNSigmaPi();
    float NSigmaKaTPCPos = trackPos.tpcNSigmaKa();
    float NSigmaKaTOFPos = trackPos.tofNSigmaKa();

    float NSigmaPiTPCNeg = trackNeg.tpcNSigmaPi();
    float NSigmaPiTOFNeg = trackNeg.tofNSigmaPi();
    float NSigmaKaTPCNeg = trackNeg.tpcNSigmaKa();
    float NSigmaKaTOFNeg = trackNeg.tofNSigmaKa();

    if ((std::abs(NSigmaPiTPCPos) <= nsigmaTPCPionKaonDzero && (!trackPos.hasTOF() || std::abs(NSigmaPiTOFPos) <= nsigmaTPCPionKaonDzero)) && (std::abs(NSigmaKaTPCNeg) <= nsigmaTPCPionKaonDzero && (!trackNeg.hasTOF() || std::abs(NSigmaKaTOFNeg) <= nsigmaTPCPionKaonDzero))) {
      retValue |= BIT(0);
    }
    if ((std::abs(NSigmaPiTPCNeg) <= nsigmaTPCPionKaonDzero && (!trackNeg.hasTOF() || std::abs(NSigmaPiTOFNeg) <= nsigmaTPCPionKaonDzero)) && (std::abs(NSigmaKaTPCPos) <= nsigmaTPCPionKaonDzero && (!trackPos.hasTOF() || std::abs(NSigmaKaTOFPos) <= nsigmaTPCPionKaonDzero))) {
      retValue |= BIT(1);
    }

    return retValue;
  }

  /// Basic additional selection of D+ candidates
  /// \param trackOppositeCharge is the opposite charge track momentum
  /// \return BIT(0) for Kpipi
  template <typename T>
  int8_t isDplusPreselected(const T& trackOppositeCharge)
  {
    int8_t retValue = 0;

    // check PID of opposite charge track
    if (!isSelectedKaon4Charm3Prong(trackOppositeCharge)) {
      return retValue;
    }

    retValue |= BIT(0);
    return retValue;
  }

  /// Basic additional selection of Ds candidates
  /// \param pTrackSameChargeFirst is the first same-charge track momentum
  /// \param pTrackSameChargeSecond is the second same-charge track momentum
  /// \param pTrackOppositeCharge is the opposite charge track momentum
  /// \param trackOppositeCharge is the opposite charge track
  /// \return BIT(0) for KKpi, BIT(1) for piKK
  template <typename P, typename T>
  int8_t isDsPreselected(const P& pTrackSameChargeFirst, const P& pTrackSameChargeSecond, const P& pTrackOppositeCharge, const T& trackOppositeCharge)
  {
    int8_t retValue = 0;

    // check PID of opposite charge track
    if (!isSelectedKaon4Charm3Prong(trackOppositeCharge)) {
      return retValue;
    }

    // check delta-mass for phi resonance
    auto invMassKKFirst = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackOppositeCharge}, std::array{massK, massK});
    auto invMassKKSecond = RecoDecay::m(std::array{pTrackSameChargeSecond, pTrackOppositeCharge}, std::array{massK, massK});

    if (std::abs(invMassKKFirst - massPhi) < 0.02) {
      retValue |= BIT(0);
    }
    if (std::abs(invMassKKSecond - massPhi) < 0.02) {
      retValue |= BIT(1);
    }

    return retValue;
  }

  /// Basic additional selection of Lc->pKpi and Xic->pKpi candidates
  /// \param trackSameChargeFirst is the first same-charge track
  /// \param trackSameChargeSecond is the second same-charge track
  /// \param trackOppositeCharge is the opposite charge track
  /// \return BIT(0) for pKpi, BIT(1) for piKp
  template <typename T>
  int8_t isCharmBaryonPreselected(const T& trackSameChargeFirst, const T& trackSameChargeSecond, const T& trackOppositeCharge)
  {
    int8_t retValue = 0;
    // check PID of opposite charge track
    if (!isSelectedKaon4Charm3Prong(trackOppositeCharge)) {
      return retValue;
    }

    if (isSelectedProton4CharmBaryons(trackSameChargeFirst)) {
      retValue |= BIT(0);
    }
    if (isSelectedProton4CharmBaryons(trackSameChargeSecond)) {
      retValue |= BIT(1);
    }

    return retValue;
  }

  /// Mass selection of D0 candidates to build Bplus candidates
  /// \param pTrackPos is the positive track momentum
  /// \param pTrackNeg is the negative track momentum
  /// \param ptD is the pt of the D0 meson candidate
  /// \return 1 for D0, 2 for D0bar, 3 for both
  template <typename T>
  int8_t isSelectedD0InMassRange(const T& pTrackPos, const T& pTrackNeg, const float& ptD, int8_t isSelected)
  {
    int8_t retValue = 0;
    if (TESTBIT(isSelected, 0)) {
      auto invMassD0 = RecoDecay::m(std::array{pTrackPos, pTrackNeg}, std::array{massPi, massK});
      if (activateQA) {
        hMassVsPtC[kD0]->Fill(ptD, invMassD0);
      }
      if (std::abs(invMassD0 - massD0) < deltaMassCharmHadronForBeauty) {
        retValue |= BIT(0);
      }
    }
    if (TESTBIT(isSelected, 1)) {
      auto invMassD0bar = RecoDecay::m(std::array{pTrackPos, pTrackNeg}, std::array{massK, massPi});
      if (activateQA) {
        hMassVsPtC[kD0]->Fill(ptD, invMassD0bar);
      }
      if (std::abs(invMassD0bar - massD0) < deltaMassCharmHadronForBeauty) {
        retValue |= BIT(1);
      }
    }

    return retValue;
  }

  /// Mass selection of D+ candidates to build B0 candidates
  /// \param pTrackSameChargeFirst is the first same-charge track momentum
  /// \param pTrackSameChargeFirst is the second same-charge track momentum
  /// \param pTrackSameChargeFirst is the opposite charge track momentum
  /// \param ptD is the pt of the D+ meson candidate
  /// \return BIT(0) (==1) for D+, 0 otherwise
  template <typename T>
  int8_t isSelectedDplusInMassRange(const T& pTrackSameChargeFirst, const T& pTrackSameChargeSecond, const T& pTrackOppositeCharge, const float& ptD)
  {
    auto invMassDplus = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackSameChargeSecond, pTrackOppositeCharge}, std::array{massPi, massPi, massK});
    if (activateQA) {
      hMassVsPtC[kDplus]->Fill(ptD, invMassDplus);
    }

    if (std::abs(invMassDplus - massDPlus) > deltaMassCharmHadronForBeauty) {
      return 0;
    }

    return BIT(0);
  }

  /// Mass selection of of Ds candidates to build Bs candidates
  /// \param pTrackSameChargeFirst is the first same-charge track momentum
  /// \param pTrackSameChargeFirst is the second same-charge track momentum
  /// \param pTrackSameChargeFirst is the opposite charge track momentum
  /// \param ptD is the pt of the Ds meson candidate
  /// \return BIT(0) for KKpi, BIT(1) for piKK, BIT(2) for phipi, BIT(3) for piphi
  template <typename T>
  int8_t isSelectedDsInMassRange(const T& pTrackSameChargeFirst, const T& pTrackSameChargeSecond, const T& pTrackOppositeCharge, const float& ptD, int8_t isSelected)
  {
    int8_t retValue = 0;
    if (TESTBIT(isSelected, 0)) {
      auto invMassDsToKKPi = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackOppositeCharge, pTrackSameChargeSecond}, std::array{massK, massK, massPi});
      if (activateQA) {
        hMassVsPtC[kDs]->Fill(ptD, invMassDsToKKPi);
      }
      if (std::abs(invMassDsToKKPi - massDs) < deltaMassCharmHadronForBeauty) {
        retValue |= BIT(0);
      }
    }
    if (TESTBIT(isSelected, 1)) {
      auto invMassDsToPiKK = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackOppositeCharge, pTrackSameChargeSecond}, std::array{massPi, massK, massK});
      if (activateQA) {
        hMassVsPtC[kDs]->Fill(ptD, invMassDsToPiKK);
      }
      if (std::abs(invMassDsToPiKK - massDs) < deltaMassCharmHadronForBeauty) {
        retValue |= BIT(1);
      }
    }

    return retValue;
  }

  /// Mass selection of Lc candidates to build Lb candidates
  /// \param pTrackSameChargeFirst is the first same-charge track momentum
  /// \param pTrackSameChargeSecond is the second same-charge track momentum
  /// \param pTrackOppositeCharge is the opposite charge track momentum
  /// \param ptLc is the pt of the D0 meson candidate
  /// \return BIT(0) for pKpi with mass cut, BIT(1) for piKp with mass cut
  template <typename T>
  int8_t isSelectedLcInMassRange(const T& pTrackSameChargeFirst, const T& pTrackSameChargeSecond, const T& pTrackOppositeCharge, const float& ptLc, const int8_t isSelected)
  {
    int8_t retValue = 0;
    if (TESTBIT(isSelected, 0)) {
      auto invMassLcToPKPi = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackOppositeCharge, pTrackSameChargeSecond}, std::array{massProton, massK, massPi});
      if (activateQA) {
        hMassVsPtC[kLc]->Fill(ptLc, invMassLcToPKPi);
      }
      if (std::abs(invMassLcToPKPi - massLc) < deltaMassCharmHadronForBeauty) {
        retValue |= BIT(0);
      }
    }
    if (TESTBIT(isSelected, 1)) {
      auto invMassLcToPiKP = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackOppositeCharge, pTrackSameChargeSecond}, std::array{massPi, massK, massProton});
      if (activateQA) {
        hMassVsPtC[kLc]->Fill(ptLc, invMassLcToPiKP);
      }
      if (std::abs(invMassLcToPiKP - massLc) < deltaMassCharmHadronForBeauty) {
        retValue |= BIT(1);
      }
    }

    return retValue;
  }

  /// Mass selection of Xic candidates to build Lb candidates
  /// \param pTrackSameChargeFirst is the first same-charge track momentum
  /// \param pTrackSameChargeSecond is the second same-charge track momentum
  /// \param pTrackOppositeCharge is the opposite charge track momentum
  /// \param ptXic is the pt of the D0 meson candidate
  /// \return BIT(0) for pKpi with mass cut, BIT(1) for piKp with mass cut
  template <typename T>
  int8_t isSelectedXicInMassRange(const T& pTrackSameChargeFirst, const T& pTrackSameChargeSecond, const T& pTrackOppositeCharge, const float& ptXic, const int8_t isSelected)
  {
    int8_t retValue = 0;
    if (TESTBIT(isSelected, 0)) {
      auto invMassXicToPKPi = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackOppositeCharge, pTrackSameChargeSecond}, std::array{massProton, massK, massPi});
      if (activateQA) {
        hMassVsPtC[kXic]->Fill(ptXic, invMassXicToPKPi);
      }
      if (std::abs(invMassXicToPKPi - massXic) < deltaMassCharmHadronForBeauty) {
        retValue |= BIT(0);
      }
    }
    if (TESTBIT(isSelected, 1)) {
      auto invMassXicToPiKP = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackOppositeCharge, pTrackSameChargeSecond}, std::array{massPi, massK, massProton});
      if (activateQA) {
        hMassVsPtC[kXic]->Fill(ptXic, invMassXicToPiKP);
      }
      if (std::abs(invMassXicToPiKP - massXic) < deltaMassCharmHadronForBeauty) {
        retValue |= BIT(1);
      }
    }

    return retValue;
  }

  /// Single-track cuts for bachelor track of beauty candidates
  /// \param scores is a 3-element array with BDT out scores
  /// \return 0 if rejected, otherwise bitmap with BIT(RecoDecay::OriginType::Prompt) and/or BIT(RecoDecay::OriginType::NonPrompt) on
  template <typename T>
  int8_t isBDTSelected(const T& scores, const int candType)
  {
    int8_t retValue = 0;
    if (scores.size() < 3) {
      return retValue;
    }

    if (scores[0] > thresholdBDTScores[candType].get(0u, "BDTbkg")) {
      return retValue;
    }
    if (scores[1] > thresholdBDTScores[candType].get(0u, "BDTprompt")) {
      retValue |= BIT(RecoDecay::OriginType::Prompt);
    }
    if (scores[2] > thresholdBDTScores[candType].get(0u, "BDTnonprompt")) {
      retValue |= BIT(RecoDecay::OriginType::NonPrompt);
    }

    return retValue;
  }

  using HfTrackIndexProng2withColl = soa::Join<aod::Hf2Prongs, aod::Colls2Prong>;
  using HfTrackIndexProng3withColl = soa::Join<aod::Hf3Prongs, aod::Colls3Prong>;
  using BigTracksMCPID = soa::Join<aod::BigTracksExtended, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::BigTracksMC>;

  Filter trackFilter = requireGlobalTrackWoDCAInFilter();
  using BigTracksPID = soa::Filtered<soa::Join<aod::BigTracksExtended, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr>>;

  void process(aod::Collision const& collision,
               aod::BCsWithTimestamps const&,
               aod::V0Datas const& theV0s,
               HfTrackIndexProng2withColl const& cand2Prongs,
               HfTrackIndexProng3withColl const& cand3Prongs,
               BigTracksPID const& tracks)
  {
    if (applyOptimisation) {
      optimisationTreeCollisions(collision.globalIndex());
    }

    if (applyML && (loadModelsFromCCDB && timestampCCDB == 0) && !sessionML[kD0]) {
      for (auto iCharmPart{0}; iCharmPart < kNCharmParticles; ++iCharmPart) {
        if (onnxFiles[iCharmPart] != "") {
          sessionML[iCharmPart].reset(InitONNXSession(onnxFiles[iCharmPart], charmParticleNames[iCharmPart], envML[iCharmPart], sessionOptions[iCharmPart], inputShapesML[iCharmPart], dataTypeML[iCharmPart], loadModelsFromCCDB, ccdbApi, mlModelPathCCDB.value, timestampCCDB));
        }
      }
    }

    hProcessedEvents->Fill(0);

    // collision process loop
    bool keepEvent[kNtriggersHF]{false};
    //

    std::vector<std::vector<long>> indicesDau2Prong{};
    for (const auto& cand2Prong : cand2Prongs) {                                        // start loop over 2 prongs
      if (!TESTBIT(cand2Prong.hfflag(), o2::aod::hf_cand_2prong::DecayType::D0ToPiK)) { // check if it's a D0
        continue;
      }

      auto trackPos = cand2Prong.prong0_as<BigTracksPID>(); // positive daughter
      auto trackNeg = cand2Prong.prong1_as<BigTracksPID>(); // negative daughter

      auto preselD0 = isDzeroPreselected(trackPos, trackNeg);
      if (!preselD0) {
        continue;
      }

      std::array<float, 3> pVecPos = {trackPos.px(), trackPos.py(), trackPos.pz()};
      std::array<float, 3> pVecNeg = {trackNeg.px(), trackNeg.py(), trackNeg.pz()};

      bool isCharmTagged{true}, isBeautyTagged{true};

      // apply ML models
      int tagBDT = 0;
      float scoresToFill[3] = {-1., -1., -1.};
      if (applyML && onnxFiles[kD0] != "") {
        isCharmTagged = false;
        isBeautyTagged = false;

        // TODO: add more feature configurations
        std::vector<float> inputFeaturesD0{trackPos.pt(), trackPos.dcaXY(), trackPos.dcaZ(), trackNeg.pt(), trackNeg.dcaXY(), trackNeg.dcaZ()};
        std::vector<double> inputFeaturesDoD0{trackPos.pt(), trackPos.dcaXY(), trackPos.dcaZ(), trackNeg.pt(), trackNeg.dcaXY(), trackNeg.dcaZ()};

        if (dataTypeML[kD0] == 1) {
          auto scores = PredictONNX(inputFeaturesD0, sessionML[kD0], inputShapesML[kD0]);
          tagBDT = isBDTSelected(scores, kD0);
          for (int iScore{0}; iScore < 3; ++iScore) {
            scoresToFill[iScore] = scores[iScore];
          }
        } else if (dataTypeML[kD0] == 11) {
          auto scores = PredictONNX(inputFeaturesDoD0, sessionML[kD0], inputShapesML[kD0]);
          tagBDT = isBDTSelected(scores, kD0);
          for (int iScore{0}; iScore < 3; ++iScore) {
            scoresToFill[iScore] = scores[iScore];
          }
        } else {
          LOG(fatal) << "Error running model inference for D0: Unexpected input data type.";
        }

        if (applyML && activateQA) {
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
        optimisationTreeCharm(collision.globalIndex(), pdg::Code::kD0, pt2Prong, scoresToFill[0], scoresToFill[1], scoresToFill[2]);
      }

      auto selD0 = isSelectedD0InMassRange(pVecPos, pVecNeg, pt2Prong, preselD0);

      if (pt2Prong >= pTThreshold2Prong) {
        keepEvent[kHighPt2P] = true;
        if (activateQA) {
          hCharmHighPt[kD0]->Fill(pt2Prong);
        }
      } // end high-pT selection

      if (isCharmTagged) {
        indicesDau2Prong.push_back(std::vector<long>{trackPos.globalIndex(), trackNeg.globalIndex()});
      } // end multi-charm selection

      for (const auto& track : tracks) { // start loop over tracks
        if (track.globalIndex() == trackPos.globalIndex() || track.globalIndex() == trackNeg.globalIndex()) {
          continue;
        }

        std::array<float, 3> pVecThird = {track.px(), track.py(), track.pz()};

        if (!keepEvent[kBeauty3P] && isBeautyTagged) {
          int isTrackSelected = isSelectedTrackForBeauty(track, kBeauty3P);
          if (isTrackSelected && ((TESTBIT(selD0, 0) && track.signed1Pt() < 0) || (TESTBIT(selD0, 1) && track.signed1Pt() > 0))) {
            auto massCand = RecoDecay::m(std::array{pVec2Prong, pVecThird}, std::array{massD0, massPi});
            auto pVecBeauty3Prong = RecoDecay::pVec(pVec2Prong, pVecThird);
            auto ptCand = RecoDecay::pt(pVecBeauty3Prong);
            if (isTrackSelected == kRegular && std::abs(massCand - massBPlus) <= deltaMassBPlus) {
              keepEvent[kBeauty3P] = true;
              // fill optimisation tree for D0
              if (applyOptimisation) {
                optimisationTreeBeauty(collision.globalIndex(), pdg::Code::kD0, pt2Prong, scoresToFill[0], scoresToFill[1], scoresToFill[2], track.dcaXY());
              }
              if (activateQA) {
                hMassVsPtB[kBplus]->Fill(ptCand, massCand);
              }
            } else if (std::abs(massCand - massDStar) <= deltaMassDStar) { // additional check for B0->D*pi polarization studies
              if (activateQA) {
                hMassVsPtC[kNCharmParticles]->Fill(ptCand, massCand);
              }
              for (const auto& trackB : tracks) { // start loop over tracks
                if (track.signed1Pt() * trackB.signed1Pt() < 0 && isSelectedTrackForBeauty(trackB, kBeauty3P) == kRegular) {
                  std::array<float, 3> pVecFourth = {trackB.px(), trackB.py(), trackB.pz()};
                  auto massCandB0 = RecoDecay::m(std::array{pVec2Prong, pVecThird, pVecFourth}, std::array{massD0, massPi, massPi});
                  if (std::abs(massCandB0 - massB0) <= deltaMassB0) {
                    keepEvent[kBeauty3P] = true;
                    // fill optimisation tree for D0
                    if (applyOptimisation) {
                      optimisationTreeBeauty(collision.globalIndex(), 413, pt2Prong, scoresToFill[0], scoresToFill[1], scoresToFill[2], track.dcaXY()); // pdgCode of D*(2010)+: 413
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
        if (!keepEvent[kFemto2P] && isCharmTagged) {
          bool isProton = isSelectedProton4Femto(track);
          if (isProton) {
            float relativeMomentum = computeRelativeMomentum(track, pVec2Prong, massD0);
            if (applyOptimisation) {
              optimisationTreeFemto(collision.globalIndex(), pdg::Code::kD0, pt2Prong, scoresToFill[0], scoresToFill[1], scoresToFill[2], relativeMomentum, track.tpcNSigmaPr(), track.tofNSigmaPr());
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
      if (!keepEvent[kGammaCharm2P] && isCharmTagged) {
        for (auto& gamma : theV0s) {
          float V0CosinePA = gamma.v0cosPA(collision.posX(), collision.posY(), collision.posZ());
          bool isGamma = isSelectedGamma(gamma, V0CosinePA);
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
    for (const auto& cand3Prong : cand3Prongs) { // start loop over 3 prongs
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

      std::array<float, 3> pVecFirst = {trackFirst.px(), trackFirst.py(), trackFirst.pz()};
      std::array<float, 3> pVecSecond = {trackSecond.px(), trackSecond.py(), trackSecond.pz()};
      std::array<float, 3> pVecThird = {trackThird.px(), trackThird.py(), trackThird.pz()};

      if (is3Prong[0]) { // D+ preselections
        is3Prong[0] = isDplusPreselected(trackSecond);
      }
      if (is3Prong[1]) { // Ds preselections
        is3Prong[1] = isDsPreselected(pVecFirst, pVecThird, pVecSecond, trackSecond);
      }
      if (is3Prong[2] || is3Prong[3]) { // charm baryon preselections
        auto presel = isCharmBaryonPreselected(trackFirst, trackThird, trackSecond);
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
        std::vector<float> inputFeatures{trackFirst.pt(), trackFirst.dcaXY(), trackFirst.dcaZ(), trackSecond.pt(), trackSecond.dcaXY(), trackSecond.dcaZ(), trackThird.pt(), trackThird.dcaXY(), trackThird.dcaZ()};
        std::vector<double> inputFeaturesD{trackFirst.pt(), trackFirst.dcaXY(), trackFirst.dcaZ(), trackSecond.pt(), trackSecond.dcaXY(), trackSecond.dcaZ(), trackThird.pt(), trackThird.dcaXY(), trackThird.dcaZ()};
        for (auto iCharmPart{0}; iCharmPart < kNCharmParticles - 1; ++iCharmPart) {
          if (!is3Prong[iCharmPart] || onnxFiles[iCharmPart + 1] == "") {
            continue;
          }

          int tagBDT = 0;
          if (dataTypeML[iCharmPart + 1] == 1) {
            auto scores = PredictONNX(inputFeatures, sessionML[iCharmPart + 1], inputShapesML[iCharmPart + 1]);
            tagBDT = isBDTSelected(scores, iCharmPart + 1);
            for (int iScore{0}; iScore < 3; ++iScore) {
              scoresToFill[iCharmPart][iScore] = scores[iScore];
            }
          } else if (dataTypeML[iCharmPart + 1] == 11) {
            auto scores = PredictONNX(inputFeaturesD, sessionML[iCharmPart + 1], inputShapesML[iCharmPart + 1]);
            tagBDT = isBDTSelected(scores, iCharmPart + 1);
            for (int iScore{0}; iScore < 3; ++iScore) {
              scoresToFill[iCharmPart][iScore] = scores[iScore];
            }
          } else {
            LOG(error) << "Error running model inference for " << charmParticleNames[iCharmPart + 1].data() << ": Unexpected input data type.";
          }

          isCharmTagged[iCharmPart] = TESTBIT(tagBDT, RecoDecay::OriginType::Prompt);
          isBeautyTagged[iCharmPart] = TESTBIT(tagBDT, RecoDecay::OriginType::NonPrompt);

          if (activateQA) {
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
      float sign3Prong = trackFirst.signed1Pt() * trackSecond.signed1Pt() * trackThird.signed1Pt();

      std::array<int8_t, kNCharmParticles - 1> is3ProngInMass{0};
      if (is3Prong[0]) {
        is3ProngInMass[0] = isSelectedDplusInMassRange(pVecFirst, pVecThird, pVecSecond, pt3Prong);
        if (applyOptimisation) {
          optimisationTreeCharm(collision.globalIndex(), pdg::Code::kDPlus, pt3Prong, scoresToFill[0][0], scoresToFill[0][1], scoresToFill[0][2]);
        }
      }
      if (is3Prong[1]) {
        is3ProngInMass[1] = isSelectedDsInMassRange(pVecFirst, pVecThird, pVecSecond, pt3Prong, is3Prong[1]);
        if (applyOptimisation) {
          optimisationTreeCharm(collision.globalIndex(), pdg::Code::kDS, pt3Prong, scoresToFill[1][0], scoresToFill[1][1], scoresToFill[1][2]);
        }
      }
      if (is3Prong[2]) {
        is3ProngInMass[2] = isSelectedLcInMassRange(pVecFirst, pVecThird, pVecSecond, pt3Prong, is3Prong[2]);
        if (applyOptimisation) {
          optimisationTreeCharm(collision.globalIndex(), pdg::Code::kLambdaCPlus, pt3Prong, scoresToFill[2][0], scoresToFill[2][1], scoresToFill[2][2]);
        }
      }
      if (is3Prong[3]) {
        is3ProngInMass[3] = isSelectedXicInMassRange(pVecFirst, pVecThird, pVecSecond, pt3Prong, is3Prong[3]);
        if (applyOptimisation) {
          optimisationTreeCharm(collision.globalIndex(), pdg::Code::kXiCPlus, pt3Prong, scoresToFill[3][0], scoresToFill[3][1], scoresToFill[3][2]);
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

      for (const auto& track : tracks) { // start loop over tracks

        if (track.globalIndex() == trackFirst.globalIndex() || track.globalIndex() == trackSecond.globalIndex() || track.globalIndex() == trackThird.globalIndex()) {
          continue;
        }

        std::array<float, 3> pVecFourth = {track.px(), track.py(), track.pz()};

        int charmParticleID[kNBeautyParticles - 2] = {pdg::Code::kDPlus, pdg::Code::kDS, pdg::Code::kLambdaCPlus, pdg::Code::kXiCPlus};

        float massCharmHypos[kNBeautyParticles - 2] = {massDPlus, massDs, massLc, massXic};
        float massBeautyHypos[kNBeautyParticles - 2] = {massB0, massBs, massLb, massXib};
        float deltaMassHypos[kNBeautyParticles - 2] = {deltaMassB0, deltaMassBs, deltaMassLb, deltaMassXib};
        if (track.signed1Pt() * sign3Prong < 0 && isSelectedTrackForBeauty(track, kBeauty4P) == kRegular) {
          for (int iHypo{0}; iHypo < kNBeautyParticles - 2 && !keepEvent[kBeauty4P]; ++iHypo) {
            if (isBeautyTagged[iHypo] && (TESTBIT(is3ProngInMass[iHypo], 0) || TESTBIT(is3ProngInMass[iHypo], 1))) {
              auto massCandB = RecoDecay::m(std::array{pVec3Prong, pVecFourth}, std::array{massCharmHypos[iHypo], massPi});
              if (std::abs(massCandB - massBeautyHypos[iHypo]) <= deltaMassHypos[iHypo]) {
                keepEvent[kBeauty4P] = true;
                if (applyOptimisation) {
                  optimisationTreeBeauty(collision.globalIndex(), charmParticleID[iHypo], pt3Prong, scoresToFill[iHypo][0], scoresToFill[iHypo][1], scoresToFill[iHypo][2], track.dcaXY());
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
        if (isSelectedProton4Femto(track)) {
          for (int iHypo{0}; iHypo < kNCharmParticles - 1 && !keepEvent[kFemto3P]; ++iHypo) {
            if (isCharmTagged[iHypo]) {
              float relativeMomentum = computeRelativeMomentum(track, pVec3Prong, massCharmHypos[iHypo]);
              if (applyOptimisation) {
                optimisationTreeFemto(collision.globalIndex(), charmParticleID[iHypo], pt3Prong, scoresToFill[iHypo][0], scoresToFill[iHypo][1], scoresToFill[iHypo][2], relativeMomentum, track.tpcNSigmaPr(), track.tofNSigmaPr());
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
      if (!keepEvent[kGammaCharm3P] && isCharmTagged[kDs - 1]) {
        for (auto& gamma : theV0s) {
          float V0CosinePA = gamma.v0cosPA(collision.posX(), collision.posY(), collision.posZ());
          bool isGamma = isSelectedGamma(gamma, V0CosinePA);
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

  void
    processTraining(aod::Hf2Prongs const& cand2Prongs,
                    aod::Hf3Prongs const& cand3Prongs,
                    aod::McParticles const& particlesMC,
                    BigTracksMCPID const&)
  {
    for (const auto& cand2Prong : cand2Prongs) { // start loop over 2 prongs

      auto trackPos = cand2Prong.prong0_as<BigTracksMCPID>(); // positive daughter
      auto trackNeg = cand2Prong.prong1_as<BigTracksMCPID>(); // negative daughter

      std::array<float, 3> pVecPos = {trackPos.px(), trackPos.py(), trackPos.pz()};
      std::array<float, 3> pVecNeg = {trackNeg.px(), trackNeg.py(), trackNeg.pz()};
      auto pVec2Prong = RecoDecay::pVec(pVecPos, pVecNeg);
      auto pt2Prong = RecoDecay::pt(pVec2Prong);

      auto invMassD0 = RecoDecay::m(std::array{pVecPos, pVecNeg}, std::array{massPi, massK});
      auto invMassD0bar = RecoDecay::m(std::array{pVecPos, pVecNeg}, std::array{massK, massPi});

      int8_t sign = 0;
      int8_t flag = RecoDecay::OriginType::None;

      // D0(bar) → π± K∓
      auto indexRec = RecoDecay::getMatchedMCRec(particlesMC, std::array{trackPos, trackNeg}, pdg::Code::kD0, array{+kPiPlus, -kKPlus}, true, &sign);
      if (indexRec > -1) {
        auto particle = particlesMC.rawIteratorAt(indexRec);
        flag = RecoDecay::getCharmHadronOrigin(particlesMC, particle);
        if (flag < RecoDecay::OriginType::Prompt) {
          continue;
        }
      }

      double pseudoRndm = trackPos.pt() * 1000. - (long)(trackPos.pt() * 1000);
      if ((fillSignal && indexRec > -1) || (fillBackground && indexRec < 0 && pseudoRndm < donwSampleBkgFactor)) {
        train2P(invMassD0, invMassD0bar, pt2Prong, trackPos.pt(), trackPos.dcaXY(), trackPos.dcaZ(), trackPos.tpcNSigmaPi(), trackPos.tpcNSigmaKa(), trackPos.tofNSigmaPi(), trackPos.tofNSigmaKa(),
                trackNeg.pt(), trackNeg.dcaXY(), trackNeg.dcaZ(), trackNeg.tpcNSigmaPi(), trackNeg.tpcNSigmaKa(), trackNeg.tofNSigmaPi(), trackNeg.tofNSigmaKa(),
                flag);
      }
    } // end loop over 2-prong candidates

    for (const auto& cand3Prong : cand3Prongs) { // start loop over 3 prongs

      auto trackFirst = cand3Prong.prong0_as<BigTracksMCPID>();  // first daughter
      auto trackSecond = cand3Prong.prong1_as<BigTracksMCPID>(); // second daughter
      auto trackThird = cand3Prong.prong2_as<BigTracksMCPID>();  // third daughter
      auto arrayDaughters = std::array{trackFirst, trackSecond, trackThird};

      std::array<float, 3> pVecFirst = {trackFirst.px(), trackFirst.py(), trackFirst.pz()};
      std::array<float, 3> pVecSecond = {trackSecond.px(), trackSecond.py(), trackSecond.pz()};
      std::array<float, 3> pVecThird = {trackThird.px(), trackThird.py(), trackThird.pz()};

      auto pVec3Prong = RecoDecay::pVec(pVecFirst, pVecSecond, pVecThird);
      auto pt3Prong = RecoDecay::pt(pVec3Prong);

      auto invMassDplus = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massPi, massK, massPi});

      auto invMassDsToKKPi = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massK, massK, massPi});
      auto invMassDsToPiKK = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massPi, massK, massK});

      auto invMassLcToPKPi = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massProton, massK, massPi});
      auto invMassLcToPiKP = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massPi, massK, massProton});

      auto invMassXicToPKPi = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massProton, massK, massPi});
      auto invMassXicToPiKP = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massPi, massK, massProton});

      float deltaMassKKFirst = -1.f;
      float deltaMassKKSecond = -1.f;
      if (TESTBIT(cand3Prong.hfflag(), o2::aod::hf_cand_3prong::DecayType::DsToKKPi)) {
        deltaMassKKFirst = std::abs(RecoDecay::m(std::array{pVecFirst, pVecSecond}, std::array{massK, massK}) - massPhi);
        deltaMassKKSecond = std::abs(RecoDecay::m(std::array{pVecThird, pVecSecond}, std::array{massK, massK}) - massPhi);
      }
      int8_t sign = 0;
      int8_t flag = RecoDecay::OriginType::None;
      int8_t channel = -1;

      // D± → π± K∓ π±
      auto indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughters, pdg::Code::kDPlus, array{+kPiPlus, -kKPlus, +kPiPlus}, true, &sign, 2);
      if (indexRec >= 0) {
        channel = kDplus;
      }
      if (indexRec < 0) {
        // Ds± → K± K∓ π±
        indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughters, 431, array{+kKPlus, -kKPlus, +kPiPlus}, true, &sign, 2); // TODO: replace hard coded pdg code
        if (indexRec >= 0) {
          channel = kDs;
        }
      }
      if (indexRec < 0) {
        // Λc± → p± K∓ π±
        indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughters, pdg::Code::kLambdaCPlus, array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2);
        if (indexRec >= 0) {
          channel = kLc;
        }
      }
      if (indexRec < 0) {
        // Ξc± → p± K∓ π±
        indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughters, pdg::Code::kXiCPlus, array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2);
        if (indexRec >= 0) {
          channel = kXic;
        }
      }

      if (indexRec > -1) {
        auto particle = particlesMC.rawIteratorAt(indexRec);
        flag = RecoDecay::getCharmHadronOrigin(particlesMC, particle);
        if (flag < RecoDecay::OriginType::Prompt) {
          continue;
        }
      }

      double pseudoRndm = trackFirst.pt() * 1000. - (long)(trackFirst.pt() * 1000);
      if ((fillSignal && indexRec > -1) || (fillBackground && indexRec < 0 && pseudoRndm < donwSampleBkgFactor)) {
        train3P(invMassDplus, invMassDsToKKPi, invMassDsToPiKK, invMassLcToPKPi, invMassLcToPiKP, invMassXicToPKPi, invMassXicToPiKP, pt3Prong, deltaMassKKFirst, deltaMassKKSecond, trackFirst.pt(), trackFirst.dcaXY(), trackFirst.dcaZ(), trackFirst.tpcNSigmaPi(), trackFirst.tpcNSigmaKa(), trackFirst.tpcNSigmaPr(), trackFirst.tofNSigmaPi(), trackFirst.tofNSigmaKa(), trackFirst.tofNSigmaPr(),
                trackSecond.pt(), trackSecond.dcaXY(), trackSecond.dcaZ(), trackSecond.tpcNSigmaPi(), trackSecond.tpcNSigmaKa(), trackSecond.tpcNSigmaPr(), trackSecond.tofNSigmaPi(), trackSecond.tofNSigmaKa(), trackSecond.tofNSigmaPr(),
                trackThird.pt(), trackThird.dcaXY(), trackThird.dcaZ(), trackThird.tpcNSigmaPi(), trackThird.tpcNSigmaKa(), trackThird.tpcNSigmaPr(), trackThird.tofNSigmaPi(), trackThird.tofNSigmaKa(), trackThird.tofNSigmaPr(),
                flag, channel, cand3Prong.hfflag());
      }
    } // end loop over 3-prong candidates
  }

  PROCESS_SWITCH(HfFilter, processTraining, "Process MC for training samples", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{

  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<AddCollisionId>(cfg));
  workflow.push_back(adaptAnalysisTask<HfFilter>(cfg));

  return workflow;
}
