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

/// \file HFFilterCharmHadronSignals.cxx
/// \brief task for the quality control of the signals of D0, D+, Ds+, Lc+, and D*+ selected in the HFFilter.cxx task
///
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

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
#include "Common/DataModel/Multiplicity.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "EventFiltering/filterTables.h"
#include "EventFiltering/PWGHF/HFFilterHelpers.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::aod::hffilters;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct HfFilterCharmHadronSignals { // Main struct for HF triggers

  Configurable<bool> applyEventSelection{"applyEventSelection", false, "flag to enable event selection (sel8 + Zvt)"};

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

  // additional selections for D*
  Configurable<float> minPtSoftPion{"minPtSoftPion", static_cast<float>(cutsPt[0][1]), "minimum pT for soft pion tracks in D*+ -> D0pi decay"};
  Configurable<float> maxPtSoftPion{"maxPtSoftPion", static_cast<float>(cutsPt[1][1]), "maximum pT for soft pion tracks in D*+ -> D0pi decay"};
  Configurable<float> maxDeltaMassDstar{"maxDeltaMassDstar", static_cast<float>(cutsMassCharmReso[0][0]), "maximum invariant-mass delta for D*+ in GeV/c2"};
  Configurable<float> maxDeltaMassDzeroFromDstar{"maxDeltaMassDzeroFromDstar", static_cast<float>(cutsDeltaMassB[0][kNBeautyParticles]), "maximum invariant-mass delta for D0 in D*+ -> D0pi decay"};

  // CCDB configuration
  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> mlModelPathCCDB{"mlModelPathCCDB", "Analysis/PWGHF/ML/HFTrigger/", "Path on CCDB"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB. Exceptions: > 0 for the specific timestamp, 0 gets the run dependent timestamp"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};
  int currentRun{0}; // needed to get proper magnetic field

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
  std::array<std::string, kNCharmParticles> onnxFiles;
  std::array<LabeledArray<double>, kNCharmParticles> thresholdBDTScores;

  ConfigurableAxis pvContributorsAxis{"pvContributorsAxis", {250, 0.f, 250.f}, "PV contributors"};
  ConfigurableAxis multiplicityAxis{"multiplicityAxis", {100, 0.f, 100.f}, "MultFT0M"};
  ConfigurableAxis zVtxAxis{"zVtxAxis", {150, -15.f, 15.f}, "#it{z}_{vtx} (cm)"};
  ConfigurableAxis invMassDmesAxis = {"invMassDmesAxis", {300, 1.65f, 2.25f}, "inv. mass (GeV/#it{c}^{2})"};
  ConfigurableAxis invMassDstarAxis = {"invMassDstarAxis", {180, 0.f, 0.18f}, "inv. mass difference (GeV/#it{c}^{2})"};
  ConfigurableAxis invMassCbaryonAxis = {"invMassCbaryonAxis", {300, 2.05f, 2.65f}, "inv. mass (GeV/#it{c}^{2})"};
  ConfigurableAxis ptAxis = {"ptAxis", {100, 0.f, 50.f}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis yAxis = {"yAxis", {10, -1.f, 1.f}, "#it{y}"};
  ConfigurableAxis phiAxis = {"phiAxis", {90, -constants::math::PI, constants::math::PI}, "#varphi (rad)"};
  ConfigurableAxis bdtPromptAxis{"bdtPromptAxis", {100, 0.f, 1.f}, "BDT prompt"};
  ConfigurableAxis bdtNonPromptAxis{"bdtNonPromptAxis", {100, 0.f, 1.f}, "BDT nonprompt"};

  HistogramRegistry registry{"registry", {{"hCollisions", "", {HistType::kTH3F, {zVtxAxis, pvContributorsAxis, multiplicityAxis}}}, {"hDzeroToKPi", "", {HistType::kTHnSparseF, {invMassDmesAxis, ptAxis, yAxis, phiAxis, zVtxAxis, pvContributorsAxis, multiplicityAxis, bdtPromptAxis, bdtNonPromptAxis}}}, {"hDplusToKPiPi", "", {HistType::kTHnSparseF, {invMassDmesAxis, ptAxis, yAxis, phiAxis, zVtxAxis, pvContributorsAxis, multiplicityAxis, bdtPromptAxis, bdtNonPromptAxis}}}, {"hDsToKKPi", "", {HistType::kTHnSparseF, {invMassDmesAxis, ptAxis, yAxis, phiAxis, zVtxAxis, pvContributorsAxis, multiplicityAxis, bdtPromptAxis, bdtNonPromptAxis}}}, {"hDstarToDzeroPi", "", {HistType::kTHnSparseF, {invMassDstarAxis, ptAxis, yAxis, phiAxis, zVtxAxis, pvContributorsAxis, multiplicityAxis, bdtPromptAxis, bdtNonPromptAxis}}}, {"hDstarToDzeroPiForBeauty", "", {HistType::kTHnSparseF, {invMassDstarAxis, ptAxis, yAxis, phiAxis, zVtxAxis, pvContributorsAxis, multiplicityAxis, bdtPromptAxis, bdtNonPromptAxis}}}, {"hLcToPKPi", "", {HistType::kTHnSparseF, {invMassCbaryonAxis, ptAxis, yAxis, phiAxis, zVtxAxis, pvContributorsAxis, multiplicityAxis, bdtPromptAxis, bdtNonPromptAxis}}}, {"hXicPlusToPKPi", "", {HistType::kTHnSparseF, {invMassCbaryonAxis, ptAxis, yAxis, phiAxis, zVtxAxis, pvContributorsAxis, multiplicityAxis, bdtPromptAxis, bdtNonPromptAxis}}}}};

  // no material correction for track propagation
  o2::base::Propagator::MatCorrType noMatCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  // helper object
  HfFilterHelper helper;

  void init(InitContext&)
  {
    helper.setPtLimitsDstarSoftPion(minPtSoftPion, maxPtSoftPion);
    helper.setPtBinsSingleTracks(std::vector<double>{hf_cuts_single_track::vecBinsPtTrack});
    helper.setCutsSingleTrackBeauty(cutsSingleTrackDummy, cutsSingleTrackDummy);
    helper.setDeltaMassCharmHadForBeauty(maxDeltaMassDzeroFromDstar);

    ccdb->setURL(url.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());

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
  }

  using CollsWithEvSelAndMult = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>;
  using BigTracksPID = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr>;

  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;
  Preslice<aod::Hf2Prongs> hf2ProngPerCollision = aod::track_association::collisionId;
  Preslice<aod::Hf3Prongs> hf3ProngPerCollision = aod::track_association::collisionId;

  void process(CollsWithEvSelAndMult const& collisions,
               aod::BCsWithTimestamps const&,
               aod::Hf2Prongs const& cand2Prongs,
               aod::Hf3Prongs const& cand3Prongs,
               aod::TrackAssoc const& trackIndices,
               BigTracksPID const& tracks)
  {
    for (const auto& collision : collisions) {
      if (applyEventSelection && (!collision.sel8() || std::fabs(collision.posZ()) > 11.f)) { // safety margin for Zvtx
        continue;
      }

      registry.fill(HIST("hCollisions"), collision.posZ(), collision.numContrib(), collision.multFT0M());

      auto thisCollId = collision.globalIndex();
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
        currentRun = bc.runNumber();
      }

      // loop over 2-prong candidates
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

        // apply ML models
        int tagBDT = 0;
        float scoresToFill[3] = {-1., -1., -1.};
        if (applyML && onnxFiles[kD0] != "") {
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

          if (!TESTBIT(tagBDT, RecoDecay::OriginType::Prompt) && !TESTBIT(tagBDT, RecoDecay::OriginType::NonPrompt)) { // if not tagged neither as prompt nor nonprompt, we skip
            continue;
          }
        }

        auto pVec2Prong = RecoDecay::pVec(pVecPos, pVecNeg);
        auto pt2Prong = RecoDecay::pt(pVec2Prong);
        auto phi2Prong = RecoDecay::phi(pVec2Prong[0], pVec2Prong[1]);
        auto yD0 = RecoDecay::y(pVec2Prong, massD0);
        auto invMassD0 = RecoDecay::m(std::array{pVecPos, pVecNeg}, std::array{massPi, massKa});
        auto invMassD0bar = RecoDecay::m(std::array{pVecPos, pVecNeg}, std::array{massKa, massPi});
        // fill THnSparse
        if (TESTBIT(preselD0, 0)) { // D0
          registry.fill(HIST("hDzeroToKPi"), invMassD0, pt2Prong, yD0, phi2Prong, collision.posZ(), collision.numContrib(), collision.multFT0M(), scoresToFill[1], scoresToFill[2]);
        }
        if (TESTBIT(preselD0, 1)) { // D0bar
          registry.fill(HIST("hDzeroToKPi"), invMassD0bar, pt2Prong, yD0, phi2Prong, collision.posZ(), collision.numContrib(), collision.multFT0M(), scoresToFill[1], scoresToFill[2]);
        }

        // we build D* here, as done in the HFFilter.cxx
        TH2* histNullptr = nullptr;
        auto selD0 = helper.isSelectedD0InMassRange(pVecPos, pVecNeg, pt2Prong, preselD0, false, histNullptr);
        auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
        for (const auto& trackId : trackIdsThisCollision) { // start loop over tracks
          auto track = trackId.track_as<BigTracksPID>();

          if (track.globalIndex() == trackPos.globalIndex() || track.globalIndex() == trackNeg.globalIndex()) {
            continue;
          }
          if ((track.sign() > 0 && !TESTBIT(selD0, 0)) || (track.sign() < 0 && !TESTBIT(selD0, 1))) {
            continue;
          }

          auto trackParThird = getTrackPar(track);
          o2::gpu::gpustd::array<float, 2> dcaThird{track.dcaXY(), track.dcaZ()};
          std::array<float, 3> pVecThird = {track.px(), track.py(), track.pz()};
          if (track.collisionId() != thisCollId) {
            o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParThird, 2.f, noMatCorr, &dcaThird);
            getPxPyPz(trackParThird, pVecThird);
          }
          auto isTrackSelected = helper.isSelectedTrackForSoftPionOrBeauty(track, trackParThird, dcaThird, 2);
          if (TESTBIT(isTrackSelected, kSoftPion)) {
            std::array<float, 2> massDausD0{massPi, massKa};
            auto invMassD0dau = invMassD0;
            if (track.sign() < 0) {
              massDausD0[0] = massKa;
              massDausD0[1] = massPi;
              invMassD0dau = invMassD0bar;
            }
            auto invMassDstar = RecoDecay::m(std::array{pVecPos, pVecNeg, pVecThird}, std::array{massDausD0[0], massDausD0[1], massPi});
            auto massDiffDstar = invMassDstar - invMassD0dau;
            auto pVecDstar = RecoDecay::pVec(pVec2Prong, pVecThird);
            auto ptDstar = RecoDecay::pt(pVecDstar);
            auto phiDstar = RecoDecay::phi(pVecDstar[0], pVecDstar[1]);
            auto yDstar = RecoDecay::y(pVecDstar, massDStar);
            if (std::fabs(massDiffDstar - (massDStar - massD0)) <= maxDeltaMassDstar) {
              registry.fill(HIST("hDstarToDzeroPi"), massDiffDstar, ptDstar, yDstar, phiDstar, collision.posZ(), collision.numContrib(), collision.multFT0M(), scoresToFill[1], scoresToFill[2]); // for D* we store the BDT output scores of the D0
              if (TESTBIT(isTrackSelected, kSoftPionForBeauty)) {
                registry.fill(HIST("hDstarToDzeroPiForBeauty"), massDiffDstar, ptDstar, yDstar, phiDstar, collision.posZ(), collision.numContrib(), collision.multFT0M(), scoresToFill[1], scoresToFill[2]); // for D* we store the BDT output scores of the D0
              }
            }
          }
        }
      } // end 2-prong loop

      // loop over 3-prong candidates
      auto cand3ProngsThisColl = cand3Prongs.sliceBy(hf3ProngPerCollision, thisCollId);
      for (const auto& cand3Prong : cand3ProngsThisColl) {
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
        if (applyML) {
          isCharmTagged = std::array<int8_t, kNCharmParticles - 1>{0};
          isBeautyTagged = std::array<int8_t, kNCharmParticles - 1>{0};

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
          }
        }

        if (!std::accumulate(isCharmTagged.begin(), isCharmTagged.end(), 0) && !std::accumulate(isBeautyTagged.begin(), isBeautyTagged.end(), 0)) {
          continue;
        }

        auto pVec3Prong = RecoDecay::pVec(pVecFirst, pVecSecond, pVecThird);
        auto pt3Prong = RecoDecay::pt(pVec3Prong);
        auto phi3Prong = RecoDecay::phi(pVec3Prong[0], pVec3Prong[1]);
        if (is3Prong[0]) { // D+
          auto yDplus = RecoDecay::y(pVec3Prong, massDPlus);
          auto invMassDplus = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massPi, massKa, massPi});
          registry.fill(HIST("hDplusToKPiPi"), invMassDplus, pt3Prong, yDplus, phi3Prong, collision.posZ(), collision.numContrib(), collision.multFT0M(), scoresToFill[0][1], scoresToFill[0][2]);
        }
        if (is3Prong[1]) { // Ds+
          auto yDs = RecoDecay::y(pVec3Prong, massDs);
          if (TESTBIT(is3Prong[1], 0)) {
            auto invMassDsToKKPi = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massKa, massKa, massPi});
            registry.fill(HIST("hDsToKKPi"), invMassDsToKKPi, pt3Prong, yDs, phi3Prong, collision.posZ(), collision.numContrib(), collision.multFT0M(), scoresToFill[1][1], scoresToFill[1][2]);
          }
          if (TESTBIT(is3Prong[1], 1)) {
            auto invMassDsToPiKK = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massPi, massKa, massKa});
            registry.fill(HIST("hDsToKKPi"), invMassDsToPiKK, pt3Prong, yDs, phi3Prong, collision.posZ(), collision.numContrib(), collision.multFT0M(), scoresToFill[1][1], scoresToFill[1][2]);
          }
        }
        if (is3Prong[2]) { // Lc+
          auto yLc = RecoDecay::y(pVec3Prong, massLc);
          if (TESTBIT(is3Prong[2], 0)) {
            auto invMassLcToPKPi = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massProton, massKa, massPi});
            registry.fill(HIST("hLcToPKPi"), invMassLcToPKPi, pt3Prong, yLc, phi3Prong, collision.posZ(), collision.numContrib(), collision.multFT0M(), scoresToFill[2][1], scoresToFill[2][2]);
          }
          if (TESTBIT(is3Prong[2], 1)) {
            auto invMassLcToPiKP = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massPi, massKa, massProton});
            registry.fill(HIST("hLcToPKPi"), invMassLcToPiKP, pt3Prong, yLc, phi3Prong, collision.posZ(), collision.numContrib(), collision.multFT0M(), scoresToFill[2][1], scoresToFill[2][2]);
          }
        }
      } // end 3-prong loop
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{

  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<HfFilterCharmHadronSignals>(cfg));

  return workflow;
}
