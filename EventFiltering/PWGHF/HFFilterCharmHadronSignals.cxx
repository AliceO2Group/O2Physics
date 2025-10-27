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

/// \file HFFilterCharmHadronSignals.cxx
/// \brief task for the quality control of the signals of D0, D+, Ds+, Lc+, and D*+ selected in the HFFilter.cxx task
///
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

#include "EventFiltering/PWGHF/HFFilterHelpers.h"
//
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"
//
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TH2.h>

#include <Rtypes.h>

#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::analysis;
using namespace o2::aod::hffilters;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct HfFilterCharmHadronSignals { // Main struct for HF triggers

  Configurable<bool> applyEventSelection{"applyEventSelection", true, "flag to enable event selection (sel8 + Zvt)"};
  Configurable<bool> applyTimeFrameBorderCut{"applyTimeFrameBorderCut", true, "flag to enable time-frame border cut"};

  // parameters for ML application
  Configurable<std::vector<double>> pTBinsBDT{"pTBinsBDT", std::vector<double>{hf_cuts_bdt_multiclass::vecBinsPt}, "track pT bin limits for BDT cut"};

  Configurable<LabeledArray<double>> thresholdBDTScoreD0ToKPi{"thresholdBDTScoreD0ToKPi", {hf_cuts_bdt_multiclass::Cuts[0], hf_cuts_bdt_multiclass::NBinsPt, hf_cuts_bdt_multiclass::NCutBdtScores, hf_cuts_bdt_multiclass::labelsPt, hf_cuts_bdt_multiclass::labelsCutBdt}, "Threshold values for BDT output scores of D0 candidates"};
  Configurable<LabeledArray<double>> thresholdBDTScoreDPlusToPiKPi{"thresholdBDTScoreDPlusToPiKPi", {hf_cuts_bdt_multiclass::Cuts[0], hf_cuts_bdt_multiclass::NBinsPt, hf_cuts_bdt_multiclass::NCutBdtScores, hf_cuts_bdt_multiclass::labelsPt, hf_cuts_bdt_multiclass::labelsCutBdt}, "Threshold values for BDT output scores of D+ candidates"};
  Configurable<LabeledArray<double>> thresholdBDTScoreDSToPiKK{"thresholdBDTScoreDSToPiKK", {hf_cuts_bdt_multiclass::Cuts[0], hf_cuts_bdt_multiclass::NBinsPt, hf_cuts_bdt_multiclass::NCutBdtScores, hf_cuts_bdt_multiclass::labelsPt, hf_cuts_bdt_multiclass::labelsCutBdt}, "Threshold values for BDT output scores of Ds+ candidates"};
  Configurable<LabeledArray<double>> thresholdBDTScoreLcToPiKP{"thresholdBDTScoreLcToPiKP", {hf_cuts_bdt_multiclass::Cuts[0], hf_cuts_bdt_multiclass::NBinsPt, hf_cuts_bdt_multiclass::NCutBdtScores, hf_cuts_bdt_multiclass::labelsPt, hf_cuts_bdt_multiclass::labelsCutBdt}, "Threshold values for BDT output scores of Lc+ candidates"};
  Configurable<LabeledArray<double>> thresholdBDTScoreXicToPiKP{"thresholdBDTScoreXicToPiKP", {hf_cuts_bdt_multiclass::Cuts[0], hf_cuts_bdt_multiclass::NBinsPt, hf_cuts_bdt_multiclass::NCutBdtScores, hf_cuts_bdt_multiclass::labelsPt, hf_cuts_bdt_multiclass::labelsCutBdt}, "Threshold values for BDT output scores of Xic+ candidates"};
  Configurable<std::string> paramCharmMassShape{"paramCharmMassShape", "2023_pass3", "Parametrisation of charm-hadron mass shape (options: 2023_pass3)"};
  Configurable<float> numSigmaDeltaMassCharmHad{"numSigmaDeltaMassCharmHad", 2.5, "Number of sigma for charm-hadron delta mass cut in B and D resonance triggers"};

  // additional selections for D*
  Configurable<float> minPtSoftPion{"minPtSoftPion", static_cast<float>(cutsPt[0][1]), "minimum pT for soft pion tracks in D*+ -> D0pi decay"};
  Configurable<float> maxPtSoftPion{"maxPtSoftPion", static_cast<float>(cutsPt[1][1]), "maximum pT for soft pion tracks in D*+ -> D0pi decay"};
  Configurable<float> minDeltaMassDstar{"minDeltaMassDstar", static_cast<float>(cutsCharmReso[0][0]), "minimum invariant-mass delta for D*+ in GeV/c2"};
  Configurable<float> maxDeltaMassDstar{"maxDeltaMassDstar", static_cast<float>(cutsCharmReso[1][0]), "maximum invariant-mass delta for D*+ in GeV/c2"};
  Configurable<std::vector<double>> pTBinsTrack{"pTBinsTrack", std::vector<double>{hf_cuts_single_track::vecBinsPtTrack}, "track pT bin limits for DCAXY pT-dependent cut (D* from beauty)"};
  Configurable<LabeledArray<double>> cutsTrackBeauty3Prong{"cutsTrackBeauty3Prong", {hf_cuts_single_track::CutsTrack[0], hf_cuts_single_track::NBinsPtTrack, hf_cuts_single_track::NCutVarsTrack, hf_cuts_single_track::labelsPtTrack, hf_cuts_single_track::labelsCutVarTrack}, "Single-track selections per pT bin for 3-prong beauty candidates"};

  // CCDB configuration
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  int currentRun{0}; // needed to detect if the run changed and trigger update of calibrations etc.

  // array of BDT thresholds
  std::array<LabeledArray<double>, kNCharmParticles> thresholdBDTScores;

  ConfigurableAxis pvContributorsAxis{"pvContributorsAxis", {250, 0.f, 250.f}, "PV contributors"};
  ConfigurableAxis multiplicityAxis{"multiplicityAxis", {100, 0.f, 1000.f}, "MultFT0M"};
  ConfigurableAxis zVtxAxis{"zVtxAxis", {150, -15.f, 15.f}, "#it{z}_{vtx} (cm)"};
  ConfigurableAxis invMassDmesAxis = {"invMassDmesAxis", {300, 1.65f, 2.25f}, "inv. mass (GeV/#it{c}^{2})"};
  ConfigurableAxis invMassDstarAxis = {"invMassDstarAxis", {180, 0.f, 0.18f}, "inv. mass difference (GeV/#it{c}^{2})"};
  ConfigurableAxis invMassCbaryonAxis = {"invMassCbaryonAxis", {300, 2.05f, 2.65f}, "inv. mass (GeV/#it{c}^{2})"};
  ConfigurableAxis ptAxis = {"ptAxis", {100, 0.f, 50.f}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis yAxis = {"yAxis", {10, -1.f, 1.f}, "#it{y}"};
  ConfigurableAxis phiAxis = {"phiAxis", {90, 0., constants::math::TwoPI}, "#varphi (rad)"};
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
    helper.setPtBinsSingleTracks(pTBinsTrack);
    helper.setCutsSingleTrackBeauty(cutsTrackBeauty3Prong, cutsTrackBeauty3Prong, cutsTrackBeauty3Prong);
    helper.setMassResolParametrisation(paramCharmMassShape);
    helper.setNumSigmaForDeltaMassCharmHadCut(numSigmaDeltaMassCharmHad);

    ccdb->setURL(url.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());

    thresholdBDTScores = {thresholdBDTScoreD0ToKPi, thresholdBDTScoreDPlusToPiKPi, thresholdBDTScoreDSToPiKK, thresholdBDTScoreLcToPiKP, thresholdBDTScoreXicToPiKP};
  }

  using CollsWithEvSelAndMult = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>;
  using BigTracksPID = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr>;

  using Hf2ProngsWithMl = soa::Join<aod::Hf2Prongs, aod::Hf2ProngMlProbs>;
  using Hf3ProngsWithMl = soa::Join<aod::Hf3Prongs, aod::Hf3ProngMlProbs>;

  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;
  Preslice<Hf2ProngsWithMl> hf2ProngPerCollision = aod::track_association::collisionId;
  Preslice<Hf3ProngsWithMl> hf3ProngPerCollision = aod::track_association::collisionId;

  void process(CollsWithEvSelAndMult const& collisions,
               aod::BCsWithTimestamps const&,
               Hf2ProngsWithMl const& cand2Prongs,
               Hf3ProngsWithMl const& cand3Prongs,
               aod::TrackAssoc const& trackIndices,
               BigTracksPID const& /*tracks*/)
  {
    for (const auto& collision : collisions) {
      if (applyEventSelection && (!collision.sel8() || std::fabs(collision.posZ()) > 11.f || (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder) && applyTimeFrameBorderCut))) { // safety margin for Zvtx
        continue;
      }

      registry.fill(HIST("hCollisions"), collision.posZ(), collision.numContrib(), collision.multFT0M());

      auto thisCollId = collision.globalIndex();
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
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
        std::array<float, 2> dcaPos{trackPos.dcaXY(), trackPos.dcaZ()};
        std::array<float, 2> dcaNeg{trackNeg.dcaXY(), trackNeg.dcaZ()};
        std::array<float, 3> pVecPos{trackPos.pVector()};
        std::array<float, 3> pVecNeg{trackNeg.pVector()};
        if (trackPos.collisionId() != thisCollId) {
          o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParPos, 2.f, noMatCorr, &dcaPos);
          getPxPyPz(trackParPos, pVecPos);
        }
        if (trackNeg.collisionId() != thisCollId) {
          o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParNeg, 2.f, noMatCorr, &dcaNeg);
          getPxPyPz(trackParNeg, pVecNeg);
        }

        // apply ML models
        std::vector<float> scores{};
        scores.insert(scores.end(), cand2Prong.mlProbSkimD0ToKPi().begin(), cand2Prong.mlProbSkimD0ToKPi().end());
        if (scores.size() != 3) {
          scores.resize(3);
          scores[0] = 2.;
          scores[1] = -1.;
          scores[2] = -1.;
        }
        auto tagBDT = helper.isBDTSelected(scores, thresholdBDTScores[kD0]);
        if (!TESTBIT(tagBDT, RecoDecay::OriginType::None)) {
          continue;
        }

        auto pVec2Prong = RecoDecay::pVec(pVecPos, pVecNeg);
        auto pt2Prong = RecoDecay::pt(pVec2Prong);
        auto phi2Prong = RecoDecay::phi(pVec2Prong[0], pVec2Prong[1]);
        auto yD0 = RecoDecay::y(pVec2Prong, massD0);
        auto invMassD0 = RecoDecay::m(std::array{pVecPos, pVecNeg}, std::array{massPi, massKa});
        auto invMassD0bar = RecoDecay::m(std::array{pVecPos, pVecNeg}, std::array{massKa, massPi});
        // fill THnSparse
        if (TESTBIT(preselD0, 0)) { // D0
          registry.fill(HIST("hDzeroToKPi"), invMassD0, pt2Prong, yD0, phi2Prong, collision.posZ(), collision.numContrib(), collision.multFT0M(), scores[1], scores[2]);
        }
        if (TESTBIT(preselD0, 1)) { // D0bar
          registry.fill(HIST("hDzeroToKPi"), invMassD0bar, pt2Prong, yD0, phi2Prong, collision.posZ(), collision.numContrib(), collision.multFT0M(), scores[1], scores[2]);
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
          std::array<float, 2> dcaThird{track.dcaXY(), track.dcaZ()};
          std::array<float, 3> pVecThird = track.pVector();
          if (track.collisionId() != thisCollId) {
            o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParThird, 2.f, noMatCorr, &dcaThird);
            getPxPyPz(trackParThird, pVecThird);
          }
          auto isTrackSelected = helper.isSelectedTrackForSoftPionOrBeauty<kBeauty3P>(track, trackParThird, dcaThird);
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
            if (minDeltaMassDstar <= massDiffDstar && massDiffDstar <= maxDeltaMassDstar) {
              registry.fill(HIST("hDstarToDzeroPi"), massDiffDstar, ptDstar, yDstar, phiDstar, collision.posZ(), collision.numContrib(), collision.multFT0M(), scores[1], scores[2]); // for D* we store the BDT output scores of the D0
              if (TESTBIT(isTrackSelected, kSoftPionForBeauty)) {
                registry.fill(HIST("hDstarToDzeroPiForBeauty"), massDiffDstar, ptDstar, yDstar, phiDstar, collision.posZ(), collision.numContrib(), collision.multFT0M(), scores[1], scores[2]); // for D* we store the BDT output scores of the D0
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
        std::array<float, 2> dcaFirst{trackFirst.dcaXY(), trackFirst.dcaZ()};
        std::array<float, 2> dcaSecond{trackSecond.dcaXY(), trackSecond.dcaZ()};
        std::array<float, 2> dcaThird{trackThird.dcaXY(), trackThird.dcaZ()};
        std::array<float, 3> pVecFirst = trackFirst.pVector();
        std::array<float, 3> pVecSecond = trackSecond.pVector();
        std::array<float, 3> pVecThird = trackThird.pVector();
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

        std::array<int8_t, kNCharmParticles - 1> isSignalTagged = is3Prong;

        std::array<std::vector<float>, kNCharmParticles - 1> scores{};
        scores[0].insert(scores[0].end(), cand3Prong.mlProbSkimDplusToPiKPi().begin(), cand3Prong.mlProbSkimDplusToPiKPi().end());
        scores[1].insert(scores[1].end(), cand3Prong.mlProbSkimDsToKKPi().begin(), cand3Prong.mlProbSkimDsToKKPi().end());
        scores[2].insert(scores[2].end(), cand3Prong.mlProbSkimLcToPKPi().begin(), cand3Prong.mlProbSkimLcToPKPi().end());
        scores[3].insert(scores[3].end(), cand3Prong.mlProbSkimXicToPKPi().begin(), cand3Prong.mlProbSkimXicToPKPi().end());

        for (auto iCharmPart{0}; iCharmPart < kNCharmParticles - 1; ++iCharmPart) {
          if (!is3Prong[iCharmPart]) { // we immediately skip if it was not selected for a given 3-prong species
            continue;
          }

          if (scores[iCharmPart].size() != 3) {
            scores[iCharmPart].resize(3);
            scores[iCharmPart][0] = 2.;
            scores[iCharmPart][1] = -1.;
            scores[iCharmPart][2] = -1.;
          }
          auto tagBDT = helper.isBDTSelected(scores[iCharmPart], thresholdBDTScores[iCharmPart + 1]);

          isSignalTagged[iCharmPart] = TESTBIT(tagBDT, RecoDecay::OriginType::None);
        }

        if (!std::accumulate(isSignalTagged.begin(), isSignalTagged.end(), 0)) {
          continue;
        }

        auto pVec3Prong = RecoDecay::pVec(pVecFirst, pVecSecond, pVecThird);
        auto pt3Prong = RecoDecay::pt(pVec3Prong);
        auto phi3Prong = RecoDecay::phi(pVec3Prong[0], pVec3Prong[1]);
        if (is3Prong[0]) { // D+
          auto yDplus = RecoDecay::y(pVec3Prong, massDPlus);
          auto invMassDplus = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massPi, massKa, massPi});
          registry.fill(HIST("hDplusToKPiPi"), invMassDplus, pt3Prong, yDplus, phi3Prong, collision.posZ(), collision.numContrib(), collision.multFT0M(), scores[0][1], scores[0][2]);
        }
        if (is3Prong[1]) { // Ds+
          auto yDs = RecoDecay::y(pVec3Prong, massDs);
          if (TESTBIT(is3Prong[1], 0)) {
            auto invMassDsToKKPi = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massKa, massKa, massPi});
            registry.fill(HIST("hDsToKKPi"), invMassDsToKKPi, pt3Prong, yDs, phi3Prong, collision.posZ(), collision.numContrib(), collision.multFT0M(), scores[1][1], scores[1][2]);
          }
          if (TESTBIT(is3Prong[1], 1)) {
            auto invMassDsToPiKK = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massPi, massKa, massKa});
            registry.fill(HIST("hDsToKKPi"), invMassDsToPiKK, pt3Prong, yDs, phi3Prong, collision.posZ(), collision.numContrib(), collision.multFT0M(), scores[1][1], scores[1][2]);
          }
        }
        if (is3Prong[2]) { // Lc+
          auto yLc = RecoDecay::y(pVec3Prong, massLc);
          if (TESTBIT(is3Prong[2], 0)) {
            auto invMassLcToPKPi = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massProton, massKa, massPi});
            registry.fill(HIST("hLcToPKPi"), invMassLcToPKPi, pt3Prong, yLc, phi3Prong, collision.posZ(), collision.numContrib(), collision.multFT0M(), scores[2][1], scores[2][2]);
          }
          if (TESTBIT(is3Prong[2], 1)) {
            auto invMassLcToPiKP = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massPi, massKa, massProton});
            registry.fill(HIST("hLcToPKPi"), invMassLcToPiKP, pt3Prong, yLc, phi3Prong, collision.posZ(), collision.numContrib(), collision.multFT0M(), scores[2][1], scores[2][2]);
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
