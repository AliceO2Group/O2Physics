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

/// \file derivedDataCreatorCorrelationsReduced.cxx
/// \brief CharmHadrons-Hadrons correlator tree creator for data and MC-reco analyses
/// \author Marcello Di Costanzo <marcello.di.costanzo@cern.ch>, Politecnico and INFN Torino
/// \author Stefano Politanò <stefano.politano@cern.ch>, CERN
/// \author Wu Chuntai <chuntai.wu@cern.ch>, CCNU, INFN Padova, and Padova University

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/HFC/DataModel/DerivedDataCorrelationTables.h"
#include "PWGHF/HFC/Utils/utilsFlow.h"
#include "PWGHF/Utils/utilsEvSelHf.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

#include <TF1.h>
#include <TString.h>

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_centrality;
using namespace o2::hf_evsel;
using namespace o2::hf_decay;
using namespace o2::hfc_utils_flow;

/// Code to produce self-contained derived datasets for v2 analyses of charm hadrons via 2PC
struct HfDerivedDataCreatorCorrelationsReduced {
  Produces<aod::HfcRedCorrColls> rowCollisions;      // Table with reduced collision info
  Produces<aod::HfcRedSEChBases> rowSECharmHadPairs; // Table with same-event pairs info
  Produces<aod::HfcRedSEHadBases> rowSEHadHadPairs;  // Table with same-event pairs info
  Produces<aod::HfcRedAssBases> rowAssocBases;       // Table with associated candidate base info
  Produces<aod::HfcRedAssTracks> rowAssocTrkSels;    // Table with associated track selection info
  Produces<aod::HfcRedTrigBases> rowTrigBases;       // Table with base trigger candidate info
  Produces<aod::HfcRedTrigCharms> rowTrigCharms;     // Table with charm trigger candidate selection info
  Produces<aod::HfcRedTrigTracks> rowTrigHads;       // Table with hadron trigger candidate selection info

  Configurable<int> centEstimator{"centEstimator", 2, "Centrality estimation (FT0A: 1, FT0C: 2, FT0M: 3, FV0A: 4)"};
  Configurable<int> selectionFlag{"selectionFlag", 15, "Selection Flag for hadron (ML score tables are required to run the task)"};
  Configurable<bool> forceCharmInCollision{"forceCharmInCollision", true, "Flag to force charm in collision"};
  Configurable<bool> acceptOnlyPrompt{"acceptOnlyPrompt", true, "Flag to accept only prompt charm in the MC process functions"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::vector<int>> classMl{"classMl", {0, 2}, "Indexes of BDT scores to be stored. Two indexes max."};
  Configurable<float> yCandMax{"yCandMax", 0.8, "max. cand. rapidity"};
  Configurable<float> ptCandMin{"ptCandMin", 0., "min. cand. pT"};
  Configurable<float> ptCandMax{"ptCandMax", 24., "max. cand. pT"};
  Configurable<int> tpcNClsCrossedRowsMin{"tpcNClsCrossedRowsMin", 70, "min. TPC crossed rows for associated tracks"};
  Configurable<float> etaTrkMax{"etaTrkMax", 1., "max. track eta"};
  Configurable<float> ptTrkMin{"ptTrkMin", 0.2, "min. track pT"};
  Configurable<float> ptTrkMax{"ptTrkMax", 3., "max. track pT"};
  Configurable<float> dcaXYTrkMax{"dcaXYTrkMax", 1., "max. track DCA XY"};
  Configurable<float> dcaZTrkMax{"dcaZTrkMax", 1., "max. track DCA Z"};
  Configurable<bool> usePtDiffDcaXYCut{"usePtDiffDcaXYCut", false, "Use pt-differential DCAxy cut for associated tracks"};
  Configurable<float> dcaXYTrkNSigmaMax{"dcaXYTrkNSigmaMax", 7, "Cut on number of sigma deviations from expected DCA in the transverse direction"};
  Configurable<std::string> dcaXYPtPrimTrkFunc{"dcaXYPtPrimTrkFunc", "(0.0026+0.005/(x^1.01))", "Functional form of pt-dependent DCAxy cut"};
  Configurable<float> deltaEtaAbsMin{"deltaEtaAbsMin", 0.5, "min. pair delta eta"};
  Configurable<float> deltaEtaAbsMax{"deltaEtaAbsMax", 2., "max. pair delta eta"};
  Configurable<float> downSampleTrksFactor{"downSampleTrksFactor", 1., "Fraction of associated tracks to keep"};
  Configurable<float> ptMaxForDownSample{"ptMaxForDownSample", 10., "Maximum pt for the application of the downsampling factor"};
  Configurable<float> centMaxForDownSample{"centMaxForDownSample", 101., "Maximum centrality for the application of the downsampling factor"};

  HfEventSelection hfEvSel; // event selection and monitoring
  HfcFlowHelper hfcFlowHelper; // flow utilities
  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb{};
  SliceCache cache;

  std::unique_ptr<TF1> funcDcaXYPtCutPrimTrk{nullptr};

  using DataCollsWithCentMult = soa::Join<aod::Collisions, aod::EvSels, aod::PVMults, aod::FT0Mults, aod::FV0Mults, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFV0As>;
  using McCollsWithCentMult = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::PVMults, aod::FT0Mults, aod::FV0Mults, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFV0As>;

  // Data
  using CandD0Data = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfMlD0>>;
  using CandDplusData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfMlDplusToPiKPi>>;
  using CandDsData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfMlDsToKKPi>>;
  using CandLcData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfMlLcToPKPi>>;
  using TracksData = soa::Filtered<soa::Join<aod::TracksWDca, aod::TrackSelection, aod::TracksExtra>>; //, aod::pidTPCFullMu, aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr>>;
  // MC
  using CandD0Mc = soa::Filtered<soa::Join<aod::HfCand2ProngMcRec, aod::HfCand2Prong, aod::HfSelD0, aod::HfMlD0>>;
  using CandDplusMc = soa::Filtered<soa::Join<aod::HfCand3ProngMcRec, aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfMlDplusToPiKPi>>;
  using CandDsMc = soa::Filtered<soa::Join<aod::HfCand3ProngMcRec, aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfMlDsToKKPi>>;
  using CandLcMc = soa::Filtered<soa::Join<aod::HfCand3ProngMcRec, aod::HfCand3Prong, aod::HfSelLc, aod::HfMlLcToPKPi>>;
  using TracksMc = soa::Filtered<soa::Join<aod::McTrackLabels, aod::TracksWDca, aod::TrackSelection, aod::TracksExtra>>; //, aod::pidTPCFullMu, aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr>>;

  Filter filterSelectD0Candidates = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlag || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlag;
  Filter filterSelectDplusCandidates = aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlag;
  Filter filterSelectDsCandidates = aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlag || aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlag;
  Filter filterSelectLcCandidates = aod::hf_sel_candidate_lc::isSelLcToPKPi >= selectionFlag || aod::hf_sel_candidate_lc::isSelLcToPiKP >= selectionFlag;
  Filter filterSelectTrkData = (nabs(aod::track::eta) < etaTrkMax) && (aod::track::pt > ptTrkMin) && (aod::track::pt < ptTrkMax) && (nabs(aod::track::dcaXY) < dcaXYTrkMax) && (nabs(aod::track::dcaZ) < dcaZTrkMax);

  Partition<CandD0Data> selD0ToPiK = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlag;
  Partition<CandD0Data> selD0ToKPi = aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlag;
  Partition<CandDsData> selDsToKKPi = aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlag;
  Partition<CandDsData> selDsToPiKK = aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlag;
  Partition<CandLcData> selLcToPKPi = aod::hf_sel_candidate_lc::isSelLcToPKPi >= selectionFlag;
  Partition<CandLcData> selLcToPiKP = aod::hf_sel_candidate_lc::isSelLcToPiKP >= selectionFlag;

  Partition<CandD0Mc> recoSelD0ToPiK = (nabs(aod::hf_cand_mc_flag::flagMcMatchRec) == static_cast<int8_t>(hf_cand_2prong::DecayChannelMain::D0ToPiK)) && aod::hf_sel_candidate_d0::isSelD0 >= selectionFlag && ((!acceptOnlyPrompt) || (acceptOnlyPrompt && aod::hf_cand_mc_flag::originMcRec == static_cast<int>(RecoDecay::OriginType::Prompt)));
  Partition<CandD0Mc> recoSelD0ToKPi = (nabs(aod::hf_cand_mc_flag::flagMcMatchRec) == static_cast<int8_t>(hf_cand_2prong::DecayChannelMain::D0ToPiK)) && aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlag && ((!acceptOnlyPrompt) || (acceptOnlyPrompt && aod::hf_cand_mc_flag::originMcRec == static_cast<int>(RecoDecay::OriginType::Prompt)));
  Partition<CandDplusMc> recoSelDPlusToPiKPi = (nabs(aod::hf_cand_mc_flag::flagMcMatchRec) == static_cast<int8_t>(hf_cand_3prong::DecayChannelMain::DplusToPiKPi)) && aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlag && ((!acceptOnlyPrompt) || (acceptOnlyPrompt && aod::hf_cand_mc_flag::originMcRec == static_cast<int>(RecoDecay::OriginType::Prompt)));
  Partition<CandDsMc> recoSelDsToKKPi = (nabs(aod::hf_cand_mc_flag::flagMcMatchRec) == static_cast<int8_t>(hf_cand_3prong::DecayChannelMain::DsToPiKK) && nabs(aod::hf_cand_mc_flag::flagMcDecayChanRec) == static_cast<int8_t>(hf_cand_3prong::DecayChannelResonant::DsToPhiPi)) && aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlag && ((!acceptOnlyPrompt) || (acceptOnlyPrompt && aod::hf_cand_mc_flag::originMcRec == static_cast<int>(RecoDecay::OriginType::Prompt)));
  Partition<CandDsMc> recoSelDsToPiKK = (nabs(aod::hf_cand_mc_flag::flagMcMatchRec) == static_cast<int8_t>(hf_cand_3prong::DecayChannelMain::DsToPiKK) && nabs(aod::hf_cand_mc_flag::flagMcDecayChanRec) == static_cast<int8_t>(hf_cand_3prong::DecayChannelResonant::DsToPhiPi)) && aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlag && ((!acceptOnlyPrompt) || (acceptOnlyPrompt && aod::hf_cand_mc_flag::originMcRec == static_cast<int>(RecoDecay::OriginType::Prompt)));
  Partition<CandLcMc> recoSelLcToPKPi = (nabs(aod::hf_cand_mc_flag::flagMcMatchRec) == static_cast<int8_t>(hf_cand_3prong::DecayChannelMain::LcToPKPi)) && aod::hf_sel_candidate_lc::isSelLcToPKPi >= selectionFlag && ((!acceptOnlyPrompt) || (acceptOnlyPrompt && aod::hf_cand_mc_flag::originMcRec == static_cast<int>(RecoDecay::OriginType::Prompt)));
  Partition<CandLcMc> recoSelLcToPiKP = (nabs(aod::hf_cand_mc_flag::flagMcMatchRec) == static_cast<int8_t>(hf_cand_3prong::DecayChannelMain::LcToPKPi)) && aod::hf_sel_candidate_lc::isSelLcToPiKP >= selectionFlag && ((!acceptOnlyPrompt) || (acceptOnlyPrompt && aod::hf_cand_mc_flag::originMcRec == static_cast<int>(RecoDecay::OriginType::Prompt)));

  ConfigurableAxis binsTpcSignal{"binsTpcSignal", {500, 0., 500.}, "TPC signal bins"};
  ConfigurableAxis binsPtTrig{"binsPtTrig", {24, 0, 24}, ""};
  ConfigurableAxis binsPtAssoc{"binsPtAssoc", {20, 0, 10}, ""};
  ConfigurableAxis binsMult{"binsMult", {100, 0., 10000.}, "Multiplicity bins"};
  ConfigurableAxis binsCent{"binsCent", {100, 0., 100.}, "Centrality bins"};
  ConfigurableAxis binsPosZ{"binsPosZ", {100, -10., 10.}, "Primary vertex z coordinate"};
  ConfigurableAxis binsEta{"binsEta", {50, -2., 2.}, "Eta bins"};
  ConfigurableAxis binsPhi{"binsPhi", {64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}, "Phi bins"};
  ConfigurableAxis binsDca{"binsDca", {200, -1., 1.}, ""};

  HistogramRegistry registry{"registry", {}};

  void init(InitContext&)
  {

    hfEvSel.addHistograms(registry); // collision monitoring
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    std::vector<unsigned int> classesMl;
    for (const auto scoreClass : classMl.value) {
      if (scoreClass < 0) {
        LOG(fatal) << "ML class indexes must be non-negative";
        return;
      }
      classesMl.push_back(static_cast<unsigned int>(scoreClass));
    }

    if (classesMl.size() > 2) {
      LOG(fatal) << "At most two ML classes can be selected";
      return;
    }

    hfcFlowHelper.init(ccdb, registry, classesMl, yCandMax.value);

    const AxisSpec axisTpcSignal = {binsTpcSignal, "TPC Signal"};
    const AxisSpec axisCent = {binsCent, "Centrality"};
    const AxisSpec axisMult = {binsMult, "Multiplicity"};
    const AxisSpec axisPosZ = {binsPosZ, "PosZ"};
    const AxisSpec axisEta = {binsEta, "#it{#eta}"};
    const AxisSpec axisPhi = {binsPhi, "#it{#varphi}"};
    const AxisSpec axisPtTrig = {binsPtTrig, "#it{p}_{T} Trig (GeV/#it{c})"};
    const AxisSpec axisPtAssoc = {binsPtAssoc, "#it{p}_{T} Assoc (GeV/#it{c})"};
    const AxisSpec axisDcaXY = {binsDca, "DCA XY (cm)"};
    const AxisSpec axisDcaZ = {binsDca, "DCA Z (cm)"};

    // Histograms for data analysis
    registry.add("hPhiVsPtTrig", "Trigger candidates phiVsPt", {HistType::kTH2F, {{axisPhi}, {axisPtTrig}}});
    registry.add("hEtaVsPtTrig", "Trigger candidates etaVsPt", {HistType::kTH2F, {{axisEta}, {axisPtTrig}}});
    registry.add("hPhiVsPtTrigAssoc", "Associated particles phiVsPt", {HistType::kTH3F, {{axisPhi}, {axisPtTrig}, {axisPtAssoc}}});
    registry.add("hEtaVsPtTrigAssoc", "Associated particles etaVsPt", {HistType::kTH3F, {{axisEta}, {axisPtTrig}, {axisPtAssoc}}});
    registry.add("hPhiVsPtAssoc", "Associated particles phiVsPt", {HistType::kTH2F, {{axisPhi}, {axisPtAssoc}}});
    registry.add("hEtaVsPtAssoc", "Associated particles etaVsPt", {HistType::kTH2F, {{axisEta}, {axisPtAssoc}}});
    registry.add("hDcaXYVsPtVsCentAssoc", "Associated particles DCAxyVsPtVsCent", {HistType::kTH3F, {{axisPtAssoc}, {axisDcaXY}, {axisCent}}});
    registry.add("hDcaZVsPtVsCentAssoc", "Associated particles DCAzVsPtVsCent", {HistType::kTH3F, {{axisPtAssoc}, {axisDcaZ}, {axisCent}}});
    registry.add("hDcaXYVsPtVsMultAssoc", "Associated particles DCAxyVsPtVsMult", {HistType::kTH3F, {{axisPtAssoc}, {axisDcaXY}, {axisMult}}});
    registry.add("hDcaZVsPtVsMultAssoc", "Associated particles DCAzVsPtVsMult", {HistType::kTH3F, {{axisPtAssoc}, {axisDcaZ}, {axisMult}}});

    // Setup pt-dependent DCAxy cut function
    if (usePtDiffDcaXYCut) {
      funcDcaXYPtCutPrimTrk = std::make_unique<TF1>("funcDcaXYPtCutPrimTrk", Form("[0]*%s", dcaXYPtPrimTrkFunc.value.data()), 0.001, 100.);
      funcDcaXYPtCutPrimTrk->SetParameter(0, dcaXYTrkNSigmaMax);
      LOGF(info, "DCAxy pt-dependence function: %s", Form("[0]*%s", dcaXYPtPrimTrkFunc.value.data()));
    }
  }; // end init

  /// Check event selections for collision and fill the collision table
  /// \param collision is the collision
  template <typename Coll>
  bool checkCollision(Coll const& collision, float& cent, float& mult)
  {
    o2::hf_evsel::HfCollisionRejectionMask collRejMask{};
    if (centEstimator == CentralityEstimator::FT0A) {
      collRejMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::FT0A, aod::BCsWithTimestamps>(collision, cent, ccdb, registry);
      mult = collision.multFT0A();
    } else if (centEstimator == CentralityEstimator::FT0C) {
      collRejMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::FT0C, aod::BCsWithTimestamps>(collision, cent, ccdb, registry);
      mult = collision.multFT0C();
    } else if (centEstimator == CentralityEstimator::FT0M) {
      collRejMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::FT0M, aod::BCsWithTimestamps>(collision, cent, ccdb, registry);
      mult = collision.multFT0M();
    } else if (centEstimator == CentralityEstimator::FV0A) {
      collRejMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::FV0A, aod::BCsWithTimestamps>(collision, cent, ccdb, registry);
      mult = collision.multFV0A();
    } else {
      LOG(fatal) << "Centrality estimator not recognized for collision selection";
      std::abort();
    }
    hfEvSel.fillHistograms(collision, collRejMask, cent);
    return collRejMask == 0;
  }

  /// Checks if the trigger cand-associated track pair can be accepted for SE correlation
  /// \param assTrk is the associated track
  /// \param cand is the trigger candidate
  template <TriggerChannels CandType, typename TCand, typename TAssocTrk>
  bool acceptSameEvtPair(TAssocTrk const& assTrk, TCand const& cand, double deltaEta)
  {
    if (std::abs(deltaEta) <= deltaEtaAbsMin || std::abs(deltaEta) > deltaEtaAbsMax) {
      return false;
    }

    if (!assTrk.isGlobalTrackWoDCA() || assTrk.tpcNClsCrossedRows() < tpcNClsCrossedRowsMin) {
      return false;
    }

    int const trackGlobalIndex = assTrk.globalIndex();
    if constexpr (CandType == TriggerChannels::Hadron) {
      if (!cand.isGlobalTrackWoDCA() || cand.tpcNClsCrossedRows() < tpcNClsCrossedRowsMin) {
        return false;
      }
      if (trackGlobalIndex <= cand.globalIndex()) {
        return false; // skip self-correlation and avoid pair duplication for hadron-hadron
      }
    } else {                                           // Remove Daughter-Cand pairs for charm-hadron correlations
      if constexpr ((requires { cand.prong2Id(); })) { // Check 3-prong
        if (trackGlobalIndex == cand.prong0Id() || trackGlobalIndex == cand.prong1Id() || trackGlobalIndex == cand.prong2Id()) {
          return false;
        }
      } else { // Check 2-prong
        if (trackGlobalIndex == cand.prong0Id() || trackGlobalIndex == cand.prong1Id()) {
          return false;
        }
      }
    }
    return true;
  }

  /// Fill histograms and tables for same-event correlations
  /// \param trigCands are the trigger candidates
  /// \param assTrks are the associated tracks
  /// \param cent is the collision centrality
  /// \param mult is the collision multiplicity
  template <bool IsMc, TriggerChannels CandType, typename TTrigCands, typename TAssocTrks>
  void fillSameEvent(TTrigCands const& trigCands,
                     TAssocTrks const& assTrks,
                     const float cent,
                     const float mult)
  {
    for (const auto& trigCand : trigCands) {
      if (!hfcFlowHelper.cutCandRapidity<CandType>(trigCand)) {
        continue;
      }
      double trigCandPt = trigCand.pt();
      registry.fill(HIST("hPhiVsPtTrig"), RecoDecay::constrainAngle(trigCand.phi(), -o2::constants::math::PIHalf), trigCandPt);
      registry.fill(HIST("hEtaVsPtTrig"), trigCand.eta(), trigCandPt);
      if constexpr (CandType == TriggerChannels::Hadron) {
        rowTrigHads(rowCollisions.lastIndex(), trigCandPt, trigCand.tpcNClsCrossedRows(), trigCand.itsClusterMap(), trigCand.itsNCls(), trigCand.dcaXY(), trigCand.dcaZ());
      } else {
        std::array<float, 2> outputMl = hfcFlowHelper.getCandMlScores<CandType>(trigCand);
        rowTrigCharms(rowCollisions.lastIndex(), trigCandPt, hfcFlowHelper.getCandMass<CandType>(trigCand), outputMl[0], outputMl[1]);
        rowTrigBases(trigCand.phi(), trigCand.eta());
      }

      for (const auto& assTrk : assTrks) {
        double assTrkPt = assTrk.pt();
        if (usePtDiffDcaXYCut) {
          float const dcaXYTrkCut = funcDcaXYPtCutPrimTrk->Eval(assTrkPt);
          if (std::fabs(assTrk.dcaXY()) > dcaXYTrkCut) {
            continue;
          }
        }

        double deltaEta = assTrk.eta() - trigCand.eta();
        if (!acceptSameEvtPair<CandType>(assTrk, trigCand, deltaEta)) {
          continue;
        }
        if (downSampleTrksFactor < 1.) {
          float const pseudoRndm = assTrkPt * 1000. - static_cast<int64_t>(assTrkPt * 1000);
          if (assTrkPt < ptMaxForDownSample && cent < centMaxForDownSample && pseudoRndm >= downSampleTrksFactor) {
            continue;
          }
        }
        registry.fill(HIST("hPhiVsPtTrigAssoc"), RecoDecay::constrainAngle(assTrk.phi(), -o2::constants::math::PIHalf), trigCandPt, assTrkPt);
        registry.fill(HIST("hEtaVsPtTrigAssoc"), assTrk.eta(), trigCandPt, assTrkPt);
        registry.fill(HIST("hPhiVsPtAssoc"), RecoDecay::constrainAngle(assTrk.phi(), -o2::constants::math::PIHalf), assTrkPt);
        registry.fill(HIST("hEtaVsPtAssoc"), assTrk.eta(), assTrkPt);
        registry.fill(HIST("hDcaXYVsPtVsCentAssoc"), assTrkPt, assTrk.dcaXY(), cent);
        registry.fill(HIST("hDcaXYVsPtVsMultAssoc"), assTrkPt, assTrk.dcaXY(), mult);
        registry.fill(HIST("hDcaZVsPtVsCentAssoc"), assTrkPt, assTrk.dcaZ(), cent);
        registry.fill(HIST("hDcaZVsPtVsMultAssoc"), assTrkPt, assTrk.dcaZ(), mult);

        double deltaPhi = RecoDecay::constrainAngle(assTrk.phi() - trigCand.phi(), -o2::constants::math::PIHalf);
        rowAssocBases(rowCollisions.lastIndex(), assTrk.phi(), assTrk.eta(), assTrkPt);
        rowAssocTrkSels(assTrk.tpcNClsCrossedRows(), assTrk.itsClusterMap(), assTrk.itsNCls(), assTrk.dcaXY(), assTrk.dcaZ());
        if constexpr (CandType == TriggerChannels::Hadron) {
          rowSEHadHadPairs(rowCollisions.lastIndex(), rowTrigHads.lastIndex(), assTrkPt, deltaEta, deltaPhi);
        } else {
          rowSECharmHadPairs(rowCollisions.lastIndex(), rowTrigCharms.lastIndex(), assTrkPt, deltaEta, deltaPhi);
        }
      }
    }
  }

  /// Fill charm hadron tables for mixed-event
  /// \param trigCands are the charm trigger candidates
  template <bool IsMc, TriggerChannels CandType, typename TTrigCands>
  void fillCharmMixedEvent(TTrigCands const& trigCands)
  {
    for (const auto& trigCand : trigCands) {
      if (!hfcFlowHelper.cutCandRapidity<CandType>(trigCand)) {
        continue;
      }
      registry.fill(HIST("hPhiVsPtTrig"), RecoDecay::constrainAngle(trigCand.phi(), -o2::constants::math::PIHalf), trigCand.pt());
      registry.fill(HIST("hEtaVsPtTrig"), trigCand.eta(), trigCand.pt());

      std::array<float, 2> outputMl = hfcFlowHelper.getCandMlScores<CandType>(trigCand);
      rowTrigBases(trigCand.phi(), trigCand.eta());
      rowTrigCharms(rowCollisions.lastIndex(), trigCand.pt(), hfcFlowHelper.getCandMass<CandType>(trigCand), outputMl[0], outputMl[1]);
    }
  }

  /// Fill track tables for mixed-event
  /// \param assTrks are the associated tracks
  /// \param cent is the collision centrality
  /// \param mult is the collision multiplicity
  template <bool IsMc, typename TAssocTrks>
  void fillTrkMixedEvent(TAssocTrks const& assTrks,
                         const float cent,
                         const float mult)
  {
    bool first = true;
    for (const auto& assTrk : assTrks) {
      if (!assTrk.isGlobalTrackWoDCA() || assTrk.tpcNClsCrossedRows() < tpcNClsCrossedRowsMin) {
        continue;
      }
      double assTrkPt = assTrk.pt();
      if (usePtDiffDcaXYCut) {
        float const dcaXYTrkCut = funcDcaXYPtCutPrimTrk->Eval(assTrkPt);
        if (std::fabs(assTrk.dcaXY()) > dcaXYTrkCut) {
          continue;
        }
      }
      if (!first && downSampleTrksFactor < 1.) { // skip downsampling for the first track to avoid empty tables
        float const pseudoRndm = assTrkPt * 1000. - static_cast<int64_t>(assTrkPt * 1000);
        if (assTrkPt < ptMaxForDownSample && cent < centMaxForDownSample && pseudoRndm >= downSampleTrksFactor) {
          continue;
        }
      }
      first = false;
      registry.fill(HIST("hPhiVsPtAssoc"), RecoDecay::constrainAngle(assTrk.phi(), -o2::constants::math::PIHalf), assTrkPt);
      registry.fill(HIST("hEtaVsPtAssoc"), assTrk.eta(), assTrkPt);
      registry.fill(HIST("hDcaXYVsPtVsCentAssoc"), assTrkPt, assTrk.dcaXY(), cent);
      registry.fill(HIST("hDcaZVsPtVsCentAssoc"), assTrkPt, assTrk.dcaZ(), cent);
      registry.fill(HIST("hDcaXYVsPtVsMultAssoc"), assTrkPt, assTrk.dcaXY(), mult);
      registry.fill(HIST("hDcaZVsPtVsMultAssoc"), assTrkPt, assTrk.dcaZ(), mult);
      rowAssocBases(rowCollisions.lastIndex(), assTrk.phi(), assTrk.eta(), assTrkPt);
      rowAssocTrkSels(assTrk.tpcNClsCrossedRows(), assTrk.itsClusterMap(), assTrk.itsNCls(), assTrk.dcaXY(), assTrk.dcaZ());
    }
  }

  // CHARM-HADRON CORRELATIONS PROCESS FUNCTIONS --- DATA
  // Dplus data
  void processDataDplusSE(DataCollsWithCentMult::iterator const& coll,
                          CandDplusData const& candsDplus,
                          TracksData const& tracks)
  {
    if (forceCharmInCollision && candsDplus.size() < 1) {
      return;
    }
    float cent{-1.}, mult{-1.};
    if (!checkCollision(coll, cent, mult)) {
      return;
    }
    rowCollisions(mult, coll.numContrib(), cent, coll.posZ());
    fillSameEvent<false, TriggerChannels::DplusToPiKPi>(candsDplus, tracks, cent, mult);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorCorrelationsReduced, processDataDplusSE, "Process Same Event for Dplus candidates", true);

  // Dplus data
  void processDataDplusME(DataCollsWithCentMult::iterator const& coll,
                          CandDplusData const& candsDplus,
                          TracksData const& tracks)
  {
    if (forceCharmInCollision && candsDplus.size() < 1) {
      return;
    }
    float cent{-1.}, mult{-1.};
    if (!checkCollision(coll, cent, mult)) {
      return;
    }
    rowCollisions(mult, coll.numContrib(), cent, coll.posZ());
    fillCharmMixedEvent<false, TriggerChannels::DplusToPiKPi>(candsDplus);
    fillTrkMixedEvent<false>(tracks, cent, mult);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorCorrelationsReduced, processDataDplusME, "Process Mixed Event for Dplus candidates", false);

  // Ds data
  void processDataDsSE(DataCollsWithCentMult::iterator const& coll,
                       TracksData const& tracks,
                       CandDsData const&)
  {
    auto candsDsToPiKK = selDsToPiKK->sliceByCached(aod::hf_cand::collisionId, coll.globalIndex(), cache);
    auto candsDsToKKPi = selDsToKKPi->sliceByCached(aod::hf_cand::collisionId, coll.globalIndex(), cache);
    if (forceCharmInCollision && candsDsToPiKK.size() < 1 && candsDsToKKPi.size() < 1) {
      return;
    }
    float cent{-1.}, mult{-1.};
    if (!checkCollision(coll, cent, mult)) {
      return;
    }
    rowCollisions(mult, coll.numContrib(), cent, coll.posZ());
    fillSameEvent<false, TriggerChannels::DsToPiKK>(candsDsToPiKK, tracks, cent, mult);
    fillSameEvent<false, TriggerChannels::DsToKKPi>(candsDsToKKPi, tracks, cent, mult);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorCorrelationsReduced, processDataDsSE, "Process Same Event for Ds candidates", false);

  // Ds data
  void processDataDsME(DataCollsWithCentMult::iterator const& coll,
                       TracksData const& tracks,
                       CandDsData const&)
  {
    auto candsDsToPiKK = selDsToPiKK->sliceByCached(aod::hf_cand::collisionId, coll.globalIndex(), cache);
    auto candsDsToKKPi = selDsToKKPi->sliceByCached(aod::hf_cand::collisionId, coll.globalIndex(), cache);
    if (forceCharmInCollision && candsDsToPiKK.size() < 1 && candsDsToKKPi.size() < 1) {
      return;
    }
    float cent{-1.}, mult{-1.};
    if (!checkCollision(coll, cent, mult)) {
      return;
    }
    rowCollisions(mult, coll.numContrib(), cent, coll.posZ());
    fillCharmMixedEvent<false, TriggerChannels::DsToPiKK>(candsDsToPiKK);
    fillCharmMixedEvent<false, TriggerChannels::DsToKKPi>(candsDsToKKPi);
    fillTrkMixedEvent<false>(tracks, cent, mult);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorCorrelationsReduced, processDataDsME, "Process Mixed Event for Ds candidates", false);

  // D0 data
  void processDataD0SE(DataCollsWithCentMult::iterator const& coll,
                       TracksData const& tracks,
                       CandD0Data const&)
  {
    auto candsD0ToPiK = selD0ToPiK->sliceByCached(aod::hf_cand::collisionId, coll.globalIndex(), cache);
    auto candsD0ToKPi = selD0ToKPi->sliceByCached(aod::hf_cand::collisionId, coll.globalIndex(), cache);
    if (forceCharmInCollision && candsD0ToPiK.size() < 1 && candsD0ToKPi.size() < 1) {
      return;
    }
    float cent{-1.}, mult{-1.};
    if (!checkCollision(coll, cent, mult)) {
      return;
    }
    rowCollisions(mult, coll.numContrib(), cent, coll.posZ());
    fillSameEvent<false, TriggerChannels::D0ToPiK>(candsD0ToPiK, tracks, cent, mult);
    fillSameEvent<false, TriggerChannels::D0ToKPi>(candsD0ToKPi, tracks, cent, mult);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorCorrelationsReduced, processDataD0SE, "Process Same Event for D0 candidates", false);

  // D0 data
  void processDataD0ME(DataCollsWithCentMult::iterator const& coll,
                       TracksData const& tracks,
                       CandD0Data const&)
  {
    auto candsD0ToPiK = selD0ToPiK->sliceByCached(aod::hf_cand::collisionId, coll.globalIndex(), cache);
    auto candsD0ToKPi = selD0ToKPi->sliceByCached(aod::hf_cand::collisionId, coll.globalIndex(), cache);
    if (forceCharmInCollision && candsD0ToPiK.size() < 1 && candsD0ToKPi.size() < 1) {
      return;
    }
    float cent{-1.}, mult{-1.};
    if (!checkCollision(coll, cent, mult)) {
      return;
    }
    rowCollisions(mult, coll.numContrib(), cent, coll.posZ());
    fillCharmMixedEvent<false, TriggerChannels::D0ToPiK>(candsD0ToPiK);
    fillCharmMixedEvent<false, TriggerChannels::D0ToKPi>(candsD0ToKPi);
    fillTrkMixedEvent<false>(tracks, cent, mult);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorCorrelationsReduced, processDataD0ME, "Process Mixed Event for D0 candidates", false);

  // Lc data
  void processDataLcSE(DataCollsWithCentMult::iterator const& coll,
                       TracksData const& tracks,
                       CandLcData const&)
  {
    auto candsLcToPKPi = selLcToPKPi->sliceByCached(aod::hf_cand::collisionId, coll.globalIndex(), cache);
    auto candsLcToPiKP = selLcToPiKP->sliceByCached(aod::hf_cand::collisionId, coll.globalIndex(), cache);
    if (forceCharmInCollision && candsLcToPKPi.size() < 1 && candsLcToPiKP.size() < 1) {
      return;
    }
    float cent{-1.}, mult{-1.};
    if (!checkCollision(coll, cent, mult)) {
      return;
    }
    rowCollisions(mult, coll.numContrib(), cent, coll.posZ());
    fillSameEvent<false, TriggerChannels::LcToPKPi>(candsLcToPKPi, tracks, cent, mult);
    fillSameEvent<false, TriggerChannels::LcToPiKP>(candsLcToPiKP, tracks, cent, mult);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorCorrelationsReduced, processDataLcSE, "Process Same Event for Lc candidates", false);

  // Lc data
  void processDataLcME(DataCollsWithCentMult::iterator const& coll,
                       TracksData const& tracks,
                       CandLcData const&)
  {
    auto candsLcToPKPi = selLcToPKPi->sliceByCached(aod::hf_cand::collisionId, coll.globalIndex(), cache);
    auto candsLcToPiKP = selLcToPiKP->sliceByCached(aod::hf_cand::collisionId, coll.globalIndex(), cache);
    if (forceCharmInCollision && candsLcToPKPi.size() < 1 && candsLcToPiKP.size() < 1) {
      return;
    }
    float cent{-1.}, mult{-1.};
    if (!checkCollision(coll, cent, mult)) {
      return;
    }
    rowCollisions(mult, coll.numContrib(), cent, coll.posZ());
    fillCharmMixedEvent<false, TriggerChannels::LcToPKPi>(candsLcToPKPi);
    fillCharmMixedEvent<false, TriggerChannels::LcToPiKP>(candsLcToPiKP);
    fillTrkMixedEvent<false>(tracks, cent, mult);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorCorrelationsReduced, processDataLcME, "Process Mixed Event for Lc candidates", false);

  // CHARM-HADRON CORRELATIONS PROCESS FUNCTIONS --- MC
  // Dplus MC
  void processMcDplusSE(McCollsWithCentMult::iterator const& coll,
                        TracksMc const& tracks,
                        CandDplusMc const&)
  {
    auto candsRecoDplus = recoSelDPlusToPiKPi->sliceByCached(aod::hf_cand::collisionId, coll.globalIndex(), cache);
    if (forceCharmInCollision && candsRecoDplus.size() < 1) {
      return;
    }
    float cent{-1.}, mult{-1.};
    if (!checkCollision(coll, cent, mult)) {
      return;
    }
    rowCollisions(mult, coll.numContrib(), cent, coll.posZ());
    fillSameEvent<true, TriggerChannels::DplusToPiKPi>(candsRecoDplus, tracks, cent, mult);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorCorrelationsReduced, processMcDplusSE, "Process Same Event for Dplus candidates", false);

  // Dplus MC
  void processMcDplusME(McCollsWithCentMult::iterator const& coll,
                        TracksMc const& tracks,
                        CandDplusMc const&)
  {
    auto candsRecoDplus = recoSelDPlusToPiKPi->sliceByCached(aod::hf_cand::collisionId, coll.globalIndex(), cache);
    if (forceCharmInCollision && candsRecoDplus.size() < 1) {
      return;
    }
    float cent{-1.}, mult{-1.};
    if (!checkCollision(coll, cent, mult)) {
      return;
    }
    rowCollisions(mult, coll.numContrib(), cent, coll.posZ());
    fillCharmMixedEvent<true, TriggerChannels::DplusToPiKPi>(candsRecoDplus);
    fillTrkMixedEvent<true>(tracks, cent, mult);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorCorrelationsReduced, processMcDplusME, "Process Mixed Event for Dplus candidates", false);

  // Ds MC
  void processMcDsSE(McCollsWithCentMult::iterator const& coll,
                     TracksMc const& tracks,
                     CandDsMc const&)
  {
    auto candsRecoDsToPiKK = recoSelDsToPiKK->sliceByCached(aod::hf_cand::collisionId, coll.globalIndex(), cache);
    auto candsRecoDsToKKPi = recoSelDsToKKPi->sliceByCached(aod::hf_cand::collisionId, coll.globalIndex(), cache);
    if (forceCharmInCollision && candsRecoDsToPiKK.size() < 1 && candsRecoDsToKKPi.size() < 1) {
      return;
    }
    float cent{-1.}, mult{-1.};
    if (!checkCollision(coll, cent, mult)) {
      return;
    }
    rowCollisions(mult, coll.numContrib(), cent, coll.posZ());
    fillSameEvent<true, TriggerChannels::DsToPiKK>(candsRecoDsToPiKK, tracks, cent, mult);
    fillSameEvent<true, TriggerChannels::DsToKKPi>(candsRecoDsToKKPi, tracks, cent, mult);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorCorrelationsReduced, processMcDsSE, "Process Same Event for Ds candidates", false);

  // Ds MC
  void processMcDsME(McCollsWithCentMult::iterator const& coll,
                     TracksMc const& tracks,
                     CandDsMc const&)
  {
    auto candsRecoDsToPiKK = recoSelDsToPiKK->sliceByCached(aod::hf_cand::collisionId, coll.globalIndex(), cache);
    auto candsRecoDsToKKPi = recoSelDsToKKPi->sliceByCached(aod::hf_cand::collisionId, coll.globalIndex(), cache);
    if (forceCharmInCollision && candsRecoDsToPiKK.size() < 1 && candsRecoDsToKKPi.size() < 1) {
      return;
    }
    float cent{-1.}, mult{-1.};
    if (!checkCollision(coll, cent, mult)) {
      return;
    }
    rowCollisions(mult, coll.numContrib(), cent, coll.posZ());
    fillCharmMixedEvent<true, TriggerChannels::DsToPiKK>(candsRecoDsToPiKK);
    fillCharmMixedEvent<true, TriggerChannels::DsToKKPi>(candsRecoDsToKKPi);
    fillTrkMixedEvent<true>(tracks, cent, mult);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorCorrelationsReduced, processMcDsME, "Process Mixed Event for Ds candidates", false);

  // D0 MC
  void processMcD0SE(McCollsWithCentMult::iterator const& coll,
                     TracksMc const& tracks,
                     CandD0Mc const&)
  {
    auto candsRecoD0ToPiK = recoSelD0ToPiK->sliceByCached(aod::hf_cand::collisionId, coll.globalIndex(), cache);
    auto candsRecoD0ToKPi = recoSelD0ToKPi->sliceByCached(aod::hf_cand::collisionId, coll.globalIndex(), cache);
    if (forceCharmInCollision && candsRecoD0ToPiK.size() < 1 && candsRecoD0ToKPi.size() < 1) {
      return;
    }
    float cent{-1.}, mult{-1.};
    if (!checkCollision(coll, cent, mult)) {
      return;
    }
    rowCollisions(mult, coll.numContrib(), cent, coll.posZ());
    fillSameEvent<true, TriggerChannels::D0ToPiK>(candsRecoD0ToPiK, tracks, cent, mult);
    fillSameEvent<true, TriggerChannels::D0ToKPi>(candsRecoD0ToKPi, tracks, cent, mult);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorCorrelationsReduced, processMcD0SE, "Process Same Event for D0 candidates", false);

  // D0 MC
  void processMcD0ME(McCollsWithCentMult::iterator const& coll,
                     TracksMc const& tracks,
                     CandD0Mc const&)
  {
    auto candsRecoD0ToPiK = recoSelD0ToPiK->sliceByCached(aod::hf_cand::collisionId, coll.globalIndex(), cache);
    auto candsRecoD0ToKPi = recoSelD0ToKPi->sliceByCached(aod::hf_cand::collisionId, coll.globalIndex(), cache);
    if (forceCharmInCollision && candsRecoD0ToPiK.size() < 1 && candsRecoD0ToKPi.size() < 1) {
      return;
    }
    float cent{-1.}, mult{-1.};
    if (!checkCollision(coll, cent, mult)) {
      return;
    }
    rowCollisions(mult, coll.numContrib(), cent, coll.posZ());
    fillCharmMixedEvent<true, TriggerChannels::D0ToPiK>(candsRecoD0ToPiK);
    fillCharmMixedEvent<true, TriggerChannels::D0ToKPi>(candsRecoD0ToKPi);
    fillTrkMixedEvent<true>(tracks, cent, mult);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorCorrelationsReduced, processMcD0ME, "Process Mixed Event for D0 candidates", false);

  // Lc MC
  void processMcLcSE(McCollsWithCentMult::iterator const& coll,
                     TracksMc const& tracks,
                     CandLcMc const&)
  {
    auto candsRecoLcToPKPi = recoSelLcToPKPi->sliceByCached(aod::hf_cand::collisionId, coll.globalIndex(), cache);
    auto candsRecoLcToPiKP = recoSelLcToPiKP->sliceByCached(aod::hf_cand::collisionId, coll.globalIndex(), cache);
    if (forceCharmInCollision && candsRecoLcToPKPi.size() < 1 && candsRecoLcToPiKP.size() < 1) {
      return;
    }
    float cent{-1.}, mult{-1.};
    if (!checkCollision(coll, cent, mult)) {
      return;
    }
    rowCollisions(mult, coll.numContrib(), cent, coll.posZ());
    fillSameEvent<true, TriggerChannels::LcToPKPi>(candsRecoLcToPKPi, tracks, cent, mult);
    fillSameEvent<true, TriggerChannels::LcToPiKP>(candsRecoLcToPiKP, tracks, cent, mult);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorCorrelationsReduced, processMcLcSE, "Process Same Event for Lc candidates", false);

  // Lc MC
  void processMcLcME(McCollsWithCentMult::iterator const& coll,
                     TracksMc const& tracks,
                     CandLcMc const&)
  {
    auto candsRecoLcToPKPi = recoSelLcToPKPi->sliceByCached(aod::hf_cand::collisionId, coll.globalIndex(), cache);
    auto candsRecoLcToPiKP = recoSelLcToPiKP->sliceByCached(aod::hf_cand::collisionId, coll.globalIndex(), cache);
    if (forceCharmInCollision && candsRecoLcToPKPi.size() < 1 && candsRecoLcToPiKP.size() < 1) {
      return;
    }
    float cent{-1.}, mult{-1.};
    if (!checkCollision(coll, cent, mult)) {
      return;
    }
    rowCollisions(mult, coll.numContrib(), cent, coll.posZ());
    fillCharmMixedEvent<true, TriggerChannels::LcToPKPi>(candsRecoLcToPKPi);
    fillCharmMixedEvent<true, TriggerChannels::LcToPiKP>(candsRecoLcToPiKP);
    fillTrkMixedEvent<true>(tracks, cent, mult);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorCorrelationsReduced, processMcLcME, "Process Mixed Event for Lc candidates", false);

  // HADRON-HADRON CORRELATIONS PROCESS FUNCTIONS --- DATA
  // Hadron Hadron Same Event
  void processDataHadronHadronSE(DataCollsWithCentMult::iterator const& coll,
                                 TracksData const& tracks)
  {
    float cent{-1.}, mult{-1.};
    if (!checkCollision(coll, cent, mult)) {
      return;
    }
    rowCollisions(mult, coll.numContrib(), cent, coll.posZ());
    fillSameEvent<false, TriggerChannels::Hadron>(tracks, tracks, cent, mult);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorCorrelationsReduced, processDataHadronHadronSE, "Process Same Event for hadron candidates", false);

  // Hadron Hadron Mixed Event
  void processDataHadronHadronME(DataCollsWithCentMult::iterator const& coll,
                                 TracksData const& tracks)
  {
    float cent{-1.}, mult{-1.};
    if (!checkCollision(coll, cent, mult)) {
      return;
    }
    rowCollisions(mult, coll.numContrib(), cent, coll.posZ());
    fillTrkMixedEvent<false>(tracks, cent, mult);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorCorrelationsReduced, processDataHadronHadronME, "Process Mixed Event for hadron candidates", false);

  // HADRON-HADRON CORRELATIONS PROCESS FUNCTIONS --- MC
  // Hadron Hadron Same Event
  void processMcHadronHadronSE(McCollsWithCentMult::iterator const& coll,
                               TracksMc const& tracks)
  {
    float cent{-1.}, mult{-1.};
    if (!checkCollision(coll, cent, mult)) {
      return;
    }
    rowCollisions(mult, coll.numContrib(), cent, coll.posZ());
    fillSameEvent<true, TriggerChannels::Hadron>(tracks, tracks, cent, mult);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorCorrelationsReduced, processMcHadronHadronSE, "Process Same Event for hadron candidates", false);

  // Hadron Hadron Mixed Event
  void processMcHadronHadronME(McCollsWithCentMult::iterator const& coll,
                               TracksMc const& tracks)
  {
    float cent{-1.}, mult{-1.};
    if (!checkCollision(coll, cent, mult)) {
      return;
    }
    rowCollisions(mult, coll.numContrib(), cent, coll.posZ());
    fillTrkMixedEvent<true>(tracks, cent, mult);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorCorrelationsReduced, processMcHadronHadronME, "Process Mixed Event for hadron candidates", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfDerivedDataCreatorCorrelationsReduced>(cfgc)};
}
