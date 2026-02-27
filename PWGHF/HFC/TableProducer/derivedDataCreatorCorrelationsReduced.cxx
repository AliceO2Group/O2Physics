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
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/HFC/DataModel/DerivedDataCorrelationTables.h"
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

#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_centrality;
using namespace o2::hf_evsel;

enum CandidateType {
  DplusToPiKPi = 0,
  DsToKKPi,
  DsToPiKK,
  D0ToPiK,
  D0ToKPi,
  LcToPKPi,
  Hadron
};

/// Code to select collisions with at least one Ds meson
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
  Configurable<std::vector<double>> binsPtTrig{"binsPtTrig", std::vector<double>{0., 1., 2., 3., 5., 8., 12., 24., 36.}, "pT bin limits for trigger candidates"};
  Configurable<std::vector<double>> binsPtAssoc{"binsPtAssoc", std::vector<double>{0.2, 1., 2., 50.}, "pT bin limits for associated particles"};

  HfEventSelection hfEvSel; // event selection and monitoring
  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb;
  SliceCache cache;

  double massCharm{0.};
  TF1* funcDcaXYPtCutPrimTrk = nullptr;

  using CollsWithCentMult = soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFV0As>;
  using CandDsData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfMlDsToKKPi>>;
  using CandDplusData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfMlDplusToPiKPi>>;
  using CandD0Data = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfMlD0>>;
  using CandLcData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfMlLcToPKPi>>;
  using TracksData = soa::Filtered<soa::Join<aod::TracksWDca, aod::TrackSelection, aod::TracksExtra>>;

  Filter filterSelectDsCandidates = aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlag || aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlag;
  Filter filterSelectDplusCandidates = aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlag;
  Filter filterSelectD0Candidates = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlag || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlag;
  Filter filterSelectLcCandidates = aod::hf_sel_candidate_lc::isSelLcToPKPi >= selectionFlag;
  Filter filterSelectTrkData = (nabs(aod::track::eta) < etaTrkMax) && (aod::track::pt > ptTrkMin) && (aod::track::pt < ptTrkMax) && (nabs(aod::track::dcaXY) < dcaXYTrkMax) && (nabs(aod::track::dcaZ) < dcaZTrkMax);

  Preslice<CandDsData> candsDsPerColl = aod::hf_cand::collisionId;
  Preslice<CandDplusData> candsDplusPerColl = aod::hf_cand::collisionId;
  Preslice<CandD0Data> candsD0PerColl = aod::hf_cand::collisionId;
  Preslice<CandLcData> candsLcPerColl = aod::hf_cand::collisionId;
  Preslice<TracksData> trackIndicesPerColl = aod::track::collisionId;

  Partition<CandDsData> selectedDsToKKPi = aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlag;
  Partition<CandDsData> selectedDsToPiKK = aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlag;
  Partition<CandD0Data> selectedD0ToPiK = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlag;
  Partition<CandD0Data> selectedD0ToKPi = aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlag;
  Partition<CandLcData> selectedLcToPKPi = aod::hf_sel_candidate_lc::isSelLcToPKPi >= selectionFlag;

  ConfigurableAxis binsInvMass{"binsInvMass", {300, 1.6, 2.2}, ""};
  ConfigurableAxis binsMultFT0M{"binsMultFT0M", {100, 0., 10000.}, "Multiplicity as FT0M signal amplitude"};
  ConfigurableAxis binsCent{"binsCent", {100, 0., 100.}, "Centrality bins"};
  ConfigurableAxis binsPosZ{"binsPosZ", {100, -10., 10.}, "Primary vertex z coordinate"};
  ConfigurableAxis binsEta{"binsEta", {50, -2., 2.}, "Eta bins"};
  ConfigurableAxis binsPhi{"binsPhi", {64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}, "Phi bins"};
  ConfigurableAxis binsDeltaEta{"binsDeltaEta", {100, -2., 2.}, "Delta Eta bins"};
  ConfigurableAxis binsDeltaPhi{"binsDeltaPhi", {64, -3., 3.}, "Delta Phi bins"};
  ConfigurableAxis binsMlOne{"binsMlOne", {100, 0., 1.}, ""};
  ConfigurableAxis binsMlTwo{"binsMlTwo", {100, 0., 1.}, ""};
  ConfigurableAxis binsDca{"binsDca", {200, -1., 1.}, ""};

  HistogramRegistry registry{"registry", {}};

  void init(InitContext&)
  {
    if (doprocessDplusSameEvent || doprocessDplusMixedEvent) {
      massCharm = o2::constants::physics::MassDPlus;
    } else if (doprocessDsSameEvent || doprocessDsMixedEvent) {
      massCharm = o2::constants::physics::MassDS;
    } else if (doprocessD0SameEvent || doprocessD0MixedEvent) {
      massCharm = o2::constants::physics::MassD0;
    } else if (doprocessLcSameEvent || doprocessLcMixedEvent) {
      massCharm = o2::constants::physics::MassLambdaCPlus;
    } else if (doprocessHadronHadronSameEvent || doprocessHadronHadronMixedEvent) {
      LOG(info) << "Charm mass not set, processing Hadron-Hadron case";
    } else {
      LOG(fatal) << "No decay channel selected to process";
    }

    hfEvSel.addHistograms(registry); // collision monitoring
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    const AxisSpec axisCent = {binsCent, "Centrality"};
    const AxisSpec axisMultFT0M = {binsMultFT0M, "MultiplicityFT0M"};
    const AxisSpec axisPosZ = {binsPosZ, "PosZ"};
    const AxisSpec axisEta = {binsEta, "#it{#eta}"};
    const AxisSpec axisPhi = {binsPhi, "#it{#varphi}"};
    const AxisSpec axisPtTrig = {(std::vector<double>)binsPtTrig, "#it{p}_{T} Trig (GeV/#it{c})"};
    const AxisSpec axisPtAssoc = {(std::vector<double>)binsPtAssoc, "#it{p}_{T} Assoc (GeV/#it{c})"};
    const AxisSpec axisDcaXY = {binsDca, "DCA XY (cm)"};
    const AxisSpec axisDcaZ = {binsDca, "DCA Z (cm)"};

    // Histograms for data analysis
    registry.add("hPhiVsPtTrig", "Trigger candidates phiVsPt", {HistType::kTH2F, {{axisPhi}, {axisPtTrig}}});
    registry.add("hEtaVsPtTrig", "Trigger candidates etaVsPt", {HistType::kTH2F, {{axisEta}, {axisPtTrig}}});
    registry.add("hPhiVsPtTrigAssoc", "Associated particles phiVsPt", {HistType::kTH3F, {{axisPhi}, {axisPtTrig}, {axisPtAssoc}}});
    registry.add("hEtaVsPtTrigAssoc", "Associated particles etaVsPt", {HistType::kTH3F, {{axisEta}, {axisPtTrig}, {axisPtAssoc}}});
    registry.add("hPhiVsPtAssoc", "Associated particles phiVsPt", {HistType::kTH2F, {{axisPhi}, {axisPtAssoc}}});
    registry.add("hEtaVsPtAssoc", "Associated particles etaVsPt", {HistType::kTH2F, {{axisEta}, {axisPtAssoc}}});
    registry.add("hDcaXYVsPtAssoc", "Associated particles DCAxyVsPt", {HistType::kTH2F, {{axisDcaXY}, {axisPtAssoc}}});
    registry.add("hDcaZVsPtAssoc", "Associated particles DCAzVsPt", {HistType::kTH2F, {{axisDcaZ}, {axisPtAssoc}}});

    // Setup pt-dependent DCAxy cut function
    if (usePtDiffDcaXYCut) {
      funcDcaXYPtCutPrimTrk = new TF1("funcDcaXYPtCutPrimTrk", Form("[0]*%s", dcaXYPtPrimTrkFunc.value.data()), 0.001, 100);
      funcDcaXYPtCutPrimTrk->SetParameter(0, dcaXYTrkNSigmaMax);
      LOGF(info, "DCAxy pt-dependence function: %s", Form("[0]*%s", dcaXYPtPrimTrkFunc.value.data()));
    }
  }; // end init

  /// Get charm hadron candidate mass
  /// \param candidate is the charm hadron candidate
  template <CandidateType CandType, typename TCand>
  double getCandMass(const TCand& candidate)
  {
    if constexpr (CandType == CandidateType::DsToKKPi) {
      return HfHelper::invMassDsToKKPi(candidate);
    }
    if constexpr (CandType == CandidateType::DsToPiKK) {
      return HfHelper::invMassDsToPiKK(candidate);
    }
    if constexpr (CandType == CandidateType::DplusToPiKPi) {
      return HfHelper::invMassDplusToPiKPi(candidate);
    }
    if constexpr (CandType == CandidateType::D0ToPiK) {
      return HfHelper::invMassD0ToPiK(candidate);
    }
    if constexpr (CandType == CandidateType::D0ToKPi) {
      return HfHelper::invMassD0barToKPi(candidate);
    }
    if constexpr (CandType == CandidateType::LcToPKPi) {
      return HfHelper::invMassLcToPKPi(candidate);
    }
    return -1.;
  }

  /// Get charm hadron bdt scores
  /// \param candidate is the charm hadron candidate
  template <CandidateType CandType, typename TCand>
  std::array<float, 2> getCandMlScores(const TCand& candidate)
  {
    std::array<float, 2> outputMl{-999.f, -999.f};
    if constexpr (CandType == CandidateType::DsToKKPi) {
      for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
        outputMl[iclass] = candidate.mlProbDsToKKPi()[classMl->at(iclass)];
      }
    }
    if constexpr (CandType == CandidateType::DsToPiKK) {
      for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
        outputMl[iclass] = candidate.mlProbDsToPiKK()[classMl->at(iclass)];
      }
    }
    if constexpr (CandType == CandidateType::DplusToPiKPi) {
      for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
        outputMl[iclass] = candidate.mlProbDplusToPiKPi()[classMl->at(iclass)];
      }
    }
    if constexpr (CandType == CandidateType::D0ToPiK) {
      for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
        outputMl[iclass] = candidate.mlProbD0()[classMl->at(iclass)];
      }
    }
    if constexpr (CandType == CandidateType::D0ToKPi) {
      for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
        outputMl[iclass] = candidate.mlProbD0bar()[classMl->at(iclass)];
      }
    }
    if constexpr (CandType == CandidateType::LcToPKPi) {
      for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
        outputMl[iclass] = candidate.mlProbLcToPKPi()[classMl->at(iclass)];
      }
    }
    return outputMl;
  }

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
  template <CandidateType CandType, typename TCand, typename TAssocTrk>
  bool acceptSameEvtPair(TAssocTrk const& assTrk, TCand const& cand, double deltaEta)
  {
    if (std::abs(deltaEta) <= deltaEtaAbsMin || std::abs(deltaEta) > deltaEtaAbsMax) {
      return false;
    }

    if (!assTrk.isGlobalTrackWoDCA() || assTrk.tpcNClsCrossedRows() < tpcNClsCrossedRowsMin) {
      return false;
    }

    int const trackGlobalIndex = assTrk.globalIndex();
    if constexpr (CandType == CandidateType::Hadron) {
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
  /// \param collCentrality is the collision centrality
  template <CandidateType CandType, typename TTrigCands, typename TAssocTrks>
  void fillSameEvent(TTrigCands const& trigCands,
                     TAssocTrks const& assTrks,
                     const float collCentrality)
  {
    for (const auto& trigCand : trigCands) {
      double trigCandPt = trigCand.pt();
      registry.fill(HIST("hPhiVsPtTrig"), RecoDecay::constrainAngle(trigCand.phi(), -o2::constants::math::PIHalf), trigCandPt);
      registry.fill(HIST("hEtaVsPtTrig"), trigCand.eta(), trigCandPt);
      if constexpr (CandType == CandidateType::Hadron) {
        rowTrigHads(rowCollisions.lastIndex(), trigCandPt, trigCand.tpcNClsCrossedRows(), trigCand.itsClusterMap(), trigCand.itsNCls(), trigCand.dcaXY(), trigCand.dcaZ());
      } else {
        std::array<float, 2> outputMl = getCandMlScores<CandType>(trigCand);
        rowTrigCharms(rowCollisions.lastIndex(), trigCandPt, getCandMass<CandType>(trigCand), outputMl[0], outputMl[1]);
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
          if (assTrkPt < ptMaxForDownSample && collCentrality < centMaxForDownSample && pseudoRndm >= downSampleTrksFactor) {
            continue;
          }
        }
        registry.fill(HIST("hPhiVsPtTrigAssoc"), RecoDecay::constrainAngle(assTrk.phi(), -o2::constants::math::PIHalf), trigCandPt, assTrkPt);
        registry.fill(HIST("hEtaVsPtTrigAssoc"), assTrk.eta(), trigCandPt, assTrkPt);
        registry.fill(HIST("hPhiVsPtAssoc"), RecoDecay::constrainAngle(assTrk.phi(), -o2::constants::math::PIHalf), assTrkPt);
        registry.fill(HIST("hEtaVsPtAssoc"), assTrk.eta(), assTrkPt);
        registry.fill(HIST("hDcaXYVsPtAssoc"), assTrk.dcaXY(), assTrkPt);
        registry.fill(HIST("hDcaZVsPtAssoc"), assTrk.dcaZ(), assTrkPt);

        double deltaPhi = RecoDecay::constrainAngle(assTrk.phi() - trigCand.phi(), -o2::constants::math::PIHalf);
        rowAssocTrkSels(assTrk.tpcNClsCrossedRows(), assTrk.itsClusterMap(), assTrk.itsNCls(), assTrk.dcaXY(), assTrk.dcaZ());
        if constexpr (CandType == CandidateType::Hadron) {
          rowSEHadHadPairs(rowCollisions.lastIndex(), rowTrigHads.lastIndex(), assTrkPt, deltaEta, deltaPhi);
        } else {
          rowSECharmHadPairs(rowCollisions.lastIndex(), rowTrigCharms.lastIndex(), assTrkPt, deltaEta, deltaPhi);
        }
      }
    }
  }

  /// Fill charm hadron tables for mixed-event
  /// \param trigCands are the charm trigger candidates
  template <CandidateType CandType, typename TTrigCands>
  void fillCharmMixedEvent(TTrigCands const& trigCands)
  {
    for (const auto& trigCand : trigCands) {
      registry.fill(HIST("hPhiVsPtTrig"), RecoDecay::constrainAngle(trigCand.phi(), -o2::constants::math::PIHalf), trigCand.pt());
      registry.fill(HIST("hEtaVsPtTrig"), trigCand.eta(), trigCand.pt());

      std::array<float, 2> outputMl = getCandMlScores<CandType>(trigCand);
      rowTrigBases(trigCand.phi(), trigCand.eta());
      rowTrigCharms(rowCollisions.lastIndex(), trigCand.pt(), getCandMass<CandType>(trigCand), outputMl[0], outputMl[1]);
    }
  }

  /// Fill track tables for mixed-event
  /// \param assTrks are the associated tracks
  /// \param collCentrality is the collision centrality
  template <typename TAssocTrks>
  void fillTrkMixedEvent(TAssocTrks const& assTrks,
                         const float collCentrality)
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
        if (assTrkPt < ptMaxForDownSample && collCentrality < centMaxForDownSample && pseudoRndm >= downSampleTrksFactor) {
          continue;
        }
      }
      first = false;
      registry.fill(HIST("hPhiVsPtAssoc"), RecoDecay::constrainAngle(assTrk.phi(), -o2::constants::math::PIHalf), assTrkPt);
      registry.fill(HIST("hEtaVsPtAssoc"), assTrk.eta(), assTrkPt);
      registry.fill(HIST("hDcaXYVsPtAssoc"), assTrk.dcaXY(), assTrkPt);
      registry.fill(HIST("hDcaZVsPtAssoc"), assTrk.dcaZ(), assTrkPt);
      rowAssocBases(rowCollisions.lastIndex(), assTrk.phi(), assTrk.eta(), assTrkPt);
      rowAssocTrkSels(assTrk.tpcNClsCrossedRows(), assTrk.itsClusterMap(), assTrk.itsNCls(), assTrk.dcaXY(), assTrk.dcaZ());
    }
  }

  // Dplus with ML selections
  void processDplusSameEvent(CollsWithCentMult::iterator const& coll,
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
    fillSameEvent<CandidateType::DplusToPiKPi>(candsDplus, tracks, cent);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorCorrelationsReduced, processDplusSameEvent, "Process Same Event for Dplus candidates", true);

  // Dplus with ML selections
  void processDplusMixedEvent(CollsWithCentMult::iterator const& coll,
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
    fillCharmMixedEvent<CandidateType::DplusToPiKPi>(candsDplus);
    fillTrkMixedEvent(tracks, cent);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorCorrelationsReduced, processDplusMixedEvent, "Process Mixed Event for Dplus candidates", false);

  // Ds with ML selections
  void processDsSameEvent(CollsWithCentMult::iterator const& coll,
                          TracksData const& tracks,
                          CandDsData const&)
  {
    auto candsDsToPiKK = selectedDsToPiKK->sliceByCached(aod::hf_cand::collisionId, coll.globalIndex(), cache);
    auto candsDsToKKPi = selectedDsToKKPi->sliceByCached(aod::hf_cand::collisionId, coll.globalIndex(), cache);
    if (forceCharmInCollision && candsDsToPiKK.size() < 1 && candsDsToKKPi.size() < 1) {
      return;
    }
    float cent{-1.}, mult{-1.};
    if (!checkCollision(coll, cent, mult)) {
      return;
    }
    rowCollisions(mult, coll.numContrib(), cent, coll.posZ());
    fillSameEvent<CandidateType::DsToPiKK>(candsDsToPiKK, tracks, cent);
    fillSameEvent<CandidateType::DsToKKPi>(candsDsToKKPi, tracks, cent);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorCorrelationsReduced, processDsSameEvent, "Process Same Event for Ds candidates", false);

  // Ds with ML selections
  void processDsMixedEvent(CollsWithCentMult::iterator const& coll,
                           TracksData const& tracks,
                           CandDsData const&)
  {
    auto candsDsToPiKK = selectedDsToPiKK->sliceByCached(aod::hf_cand::collisionId, coll.globalIndex(), cache);
    auto candsDsToKKPi = selectedDsToKKPi->sliceByCached(aod::hf_cand::collisionId, coll.globalIndex(), cache);
    if (forceCharmInCollision && candsDsToPiKK.size() < 1 && candsDsToKKPi.size() < 1) {
      return;
    }
    float cent{-1.}, mult{-1.};
    if (!checkCollision(coll, cent, mult)) {
      return;
    }
    rowCollisions(mult, coll.numContrib(), cent, coll.posZ());
    fillCharmMixedEvent<CandidateType::DsToPiKK>(candsDsToPiKK);
    fillCharmMixedEvent<CandidateType::DsToKKPi>(candsDsToKKPi);
    fillTrkMixedEvent(tracks, cent);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorCorrelationsReduced, processDsMixedEvent, "Process Mixed Event for Ds candidates", false);

  // D0 with ML selections
  void processD0SameEvent(CollsWithCentMult::iterator const& coll,
                          TracksData const& tracks,
                          CandD0Data const&)
  {
    auto candsD0ToPiK = selectedD0ToPiK->sliceByCached(aod::hf_cand::collisionId, coll.globalIndex(), cache);
    auto candsD0ToKPi = selectedD0ToKPi->sliceByCached(aod::hf_cand::collisionId, coll.globalIndex(), cache);
    if (forceCharmInCollision && candsD0ToPiK.size() < 1 && candsD0ToKPi.size() < 1) {
      return;
    }
    float cent{-1.}, mult{-1.};
    if (!checkCollision(coll, cent, mult)) {
      return;
    }
    rowCollisions(mult, coll.numContrib(), cent, coll.posZ());
    fillSameEvent<CandidateType::D0ToPiK>(candsD0ToPiK, tracks, cent);
    fillSameEvent<CandidateType::D0ToKPi>(candsD0ToKPi, tracks, cent);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorCorrelationsReduced, processD0SameEvent, "Process Same Event for D0 candidates", false);

  // D0 with ML selections
  void processD0MixedEvent(CollsWithCentMult::iterator const& coll,
                           TracksData const& tracks,
                           CandD0Data const&)
  {
    auto candsD0ToPiK = selectedD0ToPiK->sliceByCached(aod::hf_cand::collisionId, coll.globalIndex(), cache);
    auto candsD0ToKPi = selectedD0ToKPi->sliceByCached(aod::hf_cand::collisionId, coll.globalIndex(), cache);
    if (forceCharmInCollision && candsD0ToPiK.size() < 1 && candsD0ToKPi.size() < 1) {
      return;
    }
    float cent{-1.}, mult{-1.};
    if (!checkCollision(coll, cent, mult)) {
      return;
    }
    rowCollisions(mult, coll.numContrib(), cent, coll.posZ());
    fillCharmMixedEvent<CandidateType::D0ToPiK>(candsD0ToPiK);
    fillCharmMixedEvent<CandidateType::D0ToKPi>(candsD0ToKPi);
    fillTrkMixedEvent(tracks, cent);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorCorrelationsReduced, processD0MixedEvent, "Process Mixed Event for D0 candidates", false);

  // Lc with ML selections
  void processLcSameEvent(CollsWithCentMult::iterator const& coll,
                           TracksData const& tracks,
                           CandLcData const&)
  {
    auto candsLcToPKPi = selectedLcToPKPi->sliceByCached(aod::hf_cand::collisionId, coll.globalIndex(), cache);
    if (forceCharmInCollision && candsLcToPKPi.size() < 1) {
      return;
    }
    float cent{-1.}, mult{-1.};
    if (!checkCollision(coll, cent, mult)) {
      return;
    }
    rowCollisions(mult, coll.numContrib(), cent, coll.posZ());
    fillSameEvent<CandidateType::LcToPKPi>(candsLcToPKPi, tracks, cent);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorCorrelationsReduced, processLcSameEvent, "Process Same Event for Lc candidates", false);

  // Lc with ML selections
  void processLcMixedEvent(CollsWithCentMult::iterator const& coll,
                           TracksData const& tracks,
                           CandLcData const&)
  {
    auto candsLcToPKPi = selectedLcToPKPi->sliceByCached(aod::hf_cand::collisionId, coll.globalIndex(), cache);
    if (forceCharmInCollision && candsLcToPKPi.size() < 1) {
      return;
    }
    float cent{-1.}, mult{-1.};
    if (!checkCollision(coll, cent, mult)) {
      return;
    }
    rowCollisions(mult, coll.numContrib(), cent, coll.posZ());
    fillCharmMixedEvent<CandidateType::LcToPKPi>(candsLcToPKPi);
    fillTrkMixedEvent(tracks, cent);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorCorrelationsReduced, processLcMixedEvent, "Process Mixed Event for Lc candidates", false);

  // Hadron Hadron Same Event
  void processHadronHadronSameEvent(CollsWithCentMult::iterator const& coll,
                                    TracksData const& tracks)
  {
    float cent{-1.}, mult{-1.};
    if (!checkCollision(coll, cent, mult)) {
      return;
    }
    rowCollisions(mult, coll.numContrib(), cent, coll.posZ());
    fillSameEvent<CandidateType::Hadron>(tracks, tracks, cent);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorCorrelationsReduced, processHadronHadronSameEvent, "Process Same Event for hadron candidates", false);

  // Hadron Hadron Mixed Event
  void processHadronHadronMixedEvent(CollsWithCentMult::iterator const& coll,
                                     TracksData const& tracks)
  {
    float cent{-1.}, mult{-1.};
    if (!checkCollision(coll, cent, mult)) {
      return;
    }
    rowCollisions(mult, coll.numContrib(), cent, coll.posZ());
    fillTrkMixedEvent(tracks, cent);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorCorrelationsReduced, processHadronHadronMixedEvent, "Process Mixed Event for hadron candidates", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfDerivedDataCreatorCorrelationsReduced>(cfgc)};
}
