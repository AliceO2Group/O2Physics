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

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/HFC/DataModel/DerivedDataCorrelationTables.h"
#include "PWGHF/Utils/utilsEvSelHf.h"

#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_centrality;
using namespace o2::hf_evsel;

using BinningTypeDerivedCent = ColumnBinningPolicy<aod::hf_collisions_reduced::PosZ, aod::hf_collisions_reduced::Centrality>;
using BinningTypeDerivedMult = ColumnBinningPolicy<aod::hf_collisions_reduced::PosZ, aod::hf_collisions_reduced::Multiplicity>;

enum CandType {
  DplusToPiKPi = 0,
  DsToKKPi,
  DsToPiKK,
  D0ToPiK,
  D0ToKPi,
  Hadron
};

enum PoolBinningPolicy {
  Centrality = 0,
  Multiplicity
};

/// Code to select collisions with at least one Ds meson
struct HfDerivedDataCreatorCorrelationsReduced {
  Produces<aod::HfcRedFlowColls> rowCollisions;
  Produces<aod::HfcRedCharmTrigs> rowCharmCandidates;
  Produces<aod::HfcRedTrkAssocs> rowAssocTrackReduced;
  Produces<aod::HfcRedSEChHads> entryCharmHadSEPair;
  Produces<aod::HfcRedSEHadHads> entryHadHadSEPair;

  Configurable<int> centEstimator{"centEstimator", 2, "Centrality estimation (FT0A: 1, FT0C: 2, FT0M: 3, FV0A: 4)"};
  Configurable<int> selectionFlag{"selectionFlag", 15, "Selection Flag for hadron (ML score tables are required to run the task)"};
  Configurable<bool> forceCharmInCollision{"forceCharmInCollision", true, "Flag to force charm in collision"};
  Configurable<bool> useCentMixing{"useCentMixing", true, "Flag to use centrality mixing"};
  Configurable<bool> useMultMixing{"useMultMixing", false, "Flag to use multiplicity mixing"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::vector<int>> classMl{"classMl", {0, 2}, "Indexes of BDT scores to be stored. Two indexes max."};
  Configurable<float> yCandMax{"yCandMax", 0.8, "max. cand. rapidity"};
  Configurable<float> ptCandMin{"ptCandMin", 1., "min. cand. pT"};
  Configurable<float> ptCandMax{"ptCandMax", 24., "max. cand. pT"};
  Configurable<int> tpcNClsCrossedRowsMin{"tpcNClsCrossedRowsMin", 70, "min. TPC crossed rows for associated tracks"};
  Configurable<float> etaTrackMax{"etaTrackMax", 1., "max. track eta"};
  Configurable<float> ptTrackMin{"ptTrackMin", 0.15, "min. track pT"};
  Configurable<float> ptTrackMax{"ptTrackMax", 5., "max. track pT"};
  Configurable<float> dcaXYTrackMax{"dcaXYTrackMax", 1., "max. track DCA XY"};
  Configurable<float> dcaZTrackMax{"dcaZTrackMax", 1., "max. track DCA Z"};
  Configurable<float> deltaEtaAbsMin{"deltaEtaAbsMin", 0.5, "min. pair delta eta"};
  Configurable<float> deltaEtaAbsMax{"deltaEtaAbsMax", 2., "max. pair delta eta"};
  Configurable<float> downSampleTracksFactorSE{"downSampleTracksFactorSE", 1., "Fraction of associated tracks to keep"};
  Configurable<float> ptMaxForDownSampleSE{"ptMaxForDownSampleSE", 10., "Maximum pt for the application of the downsampling factor"};
  Configurable<float> centMaxForDownSampleSE{"centMaxForDownSampleSE", 10., "Maximum centrality for the application of the downsampling factor"};
  Configurable<float> downSampleTracksFactorME{"downSampleTracksFactorME", 1., "Fraction of associated tracks to keep"};
  Configurable<float> ptMaxForDownSampleME{"ptMaxForDownSampleME", 10., "Maximum pt for the application of the downsampling factor"};
  Configurable<float> centMaxForDownSampleME{"centMaxForDownSampleME", 10., "Maximum centrality for the application of the downsampling factor"};
  Configurable<std::vector<double>> binsPtTrig{"binsPtTrig", std::vector<double>{1., 3., 5., 8., 16., 36.}, "pT bin limits for trigger candidates"};
  Configurable<std::vector<double>> binsPtAssoc{"binsPtAssoc", std::vector<double>{0.3, 1., 2., 50.}, "pT bin limits for associated particles"};
  Configurable<bool> fillSparses{"fillSparses", true, "Fill sparse histograms"};
  Configurable<bool> fillTables{"fillTables", false, "Fill tables"};

  HfHelper hfHelper;
  HfEventSelection hfEvSel; // event selection and monitoring
  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb;
  SliceCache cache;

  double massCharm{0.};

  using CollsWithCentMult = soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFV0As>;
  using CandDsData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfMlDsToKKPi>>;
  using CandDplusData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfMlDplusToPiKPi>>;
  using CandD0Data = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfMlD0>>;
  using TracksData = soa::Filtered<soa::Join<aod::TracksWDca, aod::TrackSelection, aod::TracksExtra>>;

  Filter filterSelectDsCandidates = aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlag || aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlag;
  Filter filterSelectDplusCandidates = aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlag;
  Filter filterSelectD0Candidates = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlag || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlag;
  Filter filterSelectTrackData = (nabs(aod::track::eta) < etaTrackMax) && (aod::track::pt > ptTrackMin) && (aod::track::pt < ptTrackMax) && (nabs(aod::track::dcaXY) < dcaXYTrackMax) && (nabs(aod::track::dcaZ) < dcaZTrackMax);

  Preslice<CandDsData> candsDsPerColl = aod::hf_cand::collisionId;
  Preslice<CandDplusData> candsDplusPerColl = aod::hf_cand::collisionId;
  Preslice<CandD0Data> candsD0PerColl = aod::hf_cand::collisionId;
  Preslice<TracksData> trackIndicesPerColl = aod::track::collisionId;

  Partition<CandDsData> selectedDsToKKPi = aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlag;
  Partition<CandDsData> selectedDsToPiKK = aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlag;
  Partition<CandD0Data> selectedD0ToPiK = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlag;
  Partition<CandD0Data> selectedD0ToKPi = aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlag;

  ConfigurableAxis zPoolBins{"zPoolBins", {VARIABLE_WIDTH, -10.0, -2.5, 2.5, 10.0}, "Z vertex position pools"};
  ConfigurableAxis multPoolBins{"multPoolBins", {VARIABLE_WIDTH, 0., 900., 1800., 6000.}, "Event multiplicity pools (FT0M)"};
  ConfigurableAxis centPoolBins{"centPoolBins", {VARIABLE_WIDTH, 0., 10., 20., 30.}, "Event centrality pools"};
  ConfigurableAxis binsInvMass{"binsInvMass", {300, 1.6, 2.2}, ""};
  ConfigurableAxis binsMultFT0M{"binsMultFT0M", {100, 0., 10000.}, "Multiplicity as FT0M signal amplitude"};
  ConfigurableAxis binsCent{"binsCent", {100, 0., 100.}, "Centrality bins"};
  ConfigurableAxis binsPosZ{"binsPosZ", {100, -10., 10.}, "Primary vertex z coordinate"};
  ConfigurableAxis binsEta{"binsEta", {50, -2., 2.}, "Eta bins"};
  ConfigurableAxis binsPhi{"binsPhi", {64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}, "Phi bins"};
  ConfigurableAxis binsDeltaEta{"binsDeltaEta", {100, -2., 2.}, "Delta Eta bins"};
  ConfigurableAxis binsDeltaPhi{"binsDeltaPhi", {64, -3., 3.}, "Delta Phi bins"};
  ConfigurableAxis binsPoolBin{"binsPoolBin", {9, 0., 9.}, "PoolBin"};
  ConfigurableAxis binsMlOne{"binsMlOne", {100, 0., 1.}, ""};
  ConfigurableAxis binsMlTwo{"binsMlTwo", {100, 0., 1.}, ""};

  HistogramRegistry registry{"registry", {}};

  ColumnBinningPolicy<aod::hf_collisions_reduced::PosZ, aod::hf_collisions_reduced::Centrality> corrBinningCent{{zPoolBins, centPoolBins}, true};
  ColumnBinningPolicy<aod::hf_collisions_reduced::PosZ, aod::hf_collisions_reduced::Multiplicity> corrBinningMult{{zPoolBins, multPoolBins}, true};

  void init(InitContext&)
  {
    if (doprocessDplusSameEvent || doprocessDplusMixedEvent) {
      massCharm = o2::constants::physics::MassDPlus;
    } else if (doprocessDsSameEvent || doprocessDsMixedEvent) {
      massCharm = o2::constants::physics::MassDS;
    } else if (doprocessD0SameEvent || doprocessD0MixedEvent) {
      massCharm = o2::constants::physics::MassD0;
    } else {
      LOG(fatal) << "No decay channel selected to process";
    }

    hfEvSel.addHistograms(registry); // collision monitoring
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    const AxisSpec axisInvMass{binsInvMass, "Inv. mass (GeV/#it{c}^{2})"};
    const AxisSpec axisCent = {binsCent, "Centrality"};
    const AxisSpec axisMultFT0M = {binsMultFT0M, "MultiplicityFT0M"};
    const AxisSpec axisPosZ = {binsPosZ, "PosZ"};
    const AxisSpec axisPoolBin = {binsPoolBin, "PoolBin"};
    const AxisSpec axisEta = {binsEta, "#it{#eta}"};
    const AxisSpec axisPhi = {binsPhi, "#it{#varphi}"};
    const AxisSpec axisDeltaEta = {binsDeltaEta, "#Delta#it{#eta}"};
    const AxisSpec axisDeltaPhi = {binsDeltaPhi, "#Delta#it{#varphi}"};
    const AxisSpec axisPtTrig = {(std::vector<double>)binsPtTrig, "#it{p}_{T} Trig (GeV/#it{c})"};
    const AxisSpec axisPtAssoc = {(std::vector<double>)binsPtAssoc, "#it{p}_{T} Assoc (GeV/#it{c})"};
    const AxisSpec axisMlOne{binsMlOne, "bdtScore0"};
    const AxisSpec axisMlTwo{binsMlTwo, "bdtScore1"};

    // Histograms for data analysis
    if (useCentMixing) {
      registry.add("hCent", "Centrality", {HistType::kTH2F, {{axisCent}, {axisPoolBin}}});
    } else {
      registry.add("hMultFT0M", "Multiplicity FT0M", {HistType::kTH2F, {{axisMultFT0M}, {axisPoolBin}}});
    }
    registry.add("hZVtx", "z vertex", {HistType::kTH2F, {{axisPosZ}, {axisPoolBin}}});
    registry.add("hCollisionPoolBin", "Collision pool bin", {HistType::kTH1F, {axisPoolBin}});
    registry.add("hPoolBinTrig", "Trigger candidates pool bin", {HistType::kTH1F, {axisPoolBin}});
    registry.add("hPhiVsPtTrig", "Trigger candidates phiVsPt", {HistType::kTH2F, {{axisPhi}, {axisPtTrig}}});
    registry.add("hEtaVsPtTrig", "Trigger candidates etaVsPt", {HistType::kTH2F, {{axisEta}, {axisPtTrig}}});
    registry.add("hPhiVsPtTrigAssoc", "Associated particles phiVsPt", {HistType::kTH3F, {{axisPhi}, {axisPtTrig}, {axisPtAssoc}}});
    registry.add("hEtaVsPtTrigAssoc", "Associated particles etaVsPt", {HistType::kTH3F, {{axisEta}, {axisPtTrig}, {axisPtAssoc}}});
    registry.add("hPhiVsPtAssoc", "Associated particles phiVsPt", {HistType::kTH2F, {{axisPhi}, {axisPtAssoc}}});
    registry.add("hEtaVsPtAssoc", "Associated particles etaVsPt", {HistType::kTH2F, {{axisEta}, {axisPtAssoc}}});
    registry.add("hPoolBinAssoc", "Associated particles pool bin", {HistType::kTH1F, {axisPoolBin}});

    // if (fillSparses) {
    //   std::vector<AxisSpec> axes = {axisPtTrig, axisPtAssoc, axisDeltaEta, axisDeltaPhi, axisPoolBin};
    //   if (doprocessSameEventDplusWCentMix) {
    //     registry.add("hSparseCorrelationsSEHadHad", "THn for SE Had-Had correlations", HistType::kTHnSparseF, axes);
    //   } else if (doprocessMixedEventHadHadWCentMix || doprocessMixedEventHadHadWMultMix) {
    //     registry.add("hSparseCorrelationsMEHadHad", "THn for ME Had-Had correlations", HistType::kTHnSparseF, axes);
    //   } else {
    //     axes.insert(axes.end(), {axisMlOne, axisMlTwo, axisInvMass});
    //     if (doprocessSameEventDplusCharmHadWCentMix || doprocessSameEventCharmHadWMultMix) {
    //       registry.add("hSparseCorrelationsSECharmHad", "THn for SE Charm-Had correlations", HistType::kTHnSparseF, axes);
    //     } else if (doprocessMixedEventCharmHadWCentMix || doprocessMixedEventCharmHadWMultMix) {
    //       registry.add("hSparseCorrelationsMECharmHad", "THn for ME Charm-Had correlations", HistType::kTHnSparseF, axes);
    //     }
    //   }
    // }
  }; // end init

  /// Get charm hadron candidate mass
  /// \param candidate is the charm hadron candidate
  template <CandType candType, typename TCand>
  double getCandMass(const TCand& candidate)
  {
    if constexpr (candType == CandType::DsToKKPi) {
      return hfHelper.invMassDsToKKPi(candidate);
    }
    if constexpr (candType == CandType::DsToPiKK) {
      return hfHelper.invMassDsToPiKK(candidate);
    }
    if constexpr (candType == CandType::DplusToPiKPi) {
      return hfHelper.invMassDplusToPiKPi(candidate);
    }
    if constexpr (candType == CandType::D0ToPiK) {
      return hfHelper.invMassD0ToPiK(candidate);
    }
    if constexpr (candType == CandType::D0ToKPi) {
      return hfHelper.invMassD0barToKPi(candidate);
    }
    return -1.;
  }

  /// Get charm hadron bdt scores
  /// \param candidate is the charm hadron candidate
  template <CandType candType, typename TCand>
  std::vector<float> getCandMlScores(const TCand& candidate)
  {
    std::vector<float> outputMl{-999., -999.};
    if constexpr (candType == CandType::DsToKKPi) {
      for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
        outputMl[iclass] = candidate.mlProbDsToKKPi()[classMl->at(iclass)];
      }
    }
    if constexpr (candType == CandType::DsToPiKK) {
      for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
        outputMl[iclass] = candidate.mlProbDsToPiKK()[classMl->at(iclass)];
      }
    }
    if constexpr (candType == CandType::DplusToPiKPi) {
      for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
        outputMl[iclass] = candidate.mlProbDplusToPiKPi()[classMl->at(iclass)];
      }
    }
    if constexpr (candType == CandType::D0ToPiK) {
      for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
        outputMl[iclass] = candidate.mlProbD0()[classMl->at(iclass)];
      }
    }
    if constexpr (candType == CandType::D0ToKPi) {
      for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
        outputMl[iclass] = candidate.mlProbD0bar()[classMl->at(iclass)];
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
    if (collRejMask != 0) {
      return false;
    }
    return true;
  }

  template <typename TCand, typename TTrack>
  bool checkDaughterTrack(TTrack const& track, TCand const& cand)
  {
    if constexpr ((requires { cand.prong2Id(); })) {  // Check 3-prong
      return (track.globalIndex() == cand.prong0Id() || track.globalIndex() == cand.prong1Id() || track.globalIndex() == cand.prong2Id());
    } else {  // Check 2-prong
      return (track.globalIndex() == cand.prong0Id() || track.globalIndex() == cand.prong1Id());
    }
  }

  /// Fill Charm-Hadron correlation table and sparse
  /// \param trigCand is the trigger charm hadron candidate
  /// \param assTrk is the associated hadron track
  /// \param poolBin is the pool bin of the collision
  template <CandType candType, typename TTrigCand, typename TTrack>
  void fillCharmHadSE(TTrigCand const& trigCand,
                      TTrack const& assTrk,
                      const int poolBin)
  {
    double deltaEta = assTrk.eta() - trigCand.eta();
    double deltaPhi = RecoDecay::constrainAngle(assTrk.phi() - trigCand.phi(), -o2::constants::math::PIHalf);
    std::vector<float> outputMl = getCandMlScores<candType>(trigCand);
    double massCand = getCandMass<candType>(trigCand);
    if (fillTables) {
      entryCharmHadSEPair(deltaEta, deltaPhi, poolBin, trigCand.pt(), massCand, outputMl[0], outputMl[1],
                          assTrk.pt(), assTrk.tpcNClsCrossedRows(), assTrk.itsClusterMap(), assTrk.itsNCls(), assTrk.dcaXY(), assTrk.dcaZ());
    }
    // if (fillSparses) {
    //   registry.fill(HIST("hSparseCorrelationsSECharmHad"), trigCand.pt(), assTrk.pt(),
    //                 deltaEta, deltaPhi, poolBin, trigCand.bdtScore0(),
    //                 trigCand.bdtScore1(), trigCand.invMassCand());
    // }
  }

  /// Fill Hadron-Hadron correlation table and sparse
  /// \param trigTrack is the trigger hadron candidate
  /// \param assTrk is the associated hadron track
  /// \param poolBin is the pool bin of the collision
  template <typename TCand>
  void fillHadHadSE(TCand const& trigTrack,
                    TCand const& assTrk,
                    const int poolBin)
  {
    double deltaEta = assTrk.eta() - trigTrack.eta();
    double deltaPhi = RecoDecay::constrainAngle(assTrk.phi() - trigTrack.phi(), -o2::constants::math::PIHalf);
    if (fillTables) {
      entryHadHadSEPair(deltaEta, deltaPhi, poolBin,
                        trigTrack.pt(), assTrk.tpcNClsCrossedRows(), assTrk.itsClusterMap(), assTrk.itsNCls(), assTrk.dcaXY(), assTrk.dcaZ(),
                        assTrk.pt(), assTrk.tpcNClsCrossedRows(), assTrk.itsClusterMap(), assTrk.itsNCls(), assTrk.dcaXY(), assTrk.dcaZ());
    }
    // if (fillSparses) {
    //   registry.fill(HIST("hSparseCorrelationsSEHadHad"), trigCand.pt(), assTrk.pt(), deltaEta, deltaPhi, poolBin);
    // }
  }

  template <CandType candType, typename TTrigCands, typename TAssocTracks>
  void fillSameEvent(TTrigCands const& trigCands,
                     TAssocTracks const& assTrks,
                     const int poolBin,
                     const float collCentrality)
  {

    for (const auto& trigCand : trigCands) {
      registry.fill(HIST("hPoolBinTrig"), poolBin);
      registry.fill(HIST("hPhiVsPtTrig"), RecoDecay::constrainAngle(trigCand.phi(), -o2::constants::math::PIHalf), trigCand.pt());
      registry.fill(HIST("hEtaVsPtTrig"), trigCand.eta(), trigCand.pt());
      for (const auto& assTrk : assTrks) {
        if (assTrk.tpcNClsCrossedRows() < tpcNClsCrossedRowsMin) {
          continue;
        }
        double deltaEta = assTrk.eta() - trigCand.eta();
        if (std::abs(deltaEta) < deltaEtaAbsMin || std::abs(deltaEta) > deltaEtaAbsMax) {
          continue;
        }
        if constexpr (candType == CandType::Hadron) {
          if (assTrk.globalIndex() == trigCand.globalIndex()) {
            continue; // skip self-correlation for hadron-hadron
          }
          if (downSampleTracksFactorSE < 1.) {
            float pseudoRndm = assTrk.pt() * 1000. - static_cast<int64_t>(assTrk.pt() * 1000);
            if (assTrk.pt() < ptMaxForDownSampleSE && collCentrality < centMaxForDownSampleSE && pseudoRndm >= downSampleTracksFactorSE) {
              continue;
            }
          }
        } else {
          if (checkDaughterTrack(assTrk, trigCand)) {
            continue; // skip daughter tracks for charm-hadron
          }
        }

        registry.fill(HIST("hPoolBinAssoc"), poolBin);
        registry.fill(HIST("hPhiVsPtTrigAssoc"), RecoDecay::constrainAngle(assTrk.phi(), -o2::constants::math::PIHalf), trigCand.pt(), assTrk.pt());
        registry.fill(HIST("hEtaVsPtAssoc"), assTrk.eta(), trigCand.pt(), assTrk.pt());
        if constexpr (candType == CandType::Hadron) {
          fillHadHadSE(trigCand, assTrk, poolBin);
        } else {
          fillCharmHadSE<candType>(trigCand, assTrk, poolBin);
        }
      }
    }
  }

  template <CandType candType, typename TTrigCands>
  void fillCharmMixedEvent(TTrigCands const& trigCands,
                           const int poolBin)
  {
    for (const auto& trigCand : trigCands) {
      registry.fill(HIST("hPoolBinTrig"), poolBin);
      registry.fill(HIST("hPhiVsPtTrig"), RecoDecay::constrainAngle(trigCand.phi(), -o2::constants::math::PIHalf), trigCand.pt());
      registry.fill(HIST("hEtaVsPtTrig"), trigCand.eta(), trigCand.pt());

      rowCharmCandidates(rowCollisions.lastIndex(), trigCand.phi(), trigCand.eta(), trigCand.pt(),
                         getCandMass<candType>(trigCand), getCandMlScores<candType>(trigCand)[0],
                         getCandMlScores<candType>(trigCand)[1]);
    }
  }

  template <typename TAssocTracks>
  void fillTrackMixedEvent(TAssocTracks const& assTrks,
                           const int poolBin,
                           const float collCentrality)
  {
    bool first = true;
    for (const auto& assTrk : assTrks) {
      if (assTrk.tpcNClsCrossedRows() < tpcNClsCrossedRowsMin) {
        continue;
      }
      if (!first && downSampleTracksFactorME < 1.) {
        float pseudoRndm = assTrk.pt() * 1000. - static_cast<int64_t>(assTrk.pt() * 1000);
        if (assTrk.pt() < ptMaxForDownSampleME && collCentrality < centMaxForDownSampleME && pseudoRndm >= downSampleTracksFactorME) {
          continue;
        }
      }
      first = false;
      registry.fill(HIST("hPoolBinAssoc"), poolBin);
      registry.fill(HIST("hPhiVsPtAssoc"), RecoDecay::constrainAngle(assTrk.phi(), -o2::constants::math::PIHalf), assTrk.pt());
      registry.fill(HIST("hEtaVsPtAssoc"), assTrk.eta(), assTrk.pt());
      rowAssocTrackReduced(rowCollisions.lastIndex(), assTrk.phi(), assTrk.eta(), 
                           assTrk.pt(), assTrk.tpcNClsCrossedRows(), assTrk.itsClusterMap(), 
                           assTrk.itsNCls(), assTrk.dcaXY(), assTrk.dcaZ());
    }
  }

  /// Get the binning pool associated to the collision
  /// \param collision is the collision
  /// \param corrBinning is the binning policy for the correlation
  template <PoolBinningPolicy binningPolicy, bool fillHistos, typename TColl>
  int getPoolBin(const TColl& collision, const float cent, const float mult)
  {
    int poolBin{0};
    if constexpr (binningPolicy == PoolBinningPolicy::Centrality) {
      poolBin = corrBinningCent.getBin(std::make_tuple(collision.posZ(), cent));
      if constexpr (fillHistos) {
        registry.fill(HIST("hCent"), cent, poolBin);
      }
    } else if constexpr (binningPolicy == PoolBinningPolicy::Multiplicity) {
      poolBin = corrBinningMult.getBin(std::make_tuple(collision.posZ(), mult));
      if constexpr (fillHistos) {
        registry.fill(HIST("hMultFT0M"), mult, poolBin);
      }
    }
    return poolBin;
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
    int poolBin = useCentMixing ? getPoolBin<PoolBinningPolicy::Centrality, true>(coll, cent, mult) : getPoolBin<PoolBinningPolicy::Multiplicity, true>(coll, cent, mult);
    registry.fill(HIST("hCollisionPoolBin"), poolBin);
    registry.fill(HIST("hZVtx"), coll.posZ(), poolBin);
    fillSameEvent<CandType::DplusToPiKPi>(candsDplus, tracks, poolBin, cent);
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
    int poolBin = useCentMixing ? getPoolBin<PoolBinningPolicy::Centrality, true>(coll, cent, mult) : getPoolBin<PoolBinningPolicy::Multiplicity, true>(coll, cent, mult);
    registry.fill(HIST("hCollisionPoolBin"), poolBin);
    registry.fill(HIST("hZVtx"), coll.posZ(), poolBin);
    fillCharmMixedEvent<CandType::DplusToPiKPi>(candsDplus, poolBin);
    fillTrackMixedEvent(tracks, poolBin, cent);
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
    int poolBin = useCentMixing ? getPoolBin<PoolBinningPolicy::Centrality, true>(coll, cent, mult) : getPoolBin<PoolBinningPolicy::Multiplicity, true>(coll, cent, mult);
    registry.fill(HIST("hCollisionPoolBin"), poolBin);
    registry.fill(HIST("hZVtx"), coll.posZ(), poolBin);
    fillSameEvent<CandType::DsToPiKK>(candsDsToPiKK, tracks, poolBin, cent);
    fillSameEvent<CandType::DsToKKPi>(candsDsToKKPi, tracks, poolBin, cent);
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
    int poolBin = useCentMixing ? getPoolBin<PoolBinningPolicy::Centrality, true>(coll, cent, mult) : getPoolBin<PoolBinningPolicy::Multiplicity, true>(coll, cent, mult);
    registry.fill(HIST("hCollisionPoolBin"), poolBin);
    registry.fill(HIST("hZVtx"), coll.posZ(), poolBin);
    fillCharmMixedEvent<CandType::DsToPiKK>(candsDsToPiKK, poolBin);
    fillCharmMixedEvent<CandType::DsToKKPi>(candsDsToKKPi, poolBin);
    fillTrackMixedEvent(tracks, poolBin, cent);
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
    int poolBin = useCentMixing ? getPoolBin<PoolBinningPolicy::Centrality, true>(coll, cent, mult) : getPoolBin<PoolBinningPolicy::Multiplicity, true>(coll, cent, mult);
    registry.fill(HIST("hCollisionPoolBin"), poolBin);
    registry.fill(HIST("hZVtx"), coll.posZ(), poolBin);
    fillSameEvent<CandType::D0ToPiK>(candsD0ToPiK, tracks, poolBin, cent);
    fillSameEvent<CandType::D0ToKPi>(candsD0ToKPi, tracks, poolBin, cent);
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
    int poolBin = useCentMixing ? getPoolBin<PoolBinningPolicy::Centrality, true>(coll, cent, mult) : getPoolBin<PoolBinningPolicy::Multiplicity, true>(coll, cent, mult);
    registry.fill(HIST("hCollisionPoolBin"), poolBin);
    registry.fill(HIST("hZVtx"), coll.posZ(), poolBin);
    fillCharmMixedEvent<CandType::D0ToPiK>(candsD0ToPiK, poolBin);
    fillCharmMixedEvent<CandType::D0ToKPi>(candsD0ToKPi, poolBin);
    fillTrackMixedEvent(tracks, poolBin, cent);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorCorrelationsReduced, processD0MixedEvent, "Process Mixed Event for D0 candidates", false);

  // Hadron Hadron Same Event
  void processHadronHadronSameEvent(CollsWithCentMult::iterator const& coll,
                                    TracksData const& tracks)
  {
    float cent{-1.}, mult{-1.};
    if (!checkCollision(coll, cent, mult)) {
      return;
    }
    rowCollisions(mult, coll.numContrib(), cent, coll.posZ());
    int poolBin = useCentMixing ? getPoolBin<PoolBinningPolicy::Centrality, true>(coll, cent, mult) : getPoolBin<PoolBinningPolicy::Multiplicity, true>(coll, cent, mult);
    registry.fill(HIST("hCollisionPoolBin"), poolBin);
    registry.fill(HIST("hZVtx"), coll.posZ(), poolBin);
    fillSameEvent<CandType::Hadron>(tracks, tracks, poolBin, cent);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorCorrelationsReduced, processHadronHadronSameEvent, "Process Same Event for hadron candidates", true);

  // Hadron Hadron Mixed Event
  void processHadronHadronMixedEvent(CollsWithCentMult::iterator const& coll,
                                     TracksData const& tracks)
  {
    float cent{-1.}, mult{-1.};
    if (!checkCollision(coll, cent, mult)) {
      return;
    }
    rowCollisions(mult, coll.numContrib(), cent, coll.posZ());
    int poolBin = useCentMixing ? getPoolBin<PoolBinningPolicy::Centrality, true>(coll, cent, mult) : getPoolBin<PoolBinningPolicy::Multiplicity, true>(coll, cent, mult);
    registry.fill(HIST("hCollisionPoolBin"), poolBin);
    registry.fill(HIST("hZVtx"), coll.posZ(), poolBin);
    fillTrackMixedEvent(tracks, poolBin, cent);
  }
  PROCESS_SWITCH(HfDerivedDataCreatorCorrelationsReduced, processHadronHadronMixedEvent, "Process Mixed Event for hadron candidates", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfDerivedDataCreatorCorrelationsReduced>(cfgc)};
}
