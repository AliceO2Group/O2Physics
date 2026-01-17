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

/// \file correlatorFlowCharmHadronsReduced.cxx
/// \brief CharmHadrons-Hadrons correlator tree creator for data analyses
/// \author Marcello Di Costanzo <marcello.di.costanzo@cern.ch>, Politecnico and INFN Torino

#include "PWGHF/HFC/DataModel/DerivedDataCorrelationTables.h"

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/ASoA.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/SliceCache.h>
#include <Framework/runDataProcessing.h>

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using BinningCentPosZ = ColumnBinningPolicy<aod::hf_collisions_reduced::PosZ, aod::hf_collisions_reduced::Centrality>;
using BinningMultPosZ = ColumnBinningPolicy<aod::hf_collisions_reduced::PosZ, aod::hf_collisions_reduced::Multiplicity>;

/// Get charm candidate or hadron track pT
/// \param track is the candidate
template <typename TTrack>
double getPt(const TTrack& track)
{
  if constexpr (requires { track.ptAssoc(); }) {
    return track.ptAssoc();
  } else {
    return track.ptTrig();
  }
}

/// Get charm candidate or hadron track eta
/// \param track is the candidate
template <typename TTrack>
double getEta(const TTrack& track)
{
  if constexpr (requires { track.etaAssoc(); }) {
    return track.etaAssoc();
  } else {
    return track.etaTrig();
  }
}

/// Get charm candidate or hadron track phi
/// \param track is the candidate
template <typename TTrack>
double getPhi(const TTrack& track)
{
  if constexpr (requires { track.phiAssoc(); }) {
    return track.phiAssoc();
  } else {
    return track.phiTrig();
  }
}

struct HfCorrelatorFlowCharmHadronsReduced {
  // Produces<aod::HfcRedSEChHads> rowPairSECharmHads; //! Correlation pairs information Same Event
  // Produces<aod::HfcRedMEChHads> rowPairMECharmHads; //! Correlation pairs information Mixed Event
  // Produces<aod::HfcRedSEHadHads> rowPairSEHadHads;  //! Correlation pairs information Same Event
  // Produces<aod::HfcRedMEHadHads> rowPairMEHadHads;  //! Correlation pairs information Mixed Event
  // Produces<aod::HfcRedCollInfos> rowCollInfos;      //! Collision info

  Configurable<bool> fillSparses{"fillSparses", true, "Fill sparse histograms"};
  Configurable<bool> fillTables{"fillTables", false, "Fill tables"};
  Configurable<int> numberEventsMixed{"numberEventsMixed", 5, "Number of events mixed in ME process"};
  Configurable<std::vector<double>> binsPtTrig{"binsPtTrig", std::vector<double>{0., 3., 5., 8., 16., 36.}, "pT bin limits for trigger candidates"};
  Configurable<std::vector<double>> bdtScore0PtMaxs{"bdtScore0PtMaxs", std::vector<double>{0.1, 0.1, 0.1, 0.1, 0.1}, "pT-differential maximum score 0 for charm candidates"};
  Configurable<std::vector<double>> bdtScore1PtMins{"bdtScore1PtMins", std::vector<double>{0.1, 0.1, 0.1, 0.1, 0.1}, "pT-differential minimum score 1 for charm candidates"};
  Configurable<std::vector<double>> binsPtAssoc{"binsPtAssoc", std::vector<double>{0.3, 1., 2., 50.}, "pT bin limits for associated particles"};
  Configurable<float> centralityMin{"centralityMin", 0, "min. centrality"};
  Configurable<float> centralityMax{"centralityMax", 10., "max. centrality"};
  Configurable<float> deltaEtaAbsMin{"deltaEtaAbsMin", 0.5, "min. pair delta eta"};
  Configurable<float> deltaEtaAbsMax{"deltaEtaAbsMax", 2., "max. pair delta eta"};
  Configurable<float> dcaXYTrackMax{"dcaXYTrackMax", 1., "max. track DCA XY"};
  Configurable<float> dcaZTrackMax{"dcaZTrackMax", 1., "max. track DCA Z"};
  Configurable<int> tpcCrossedRowsMin{"tpcCrossedRowsMin", 1, "min. TPC crossed rows"};
  Configurable<int> itsNClsMin{"itsNClsMin", 1, "min. ITS clusters"};
  Configurable<float> downSamplePairs{"downSamplePairs", 1., "Fraction of pairs to keep"};
  Configurable<float> ptMaxForDownSample{"ptMaxForDownSample", 10., "Maximum pt for the application of the downsampling factor"};
  Configurable<float> centMaxForDownSample{"centMaxForDownSample", 10., "Maximum centrality for the application of the downsampling factor"};

  SliceCache cache;

  int poolBins{0};

  using SameEvtPairsChHad = soa::Filtered<soa::Join<aod::HfcRedSEChBases, aod::HfcRedAssTracks>>;
  using SameEvtPairsHadHad = soa::Filtered<soa::Join<aod::HfcRedSEHadBases, aod::HfcRedAssTracks>>;
  using AssocTracks = soa::Filtered<soa::Join<aod::HfcRedAssBases, aod::HfcRedAssTracks>>;
  using TrigCharmCands = soa::Join<aod::HfcRedTrigBases, aod::HfcRedTrigCharms>;

  Filter filterAssocTracks = (nabs(aod::hf_correl_charm_had_reduced::dcaXYAssoc) < dcaXYTrackMax) && (nabs(aod::hf_correl_charm_had_reduced::dcaZAssoc) < dcaZTrackMax) && (aod::hf_correl_charm_had_reduced::nTpcCrossedRowsAssoc > tpcCrossedRowsMin) && (aod::hf_correl_charm_had_reduced::itsNClsAssoc > itsNClsMin);
  Filter filterTrigTracks = (nabs(aod::hf_correl_charm_had_reduced::dcaXYTrig) < dcaXYTrackMax) && (nabs(aod::hf_correl_charm_had_reduced::dcaZTrig) < dcaZTrackMax) && (aod::hf_correl_charm_had_reduced::nTpcCrossedRowsTrig > tpcCrossedRowsMin) && (aod::hf_correl_charm_had_reduced::itsNClsTrig > itsNClsMin);
  Filter filterSameEvtPairs = (nabs(aod::hf_correl_charm_had_reduced::deltaEta) > deltaEtaAbsMin) && (nabs(aod::hf_correl_charm_had_reduced::deltaEta) < deltaEtaAbsMax);

  Preslice<AssocTracks> assocTracksPerCol = aod::hf_correl_charm_had_reduced::hfcRedCorrCollId;
  Preslice<TrigCharmCands> trigCharmCandsPerCol = aod::hf_correl_charm_had_reduced::hfcRedCorrCollId;

  ConfigurableAxis zPoolBins{"zPoolBins", {VARIABLE_WIDTH, -10.0, -2.5, 2.5, 10.0}, "Z vertex position pools"};
  ConfigurableAxis multPoolBins{"multPoolBins", {VARIABLE_WIDTH, 0., 900., 1800., 6000.}, "Event multiplicity pools (FT0M)"};
  ConfigurableAxis centPoolBins{"centPoolBins", {VARIABLE_WIDTH, 0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100}, "Event centrality pools"};
  ConfigurableAxis binsInvMass{"binsInvMass", {300, 1.6, 2.2}, "Invariant mass bins"};
  ConfigurableAxis binsMultFT0M{"binsMultFT0M", {100, 0., 10000.}, "Multiplicity as FT0M signal amplitude"};
  ConfigurableAxis binsCent{"binsCent", {100, 0., 100.}, "Centrality bins"};
  ConfigurableAxis binsPosZ{"binsPosZ", {100, -10., 10.}, "Primary vertex z coordinate"};
  ConfigurableAxis binsEta{"binsEta", {50, -2., 2.}, "Eta bins"};
  ConfigurableAxis binsPhi{"binsPhi", {64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}, "Phi bins"};
  ConfigurableAxis binsDeltaEta{"binsDeltaEta", {100, -2., 2.}, "Delta Eta bins"};
  ConfigurableAxis binsDeltaPhi{"binsDeltaPhi", {64, -3., 3.}, "Delta Phi bins"};
  ConfigurableAxis binsMlOne{"binsMlOne", {100, 0., 1.}, "ML score index 1 bins"};
  ConfigurableAxis binsMlTwo{"binsMlTwo", {100, 0., 1.}, "ML score index 2 bins"};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    if ((doprocessSameEventCharmHadWCentMix && doprocessMixedEventCharmHadWMultMix) ||
        (doprocessSameEventCharmHadWMultMix && doprocessMixedEventCharmHadWCentMix) ||
        (doprocessSameEventHadHadWCentMix && doprocessMixedEventHadHadWMultMix) ||
        (doprocessSameEventHadHadWMultMix && doprocessMixedEventHadHadWCentMix)) {
      LOGP(fatal, "You cannot mix centrality and multiplicity mixing in the same processing! Please check your configuration!");
    }
    if (!fillSparses && !fillTables) {
      LOGP(fatal, "At least one of fillSparses or fillTables must be true!");
    }
    if ((binsPtTrig.value.size() != (bdtScore0PtMaxs.value.size() + 1) || binsPtTrig.value.size() != (bdtScore1PtMins.value.size() + 1))) {
      LOGP(fatal, "The size of bdtScore0PtMaxs and bdtScore1PtMins must be the one of binsPtTrig minus one!");
    }

    if (doprocessSameEventCharmHadWCentMix || doprocessSameEventHadHadWCentMix || doprocessMixedEventCharmHadWCentMix || doprocessMixedEventHadHadWCentMix) {
      poolBins = (centPoolBins->size() - 2) * (zPoolBins->size() - 2);
    } else {
      poolBins = (multPoolBins->size() - 2) * (zPoolBins->size() - 2);
    }

    const AxisSpec axisInvMass{binsInvMass, "Inv. mass (GeV/#it{c}^{2})"};
    const AxisSpec axisCent = {binsCent, "Centrality"};
    const AxisSpec axisMultFT0M = {binsMultFT0M, "MultiplicityFT0M"};
    const AxisSpec axisPosZ = {binsPosZ, "PosZ"};
    const AxisSpec axisPoolBin = {poolBins, 0., static_cast<float>(poolBins), "PoolBin"};
    const AxisSpec axisDeltaEta = {binsDeltaEta, "#Delta#it{#eta}"};
    const AxisSpec axisDeltaPhi = {binsDeltaPhi, "#Delta#it{#varphi}"};
    const AxisSpec axisPtTrig = {(std::vector<double>)binsPtTrig, "#it{p}_{T} Trig (GeV/#it{c})"};
    const AxisSpec axisPtAssoc = {(std::vector<double>)binsPtAssoc, "#it{p}_{T} Assoc (GeV/#it{c})"};
    const AxisSpec axisMlOne{binsMlOne, "bdtScore0"};
    const AxisSpec axisMlTwo{binsMlTwo, "bdtScore1"};

    // Histograms for data analysis
    if (doprocessSameEventCharmHadWCentMix || doprocessSameEventHadHadWCentMix) {
      registry.add("hCentPoolBinSE", "Centrality SE", {HistType::kTH2F, {{axisCent}, {axisPoolBin}}});
    } else if (doprocessSameEventCharmHadWMultMix || doprocessSameEventHadHadWMultMix) {
      registry.add("hMultFT0MPoolBinSE", "Multiplicity FT0M SE", {HistType::kTH2F, {{axisMultFT0M}, {axisPoolBin}}});
    } else if (doprocessMixedEventCharmHadWCentMix || doprocessMixedEventHadHadWCentMix) {
      registry.add("hCentPoolBinME", "Centrality ME", {HistType::kTH2F, {{axisCent}, {axisPoolBin}}});
    } else if (doprocessMixedEventCharmHadWMultMix || doprocessMixedEventHadHadWMultMix) {
      registry.add("hMultFT0MPoolBinME", "Multiplicity FT0M ME", {HistType::kTH2F, {{axisMultFT0M}, {axisPoolBin}}});
    }
    registry.add("hZVtxPoolBinSE", "z vertex SE", {HistType::kTH2F, {{axisPosZ}, {axisPoolBin}}});
    registry.add("hZVtxPoolBinME", "z vertex ME", {HistType::kTH2F, {{axisPosZ}, {axisPoolBin}}});
    registry.add("hPoolBinTrigSE", "Trigger candidates pool bin SE", {HistType::kTH1F, {axisPoolBin}});
    registry.add("hPoolBinTrigME", "Trigger candidates pool bin ME", {HistType::kTH1F, {axisPoolBin}});
    registry.add("hPoolBinAssocSE", "Associated particles pool bin SE", {HistType::kTH1F, {axisPoolBin}});
    registry.add("hPoolBinAssocME", "Associated particles pool bin ME", {HistType::kTH1F, {axisPoolBin}});
    if (fillSparses) {
      std::vector<AxisSpec> axesTrigger = {axisInvMass, axisPtTrig, axisMlOne, axisMlTwo};
      std::vector<AxisSpec> axes = {axisPoolBin, axisPtTrig, axisPtAssoc, axisDeltaEta, axisDeltaPhi};
      if (doprocessSameEventHadHadWCentMix || doprocessSameEventHadHadWMultMix) {
        registry.add("hSparseCorrelationsSEHadHad", "THn for SE Had-Had correlations", HistType::kTHnSparseF, axes);
      } else if (doprocessMixedEventHadHadWCentMix || doprocessMixedEventHadHadWMultMix) {
        registry.add("hSparseCorrelationsMEHadHad", "THn for ME Had-Had correlations", HistType::kTHnSparseF, axes);
      } else {
        axes.insert(axes.end(), {axisInvMass});
        // axes.insert(axes.end(), {axisInvMass, axisMlOne, axisMlTwo});
        if (doprocessSameEventCharmHadWCentMix || doprocessSameEventCharmHadWMultMix || doprocessSameEventCharmHadWCentMixBase) {
          registry.add("hSparseCorrelationsSECharmHad", "THn for SE Charm-Had correlations", HistType::kTHnSparseF, axes);
        } else if (doprocessMixedEventCharmHadWCentMix || doprocessMixedEventCharmHadWMultMix || doprocessMixedEventCharmHadWCentMixBase) {
          registry.add("hSparseCorrelationsMECharmHad", "THn for ME Charm-Had correlations", HistType::kTHnSparseF, axes);
        }
        if (doprocessCharmTriggers) {
          registry.add("hSparseTrigCandsCharm", "THn for Charm trigger candidates", HistType::kTHnSparseF, axesTrigger);
        }
      }
    }
  }

  /// Get the binning pool associated to the collision
  /// \param collision is the collision
  /// \param binPolicy is the binning policy for the correlation
  template <bool IsMixedEvent, typename TColl, typename TBinningType>
  int getPoolBin(const TColl& collision, const TBinningType& binPolicy)
  {
    int poolBin{0};
    if constexpr (std::is_same_v<TBinningType, BinningCentPosZ>) {
      poolBin = binPolicy.getBin(std::make_tuple(collision.posZ(), collision.centrality()));
      if constexpr (IsMixedEvent) {
        registry.fill(HIST("hCentPoolBinME"), collision.centrality(), poolBin);
        registry.fill(HIST("hZVtxPoolBinME"), collision.posZ(), poolBin);
      } else {
        registry.fill(HIST("hCentPoolBinSE"), collision.centrality(), poolBin);
        registry.fill(HIST("hZVtxPoolBinSE"), collision.posZ(), poolBin);
      }
    } else if constexpr (std::is_same_v<TBinningType, BinningMultPosZ>) {
      poolBin = binPolicy.getBin(std::make_tuple(collision.posZ(), collision.multiplicity()));
      if constexpr (IsMixedEvent) {
        registry.fill(HIST("hMultFT0MPoolBinME"), collision.multiplicity(), poolBin);
        registry.fill(HIST("hZVtxPoolBinME"), collision.posZ(), poolBin);
      } else {
        registry.fill(HIST("hMultFT0MPoolBinSE"), collision.multiplicity(), poolBin);
        registry.fill(HIST("hZVtxPoolBinSE"), collision.posZ(), poolBin);
      }
    }
    return poolBin;
  }

  /// Apply pT-differential ML BDT bkg score cut
  /// \param ptTrig is the pT of the charm candidate
  template <typename TCand>
  bool isSelBdtScoreCut(TCand const& cand,
                        double ptTrig)
  {
    for (size_t iPt = 0; iPt < binsPtTrig.value.size() - 1; iPt++) {
      if (ptTrig >= binsPtTrig.value[iPt] && ptTrig < binsPtTrig.value[iPt + 1]) {
        return (cand.bdtScore0Trig() < bdtScore0PtMaxs.value[iPt]) && (cand.bdtScore1Trig() > bdtScore1PtMins.value[iPt]);
      }
    }
    return false;
  }

  /// Save info for Same Event pairs
  /// \param collisions are the selected collisions
  /// \param trigCands are the selected trigger candidates
  /// \param assocTracks are the selected associated tracks
  /// \param binPolicy is the binning policy for the correlation
  template <bool FillSparses, bool FillTables, typename TPair, typename TTrigCand, typename TBinningType>
  void fillSameEvent(TPair const& pair,
                     TTrigCand const& trigCand,
                     TBinningType binPolicy)
  {
    auto collision = pair.template hfcRedCorrColl_as<o2::aod::HfcRedCorrColls>();
    if (collision.centrality() < centralityMin || collision.centrality() > centralityMax) {
      return;
    }
    double const ptTrig = trigCand.ptTrig();
    if constexpr (requires { trigCand.bdtScore0Trig(); }) { // ML selection on bkg score for Charm-Had case
      if (!isSelBdtScoreCut(trigCand, ptTrig)) {
        return;
      }
    }
    if (downSamplePairs < 1.) {
      float const pseudoRndm = ptTrig * 1000. - static_cast<int64_t>(ptTrig * 1000);
      if (ptTrig < ptMaxForDownSample && collision.centrality() < centMaxForDownSample && pseudoRndm >= downSamplePairs) {
        return;
      }
    }
    int const poolBin = getPoolBin<false>(collision, binPolicy);
    registry.fill(HIST("hPoolBinTrigSE"), poolBin);
    registry.fill(HIST("hPoolBinAssocSE"), poolBin);
    // if constexpr (FillTables) {
    //   if constexpr (requires { trigCand.bdtScore0Trig(); }) { // Separate Charm-Had and Had-Had cases
    //     rowPairSECharmHads(poolBin, ptTrig, pair.ptAssoc(), pair.deltaEta(), pair.deltaPhi(),
    //                        trigCand.invMassTrig(), trigCand.bdtScore0Trig(), trigCand.bdtScore1Trig(),
    //                        pair.nTpcCrossedRowsAssoc(), pair.itsClsMapAssoc(), pair.itsNClsAssoc(), pair.dcaXYAssoc(), pair.dcaZAssoc());
    //   } else {
    //     rowPairSEHadHads(poolBin, ptTrig, pair.ptAssoc(), pair.deltaEta(), pair.deltaPhi(),
    //                      trigCand.nTpcCrossedRowsTrig(), trigCand.itsClsMapTrig(), trigCand.itsNClsTrig(), trigCand.dcaXYTrig(), trigCand.dcaZTrig(),
    //                      pair.nTpcCrossedRowsAssoc(), pair.itsClsMapAssoc(), pair.itsNClsAssoc(), pair.dcaXYAssoc(), pair.dcaZAssoc());
    //   }
    //   rowCollInfos(collision.multiplicity(), collision.numPvContrib(), collision.centrality());
    // }
    if constexpr (FillSparses) {
      if constexpr (requires { trigCand.bdtScore0Trig(); }) { // Separate Charm-Had and Had-Had cases
        registry.fill(HIST("hSparseCorrelationsSECharmHad"), poolBin, ptTrig, pair.ptAssoc(), pair.deltaEta(),
                      pair.deltaPhi(), trigCand.invMassTrig()); // , trigCand.bdtScore0Trig(), trigCand.bdtScore1Trig());
      } else {
        registry.fill(HIST("hSparseCorrelationsSEHadHad"), poolBin, ptTrig, pair.ptAssoc(), pair.deltaEta(), pair.deltaPhi());
      }
    }
  }

  /// Save info for Mixed Event pairs
  /// \param collisions are the selected collisions
  /// \param pairs are the mixed event pairs of trigger candidates and associated tracks
  /// \param binPolicy is the binning policy for the correlation
  template <bool FillSparses, bool FillTables, typename TPairs, typename TBinningType>
  void fillMixedEvent(TPairs const& pairs,
                      TBinningType binPolicy)
  {
    for (const auto& [trigColl, trigCands, assocColl, assocTracks] : pairs) {
      if (trigCands.size() == 0 || assocTracks.size() == 0) {
        continue;
      }
      if (trigColl.centrality() < centralityMin || trigColl.centrality() > centralityMax ||
          assocColl.centrality() < centralityMin || assocColl.centrality() > centralityMax) {
        continue;
      }
      int const poolBinTrig = getPoolBin<true>(trigColl, binPolicy);
      int const poolBinAssoc = getPoolBin<true>(assocColl, binPolicy);
      if (poolBinAssoc != poolBinTrig) {
        LOGF(info, "Error, poolBins are different");
        continue;
      }
      registry.fill(HIST("hPoolBinTrigME"), poolBinTrig);
      registry.fill(HIST("hPoolBinAssocME"), poolBinAssoc);

      for (const auto& [trigCand, assocTrack] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(trigCands, assocTracks))) {
        // LOGF(info, "Mixed event tracks pair: (%d, %d) from events (%d, %d), track event: (%d, %d)", trigCand.index(), assocTrack.index(), trigColl.index(), assocColl.index(), trigCand.hfcRedFlowCollId(), assocTrack.hfcRedFlowCollId());
        double const deltaEta = getEta(assocTrack) - getEta(trigCand);
        if (std::abs(deltaEta) < deltaEtaAbsMin || std::abs(deltaEta) > deltaEtaAbsMax) {
          continue;
        }
        double const ptTrig = getPt(trigCand);
        if constexpr (requires { trigCand.bdtScore0Trig(); }) { // ML selection on bkg score for Charm-Had case
          if (!isSelBdtScoreCut(trigCand, ptTrig)) {
            continue;
          }
        }
        double const ptAssoc = getPt(assocTrack);
        if (downSamplePairs < 1.) {
          float const pseudoRndm = ptAssoc * 1000. - static_cast<int64_t>(ptAssoc * 1000);
          if (ptTrig < ptMaxForDownSample && trigColl.centrality() < centMaxForDownSample &&
              assocColl.centrality() < centMaxForDownSample && pseudoRndm >= downSamplePairs) {
            continue;
          }
        }
        double const deltaPhi = RecoDecay::constrainAngle(getPhi(assocTrack) - getPhi(trigCand), -o2::constants::math::PIHalf);
        // if constexpr (FillTables) {
        //   if constexpr (requires { trigCand.bdtScore0Trig(); }) { // Separate Charm-Had and Had-Had cases
        //     rowPairMECharmHads(poolBinTrig, ptTrig, ptAssoc, deltaEta, deltaPhi,
        //                        trigCand.invMassTrig(), trigCand.bdtScore0Trig(), trigCand.bdtScore1Trig(),
        //                        assocTrack.nTpcCrossedRowsAssoc(), assocTrack.itsClsMapAssoc(), assocTrack.itsNClsAssoc(), assocTrack.dcaXYAssoc(), assocTrack.dcaZAssoc());
        //   } else {
        //     rowPairMEHadHads(poolBinTrig, ptTrig, ptAssoc, deltaEta, deltaPhi,
        //                      trigCand.nTpcCrossedRowsAssoc(), trigCand.itsClsMapAssoc(), trigCand.itsNClsAssoc(), trigCand.dcaXYAssoc(), trigCand.dcaZAssoc(),
        //                      assocTrack.nTpcCrossedRowsAssoc(), assocTrack.itsClsMapAssoc(), assocTrack.itsNClsAssoc(), assocTrack.dcaXYAssoc(), assocTrack.dcaZAssoc());
        //   }
        //   rowCollInfos(trigColl.multiplicity(), trigColl.numPvContrib(), trigColl.centrality());
        // }
        if constexpr (FillSparses) {
          if constexpr (requires { trigCand.bdtScore0Trig(); }) { // Separate Charm-Had and Had-Had cases
            registry.fill(HIST("hSparseCorrelationsMECharmHad"), poolBinTrig, ptTrig, ptAssoc, deltaEta,
                          deltaPhi, trigCand.invMassTrig()); //, trigCand.bdtScore0Trig(), trigCand.bdtScore1Trig());
          } else {
            registry.fill(HIST("hSparseCorrelationsMEHadHad"), poolBinTrig, ptTrig, ptAssoc, deltaEta, deltaPhi);
          }
        }
      }
    }
  }

  /// Correlations for Same Event pairs
  /// \param poolBin collision pool bin based on multiplicity and z-vertex position
  /// \param trigCandsThisColl are the selected trigger candidates in the collision
  /// \param assocTracksThisColl are the selected associated tracks in the collision
  template <bool FillSparses, typename TTrigCand, typename TTrackAssoc>
  void doCorrelationsSameEvent(int poolBin,
                               const TTrigCand& trigCandsThisColl,
                               const TTrackAssoc& assocTracksThisColl)
  {
    for (const auto& trigCand : trigCandsThisColl) {
      double const ptTrig = trigCand.ptTrig();
      if constexpr (requires { trigCand.bdtScore0Trig(); }) { // ML selection on bkg score for Charm-Had case
        if (!isSelBdtScoreCut(trigCand, ptTrig)) {
          continue;
        }
      }

      for (const auto& assTrk : assocTracksThisColl) {
        // TODO: Remove Ds daughters
        /*if (assTrk.originTrackId() == candidate.prong0Id() ||
            assTrk.originTrackId() == candidate.prong1Id() ||
            assTrk.originTrackId() == candidate.prong2Id()) {
          continue;
        }*/
        // TODO: DCA cut
        double deltaPhi = RecoDecay::constrainAngle(assTrk.phiAssoc() - trigCand.phiTrig(), -o2::constants::math::PIHalf);
        double deltaEta = assTrk.etaAssoc() - trigCand.etaTrig();
        if constexpr (FillSparses) {
          if constexpr (requires { trigCand.bdtScore0Trig(); }) { // Separate Charm-Had and Had-Had cases
            registry.fill(HIST("hSparseCorrelationsSECharmHad"), poolBin, ptTrig, assTrk.ptAssoc(), deltaEta,
                          deltaPhi, trigCand.invMassTrig());
          } else {
            registry.fill(HIST("hSparseCorrelationsSEHadHad"), poolBin, ptTrig, assTrk.ptAssoc(), deltaEta, deltaPhi);
          }
        }
      }
    }
  }

  /// Correlations for Mixed Event pairs
  /// \param collisions are the selected collisions
  /// \param trigCands are the trigger candidates
  /// \param assocTracks  are the associated tracks
  /// \param binPolicy is the binning policy for the correlation
  template <bool FillSparses, typename TCollision, typename TTrigCand, typename TTrackAssoc, typename TBinningType>
  void doCorrelationsMixedEvent(const TCollision& collisions,
                                const TTrigCand& trigCands,
                                const TTrackAssoc& assocTracks,
                                TBinningType binPolicy)
  {
    auto tracksTuple = std::make_tuple(trigCands, assocTracks);

    Pair<TCollision, TTrigCand, TTrackAssoc, TBinningType> pairData{binPolicy, numberEventsMixed, -1, collisions, tracksTuple, &cache};

    for (const auto& [c1, tracks1, c2, tracks2] : pairData) {
      if (tracks1.size() == 0) {
        continue;
      }

      int poolBin = binPolicy.getBin({c2.posZ(), c2.multiplicity()});
      int poolBinTrigCand = binPolicy.getBin({c1.posZ(), c1.multiplicity()});

      if (poolBin != poolBinTrigCand) {
        LOGF(info, "Error, poolBins are different");
        continue;
      }

      for (const auto& [trigCand, assTrk] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        if (!isSelBdtScoreCut(trigCand, trigCand.ptTrig())) {
          continue;
        }
        LOGF(info, "Mixed event tracks pair: (%d, %d) from events (%d, %d), track event: (%d, %d)", trigCand.index(), assTrk.index(), c1.index(), c2.index(), trigCand.hfcRedCorrCollId(), assTrk.hfcRedCorrCollId());

        double deltaPhi = RecoDecay::constrainAngle(assTrk.phiAssoc() - trigCand.phiTrig(), -o2::constants::math::PIHalf);
        double deltaEta = assTrk.etaAssoc() - trigCand.etaTrig();
        if constexpr (FillSparses) {
          if constexpr (requires { trigCand.bdtScore0Trig(); }) { // Separate Charm-Had and Had-Had cases
            registry.fill(HIST("hSparseCorrelationsMECharmHad"), poolBin, trigCand.ptTrig(), assTrk.ptAssoc(), deltaEta,
                          deltaPhi, trigCand.invMassTrig());
          } else {
            registry.fill(HIST("hSparseCorrelationsMEHadHad"), poolBin, trigCand.ptTrig(), assTrk.ptAssoc(), deltaEta, deltaPhi);
          }
        }
      }
    }
  }

  void processSameEventCharmHadWMultMix(SameEvtPairsChHad::iterator const& pair,
                                        aod::HfcRedTrigCharms const&,
                                        aod::HfcRedCorrColls const&)
  {
    BinningMultPosZ binPolicyPosZMult{{zPoolBins, multPoolBins}, true};
    auto trigCand = pair.template hfcRedTrigCharm_as<aod::HfcRedTrigCharms>();
    if (fillSparses && fillTables) {
      fillSameEvent<true, true>(pair, trigCand, binPolicyPosZMult);
    } else if (fillSparses) {
      fillSameEvent<true, false>(pair, trigCand, binPolicyPosZMult);
    } else if (fillTables) {
      fillSameEvent<false, true>(pair, trigCand, binPolicyPosZMult);
    }
  }
  PROCESS_SWITCH(HfCorrelatorFlowCharmHadronsReduced, processSameEventCharmHadWMultMix, "Process Same Event for Charm-Had with multiplicity pools", true);

  void processSameEventHadHadWMultMix(SameEvtPairsHadHad::iterator const& pair,
                                      aod::HfcRedTrigTracks const&,
                                      aod::HfcRedCorrColls const&)
  {
    BinningMultPosZ binPolicyPosZMult{{zPoolBins, multPoolBins}, true};
    auto trigCand = pair.template hfcRedTrigTrack_as<aod::HfcRedTrigTracks>();
    if (fillSparses && fillTables) {
      fillSameEvent<true, true>(pair, trigCand, binPolicyPosZMult);
    } else if (fillSparses) {
      fillSameEvent<true, false>(pair, trigCand, binPolicyPosZMult);
    } else if (fillTables) {
      fillSameEvent<false, true>(pair, trigCand, binPolicyPosZMult);
    }
  }
  PROCESS_SWITCH(HfCorrelatorFlowCharmHadronsReduced, processSameEventHadHadWMultMix, "Process Same Event for Had-Had with multiplicity pools", false);

  void processSameEventCharmHadWCentMix(SameEvtPairsChHad::iterator const& pair,
                                        aod::HfcRedTrigCharms const&,
                                        aod::HfcRedCorrColls const&)
  {
    BinningCentPosZ binPolicyPosZCent{{zPoolBins, centPoolBins}, true};
    auto trigCand = pair.template hfcRedTrigCharm_as<aod::HfcRedTrigCharms>();
    if (fillSparses && fillTables) {
      fillSameEvent<true, true>(pair, trigCand, binPolicyPosZCent);
    } else if (fillSparses) {
      fillSameEvent<true, false>(pair, trigCand, binPolicyPosZCent);
    } else if (fillTables) {
      fillSameEvent<false, true>(pair, trigCand, binPolicyPosZCent);
    }
  }
  PROCESS_SWITCH(HfCorrelatorFlowCharmHadronsReduced, processSameEventCharmHadWCentMix, "Process Same Event for Charm-Had with centrality pools", true);

  void processSameEventHadHadWCentMix(SameEvtPairsHadHad::iterator const& pair,
                                      aod::HfcRedTrigTracks const&,
                                      aod::HfcRedCorrColls const&)
  {
    BinningCentPosZ binPolicyPosZCent{{zPoolBins, centPoolBins}, true};
    auto trigCand = pair.template hfcRedTrigTrack_as<aod::HfcRedTrigTracks>();
    if (fillSparses && fillTables) {
      fillSameEvent<true, true>(pair, trigCand, binPolicyPosZCent);
    } else if (fillSparses) {
      fillSameEvent<true, false>(pair, trigCand, binPolicyPosZCent);
    } else if (fillTables) {
      fillSameEvent<false, true>(pair, trigCand, binPolicyPosZCent);
    }
  }
  PROCESS_SWITCH(HfCorrelatorFlowCharmHadronsReduced, processSameEventHadHadWCentMix, "Process Same Event for Had-Had with centrality pools", false);

  void processMixedEventCharmHadWCentMix(aod::HfcRedCorrColls const& collisions,
                                         TrigCharmCands const& candidates,
                                         AssocTracks const& tracks)
  {
    BinningCentPosZ binPolicyPosZCent{{zPoolBins, centPoolBins}, true};
    auto pairsTuple = std::make_tuple(candidates, tracks);
    Pair<aod::HfcRedCorrColls, TrigCharmCands, AssocTracks, BinningCentPosZ> const pairs{binPolicyPosZCent, numberEventsMixed, -1, collisions, pairsTuple, &cache};
    if (fillSparses && fillTables) {
      fillMixedEvent<true, true>(pairs, binPolicyPosZCent);
    } else if (fillSparses) {
      fillMixedEvent<true, false>(pairs, binPolicyPosZCent);
    } else if (fillTables) {
      fillMixedEvent<false, true>(pairs, binPolicyPosZCent);
    }
  }
  PROCESS_SWITCH(HfCorrelatorFlowCharmHadronsReduced, processMixedEventCharmHadWCentMix, "Process Mixed Event for Charm-Had with centrality pools", false);

  void processMixedEventCharmHadWMultMix(aod::HfcRedCorrColls const& collisions,
                                         TrigCharmCands const& candidates,
                                         AssocTracks const& tracks)
  {
    BinningMultPosZ binPolicyPosZMult{{zPoolBins, multPoolBins}, true};
    auto pairsTuple = std::make_tuple(candidates, tracks);
    Pair<aod::HfcRedCorrColls, TrigCharmCands, AssocTracks, BinningMultPosZ> const pairs{binPolicyPosZMult, numberEventsMixed, -1, collisions, pairsTuple, &cache};
    if (fillSparses && fillTables) {
      fillMixedEvent<true, true>(pairs, binPolicyPosZMult);
    } else if (fillSparses) {
      fillMixedEvent<true, false>(pairs, binPolicyPosZMult);
    } else if (fillTables) {
      fillMixedEvent<false, true>(pairs, binPolicyPosZMult);
    }
  }
  PROCESS_SWITCH(HfCorrelatorFlowCharmHadronsReduced, processMixedEventCharmHadWMultMix, "Process Mixed Event for Charm-Had with multiplicity pools", false);

  void processMixedEventHadHadWCentMix(aod::HfcRedCorrColls const& collisions,
                                       AssocTracks const& tracks)
  {
    BinningCentPosZ binPolicyPosZCent{{zPoolBins, centPoolBins}, true};
    auto tracksTuple = std::make_tuple(tracks);
    SameKindPair<aod::HfcRedCorrColls, AssocTracks, BinningCentPosZ> const pairs{binPolicyPosZCent, numberEventsMixed, -1, collisions, tracksTuple, &cache};
    if (fillSparses && fillTables) {
      fillMixedEvent<true, true>(pairs, binPolicyPosZCent);
    } else if (fillSparses) {
      fillMixedEvent<true, false>(pairs, binPolicyPosZCent);
    } else if (fillTables) {
      fillMixedEvent<false, true>(pairs, binPolicyPosZCent);
    }
  }
  PROCESS_SWITCH(HfCorrelatorFlowCharmHadronsReduced, processMixedEventHadHadWCentMix, "Process Mixed Event for Had-Had with centrality pools", false);

  void processMixedEventHadHadWMultMix(aod::HfcRedCorrColls const& collisions,
                                       AssocTracks const& tracks)
  {
    BinningMultPosZ binPolicyPosZMult{{zPoolBins, multPoolBins}, true};
    auto tracksTuple = std::make_tuple(tracks);
    SameKindPair<aod::HfcRedCorrColls, AssocTracks, BinningMultPosZ> const pairs{binPolicyPosZMult, numberEventsMixed, -1, collisions, tracksTuple, &cache};
    if (fillSparses && fillTables) {
      fillMixedEvent<true, true>(pairs, binPolicyPosZMult);
    } else if (fillSparses) {
      fillMixedEvent<true, false>(pairs, binPolicyPosZMult);
    } else if (fillTables) {
      fillMixedEvent<false, true>(pairs, binPolicyPosZMult);
    }
  }
  PROCESS_SWITCH(HfCorrelatorFlowCharmHadronsReduced, processMixedEventHadHadWMultMix, "Process Mixed Event for Had-Had with multiplicity pools", false);

  void processCharmTriggers(aod::HfcRedTrigCharms const& trigCands,
                            aod::HfcRedCorrColls const&)
  {
    for (const auto& trigCand : trigCands) {
      auto collision = trigCand.template hfcRedCorrColl_as<o2::aod::HfcRedCorrColls>();
      if (collision.centrality() < centralityMin || collision.centrality() > centralityMax) {
        continue;
      }
      if (!isSelBdtScoreCut(trigCand, trigCand.ptTrig())) {
        continue;
      }
      registry.fill(HIST("hSparseTrigCandsCharm"), trigCand.invMassTrig(), trigCand.ptTrig(), trigCand.bdtScore0Trig(), trigCand.bdtScore1Trig());
    }
  }
  PROCESS_SWITCH(HfCorrelatorFlowCharmHadronsReduced, processCharmTriggers, "Process charm trigger info", false);

  void processSameEventCharmHadWCentMixBase(aod::HfcRedCorrColls const& collisions,
                                            TrigCharmCands const& candidates,
                                            aod::HfcRedAssBases const& tracks)
  {
    BinningCentPosZ binPolicyPosZCent{{zPoolBins, centPoolBins}, true};

    for (const auto& collision : collisions) {
      if (collision.centrality() < centralityMin || collision.centrality() > centralityMax) {
        continue;
      }
      int poolBin = binPolicyPosZCent.getBin({collision.posZ(), collision.multiplicity()});

      auto thisCollId = collision.globalIndex();
      auto candsThisColl = candidates.sliceBy(trigCharmCandsPerCol, thisCollId);
      auto tracksThisColl = tracks.sliceBy(assocTracksPerCol, thisCollId);

      doCorrelationsSameEvent<true>(poolBin, candsThisColl, tracksThisColl);
    }
  }
  PROCESS_SWITCH(HfCorrelatorFlowCharmHadronsReduced, processSameEventCharmHadWCentMixBase, "Process Same Event base", false);

  void processMixedEventCharmHadWCentMixBase(aod::HfcRedCorrColls const& collisions,
                                             TrigCharmCands const& candidates,
                                             aod::HfcRedAssBases const& tracks)
  {
    BinningCentPosZ binPolicyPosZCent{{zPoolBins, centPoolBins}, true};

    doCorrelationsMixedEvent<true>(collisions, candidates, tracks, binPolicyPosZCent);
  }
  PROCESS_SWITCH(HfCorrelatorFlowCharmHadronsReduced, processMixedEventCharmHadWCentMixBase, "Process Mixed Event base", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCorrelatorFlowCharmHadronsReduced>(cfgc)};
}
