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

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/HFC/DataModel/DerivedDataCorrelationTables.h"

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
#include <Framework/SliceCache.h>
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

struct HfCorrelatorFlowCharmHadronsReduced {
  Produces<aod::HfcRedChHads> entryCharmHadPair;
  Produces<aod::HfcRedHadHads> entryHadHadPair;

  Configurable<bool> fillSparses{"fillSparses", true, "Fill sparse histograms"};
  Configurable<bool> fillTables{"fillTables", false, "Fill tables"};
  Configurable<int> numberEventsMixed{"numberEventsMixed", 5, "Number of events mixed in ME process"};
  Configurable<std::vector<double>> binsPtTrig{"binsPtTrig", std::vector<double>{1., 3., 5., 8., 16., 36.}, "pT bin limits for trigger candidates"};
  Configurable<std::vector<double>> binsPtAssoc{"binsPtAssoc", std::vector<double>{0.3, 1., 2., 50.}, "pT bin limits for associated particles"};
  Configurable<float> deltaEtaAbsMin{"deltaEtaAbsMin", 0.5, "min. pair delta eta"};
  Configurable<float> deltaEtaAbsMax{"deltaEtaAbsMax", 2., "max. pair delta eta"};
  Configurable<float> dcaXYTrackMax{"dcaXYTrackMax", 1., "max. track DCA XY"};
  Configurable<float> dcaZTrackMax{"dcaZTrackMax", 1., "max. track DCA Z"};
  Configurable<int> tpcCrossedRowsMin{"tpcCrossedRowsMin", 1, "min. TPC crossed rows"};
  Configurable<int> itsNClsMin{"itsNClsMin", 1, "min. ITS clusters"};

  SliceCache cache;

  using AssocTracks = soa::Filtered<soa::Join<aod::HfcRedTrkAssocs, aod::HfcRedTrkSels>>;

  Filter filterSelectTrackData = (nabs(aod::hf_assoc_track_reduced::dcaXY) < dcaXYTrackMax) && (nabs(aod::hf_assoc_track_reduced::dcaZ) < dcaZTrackMax) && (aod::hf_assoc_track_reduced::nTpcCrossedRows > tpcCrossedRowsMin) && (aod::hf_assoc_track_reduced::itsNCls > itsNClsMin);

  Preslice<AssocTracks> tracksPerCol = aod::hf_candidate_reduced::hfcRedFlowCollId;
  Preslice<aod::HfcRedCharmTrigs> candsPerCol = aod::hf_candidate_reduced::hfcRedFlowCollId;

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

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    std::array<bool, 8> doprocess{doprocessSameEventCharmHadWCentMix, doprocessSameEventCharmHadWMultMix, doprocessMixedEventCharmHadWCentMix, doprocessMixedEventCharmHadWMultMix,
                                  doprocessSameEventHadHadWCentMix, doprocessSameEventHadHadWMultMix, doprocessMixedEventHadHadWCentMix, doprocessMixedEventHadHadWMultMix};
    if ((std::accumulate(doprocess.begin(), doprocess.end(), 0)) > 1) {
      LOGP(fatal, "Only one process function should be enabled! Please check your configuration!");
      if (!((doprocessSameEventCharmHadWCentMix && doprocessMixedEventCharmHadWCentMix) || (doprocessSameEventCharmHadWMultMix && doprocessMixedEventCharmHadWMultMix))) {
        LOG(fatal) << "Different binning policies between Same Event and Mixed Event";
      }
      if (!((doprocessSameEventHadHadWCentMix && doprocessMixedEventHadHadWCentMix) || (doprocessSameEventHadHadWMultMix && doprocessMixedEventHadHadWMultMix))) {
        LOG(fatal) << "Different binning policies between Same Event and Mixed Event";
      }
    }

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
    if (doprocessSameEventCharmHadWCentMix || doprocessMixedEventCharmHadWCentMix || doprocessSameEventHadHadWCentMix || doprocessMixedEventHadHadWCentMix) {
      registry.add("hCent", "Centrality", {HistType::kTH2F, {{axisCent}, {axisPoolBin}}});
    } else {
      registry.add("hMultFT0M", "Multiplicity FT0M", {HistType::kTH2F, {{axisMultFT0M}, {axisPoolBin}}});
    }
    registry.add("hZVtx", "z vertex", {HistType::kTH2F, {{axisPosZ}, {axisPoolBin}}});
    registry.add("hCollisionPoolBin", "Collision pool bin", {HistType::kTH1F, {axisPoolBin}});
    registry.add("hPoolBinTrig", "Trigger candidates pool bin", {HistType::kTH1F, {axisPoolBin}});
    registry.add("hPhiVsPtTrig", "Trigger candidates phiVsPt", {HistType::kTH2F, {{axisPhi}, {axisPtTrig}}});
    registry.add("hEtaVsPtTrig", "Trigger candidates etaVsPt", {HistType::kTH2F, {{axisEta}, {axisPtTrig}}});
    registry.add("hPoolBinAssoc", "Associated particles pool bin", {HistType::kTH1F, {axisPoolBin}});
    registry.add("hPhiVsPtAssoc", "Associated particles phiVsPt", {HistType::kTH3F, {{axisPhi}, {axisPtTrig}, {axisPtAssoc}}});
    registry.add("hEtaVsPtAssoc", "Associated particles etaVsPt", {HistType::kTH3F, {{axisEta}, {axisPtTrig}, {axisPtAssoc}}});

    if (fillSparses) {
      std::vector<AxisSpec> axes = {axisPtTrig, axisPtAssoc, axisDeltaEta, axisDeltaPhi, axisPoolBin};
      if (doprocessSameEventHadHadWCentMix || doprocessSameEventHadHadWMultMix) {
        registry.add("hSparseCorrelationsSEHadHad", "THn for SE Had-Had correlations", HistType::kTHnSparseF, axes);
      } else if (doprocessMixedEventHadHadWCentMix || doprocessMixedEventHadHadWMultMix) {
        registry.add("hSparseCorrelationsMEHadHad", "THn for ME Had-Had correlations", HistType::kTHnSparseF, axes);
      } else {
        axes.insert(axes.end(), {axisMlOne, axisMlTwo, axisInvMass});
        if (doprocessSameEventCharmHadWCentMix || doprocessSameEventCharmHadWMultMix) {
          registry.add("hSparseCorrelationsSECharmHad", "THn for SE Charm-Had correlations", HistType::kTHnSparseF, axes);
        } else if (doprocessMixedEventCharmHadWCentMix || doprocessMixedEventCharmHadWMultMix) {
          registry.add("hSparseCorrelationsMECharmHad", "THn for ME Charm-Had correlations", HistType::kTHnSparseF, axes);
        }
      }
    }
  }

  /// Get charm candidate or hadron track pT
  /// \param track is the candidate
  template <typename TTrack>
  double getPt(const TTrack& track)
  {
    if constexpr (requires { track.ptAssocTrack(); }) {
      return track.ptAssocTrack();
    } else {
      return track.ptCand();
    }
  }

  /// Get charm candidate or hadron track eta
  /// \param track is the candidate
  template <typename TTrack>
  double getEta(const TTrack& track)
  {
    if constexpr (requires { track.etaAssocTrack(); }) {
      return track.etaAssocTrack();
    } else {
      return track.etaCand();
    }
  }

  /// Get charm candidate or hadron track phi
  /// \param track is the candidate
  template <typename TTrack>
  double getPhi(const TTrack& track)
  {
    if constexpr (requires { track.phiAssocTrack(); }) {
      return track.phiAssocTrack();
    } else {
      return track.phiCand();
    }
  }

  /// Get the binning pool associated to the collision
  /// \param collision is the collision
  /// \param corrBinning is the binning policy for the correlation
  template <bool fillHistos, typename TColl, typename TBinningType>
  int getPoolBin(const TColl& collision, const TBinningType& corrBinning)
  {
    int poolBin{0};
    if constexpr (std::is_same_v<TBinningType, BinningTypeDerivedCent>) {
      poolBin = corrBinning.getBin(std::make_tuple(collision.posZ(), collision.centrality()));
      if constexpr (fillHistos) {
        registry.fill(HIST("hCent"), collision.centrality(), poolBin);
      }
    } else if constexpr (std::is_same_v<TBinningType, BinningTypeDerivedMult>) {
      poolBin = corrBinning.getBin(std::make_tuple(collision.posZ(), collision.multiplicity()));
      if constexpr (fillHistos) {
        registry.fill(HIST("hMultFT0M"), collision.multiplicity(), poolBin);
      }
    }
    return poolBin;
  }

  /// Reject daughter-track pairs and same-track pairs
  /// \param cand is the trigger candidate
  /// \param track is the associated track
  template <typename TTrigPart, typename TTrack>
  bool rejSameEvtPair(const TTrigPart& cand, const TTrack& track)
  {
    if constexpr (requires { cand.originTrackId(); }) {
      // Remove same track pairs for Had-Had correlations
      return (cand.originTrackId() == track.originTrackId());
    } else {
      // Remove pairs with 2- and 3-prong daughters (prong2Id returns -1 for 2-prongs)
      return ((cand.prong0Id() == track.originTrackId()) || (cand.prong1Id() == track.originTrackId()) || (cand.prong2Id() == track.originTrackId()));
    }
  }

  /// Slice trigger candidates by collision
  /// \param cands are the trigger candidates
  /// \param collId is the collision index
  template <typename TTrigCands>
  auto sliceTrigCands(TTrigCands const& cands, const int collId)
  {
    if constexpr (std::is_same_v<TTrigCands, AssocTracks>) {
      return cands.sliceBy(tracksPerCol, collId);
    } else {
      return cands.sliceBy(candsPerCol, collId);
    }
  }

  /// Fill Charm-Hadron correlation table and sparse
  /// \param trigCand is the trigger charm hadron candidate
  /// \param assocTrack is the associated hadron track
  /// \param poolBin is the pool bin of the collision
  template <bool isMixedEvent, typename TTrigCand, typename TTrack>
  void fillCharmHadInfo(TTrigCand const& trigCand,
                        TTrack const& assocTrack,
                        const int poolBin)
  {
    double deltaEta = getEta(assocTrack) - getEta(trigCand);
    double deltaPhi = RecoDecay::constrainAngle(getPhi(assocTrack) - getPhi(trigCand), -o2::constants::math::PIHalf);
    if (fillTables) {
      entryCharmHadPair(trigCand.globalIndex(), assocTrack.globalIndex(), deltaEta, deltaPhi, poolBin);
    }
    if (fillSparses) {
      if constexpr (isMixedEvent) {
        registry.fill(HIST("hSparseCorrelationsMECharmHad"), getPt(trigCand), getPt(assocTrack),
                      deltaEta, deltaPhi, poolBin, trigCand.bdtScore0(),
                      trigCand.bdtScore1(), trigCand.invMassCand());
      } else {
        registry.fill(HIST("hSparseCorrelationsSECharmHad"), getPt(trigCand), getPt(assocTrack),
                      deltaEta, deltaPhi, poolBin, trigCand.bdtScore0(),
                      trigCand.bdtScore1(), trigCand.invMassCand());
      }
    }
  }

  /// Fill Hadron-Hadron correlation table and sparse
  /// \param trigCand is the trigger hadron candidate
  /// \param assocTrack is the associated hadron track
  /// \param poolBin is the pool bin of the collision
  template <bool isMixedEvent, typename TCand>
  void fillHadHadInfo(TCand const& trigCand,
                      TCand const& assocTrack,
                      const int poolBin)
  {
    double deltaEta = getEta(assocTrack) - getEta(trigCand);
    double deltaPhi = RecoDecay::constrainAngle(getPhi(assocTrack) - getPhi(trigCand), -o2::constants::math::PIHalf);
    if (fillTables) {
      entryHadHadPair(trigCand.globalIndex(), assocTrack.globalIndex(), deltaEta, deltaPhi, poolBin);
    }
    if (fillSparses) {
      if constexpr (isMixedEvent) {
        registry.fill(HIST("hSparseCorrelationsMEHadHad"), getPt(trigCand), getPt(assocTrack), deltaEta, deltaPhi, poolBin);
      } else {
        registry.fill(HIST("hSparseCorrelationsSEHadHad"), getPt(trigCand), getPt(assocTrack), deltaEta, deltaPhi, poolBin);
      }
    }
  }

  /// Save info for Same Event pairs
  /// \param collisions are the selected collisions
  /// \param trigCands are the selected trigger candidates
  /// \param assocTracks are the selected associated tracks
  /// \param corrBinning is the binning policy for the correlation
  template <typename TTrigCands, typename TAssocTracks, typename TBinningType>
  void fillSameEvent(aod::HfcRedFlowColls const& collisions,
                     TTrigCands const& trigCands,
                     TAssocTracks const& assocTracks,
                     TBinningType corrBinning)
  {
    for (const auto& collision : collisions) {
      int poolBin = getPoolBin<true>(collision, corrBinning);
      registry.fill(HIST("hCollisionPoolBin"), poolBin);
      registry.fill(HIST("hZVtx"), collision.posZ(), poolBin);

      auto thisCollId = collision.globalIndex();
      auto trigCandsThisColl = sliceTrigCands(trigCands, thisCollId);
      auto assocTracksThisColl = assocTracks.sliceBy(tracksPerCol, thisCollId);

      for (const auto& trigCand : trigCandsThisColl) {
        registry.fill(HIST("hPoolBinTrig"), poolBin);
        registry.fill(HIST("hPhiVsPtTrig"), RecoDecay::constrainAngle(getPhi(trigCand), -o2::constants::math::PIHalf), getPt(trigCand));
        registry.fill(HIST("hEtaVsPtTrig"), getEta(trigCand), getPt(trigCand));
        for (const auto& assocTrack : assocTracksThisColl) {
          if (rejSameEvtPair(trigCand, assocTrack)) {
            continue;
          }
          double deltaEta = getEta(assocTrack) - getEta(trigCand);
          if (std::abs(deltaEta) < deltaEtaAbsMin || std::abs(deltaEta) > deltaEtaAbsMax) {
            continue;
          }
          registry.fill(HIST("hPoolBinAssoc"), poolBin);
          registry.fill(HIST("hPhiVsPtAssoc"), RecoDecay::constrainAngle(getPhi(assocTrack), -o2::constants::math::PIHalf), getPt(trigCand), getPt(assocTrack));
          registry.fill(HIST("hEtaVsPtAssoc"), getEta(assocTrack), getPt(trigCand), getPt(assocTrack));

          if constexpr (std::is_same_v<TTrigCands, AssocTracks>) {
            fillHadHadInfo<false>(trigCand, assocTrack, poolBin);
          } else {
            fillCharmHadInfo<false>(trigCand, assocTrack, poolBin);
          }
        }
      }
    }
  }

  /// Save info for Mixed Event pairs
  /// \param collisions are the selected collisions
  /// \param trigCands are the selected trigger candidates
  /// \param assocTracks are the selected associated tracks
  /// \param corrBinning is the binning policy for the correlation
  template <typename TTrigCands, typename TAssocTracks, typename TBinningType>
  void fillMixedEvent(aod::HfcRedFlowColls const& collisions,
                      TTrigCands const& trigCands,
                      TAssocTracks const& assocTracks,
                      TBinningType corrBinning)
  {
    for (const auto& collision : collisions) {
      int poolBin = getPoolBin<true>(collision, corrBinning);
      registry.fill(HIST("hCollisionPoolBin"), poolBin);
      registry.fill(HIST("hZVtx"), collision.posZ(), poolBin);

      auto thisCollId = collision.globalIndex();
      auto trigCandsThisColl = sliceTrigCands(trigCands, thisCollId);
      auto assocTracksThisColl = assocTracks.sliceBy(tracksPerCol, thisCollId);
      for (const auto& trigCand : trigCandsThisColl) {
        registry.fill(HIST("hPoolBinTrig"), poolBin);
        registry.fill(HIST("hPhiVsPtTrig"), RecoDecay::constrainAngle(getPhi(trigCand), -o2::constants::math::PIHalf), getPt(trigCand));
        registry.fill(HIST("hEtaVsPtTrig"), getEta(trigCand), getPt(trigCand));
        for (const auto& assocTrack : assocTracksThisColl) {
          registry.fill(HIST("hPoolBinAssoc"), poolBin);
          registry.fill(HIST("hPhiVsPtAssoc"), RecoDecay::constrainAngle(getPhi(assocTrack), -o2::constants::math::PIHalf), getPt(trigCand), getPt(assocTrack));
          registry.fill(HIST("hEtaVsPtAssoc"), getEta(assocTrack), getPt(trigCand), getPt(assocTrack));
        }
      }
    }

    auto pairsTuple = std::make_tuple(trigCands, assocTracks);
    Pair<aod::HfcRedFlowColls, TTrigCands, TAssocTracks, TBinningType> pairData{corrBinning, numberEventsMixed, -1, collisions, pairsTuple, &cache};

    for (const auto& [trigColl, trigCands, assocColl, assocTracks] : pairData) {
      if (trigCands.size() == 0 || assocTracks.size() == 0) {
        continue;
      }
      int poolBinCharm = getPoolBin<false>(trigColl, corrBinning);
      int poolBinAssoc = getPoolBin<false>(assocColl, corrBinning);
      if (poolBinAssoc != poolBinCharm) {
        LOGF(info, "Error, poolBins are different");
        continue;
      }

      for (const auto& [trigCand, assocTrack] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(trigCands, assocTracks))) {
        double deltaEta = getEta(assocTrack) - getEta(trigCand);
        if (std::abs(deltaEta) < deltaEtaAbsMin || std::abs(deltaEta) > deltaEtaAbsMax) {
          continue;
        }
        // LOGF(info, "Mixed event tracks pair: (%d, %d) from events (%d, %d), track event: (%d, %d)", trigCand.index(), assocTrack.index(), trigColl.index(), assocColl.index(), trigCand.hfcRedFlowCollId(), assocTrack.hfcRedFlowCollId());
        if constexpr (std::is_same_v<TTrigCands, AssocTracks>) {
          fillHadHadInfo<true>(trigCand, assocTrack, poolBinCharm);
        } else {
          fillCharmHadInfo<true>(trigCand, assocTrack, poolBinCharm);
        }
      }
    }
  }

  void processSameEventCharmHadWCentMix(aod::HfcRedFlowColls const& collisions,
                                        soa::Join<aod::HfcRedCharmTrigs, aod::HfcRedCharmMls> const& candidates,
                                        AssocTracks const& tracks)
  {
    ColumnBinningPolicy<aod::hf_collisions_reduced::PosZ, aod::hf_collisions_reduced::Centrality> corrBinningCent{{zPoolBins, centPoolBins}, true};
    fillSameEvent(collisions, candidates, tracks, corrBinningCent);
  }
  PROCESS_SWITCH(HfCorrelatorFlowCharmHadronsReduced, processSameEventCharmHadWCentMix, "Process Same Event for Charm-Had with centrality pools", true);

  void processSameEventCharmHadWMultMix(aod::HfcRedFlowColls const& collisions,
                                        soa::Join<aod::HfcRedCharmTrigs, aod::HfcRedCharmMls> const& candidates,
                                        AssocTracks const& tracks)
  {
    ColumnBinningPolicy<aod::hf_collisions_reduced::PosZ, aod::hf_collisions_reduced::Multiplicity> corrBinningMult{{zPoolBins, multPoolBins}, true};
    fillSameEvent(collisions, candidates, tracks, corrBinningMult);
  }
  PROCESS_SWITCH(HfCorrelatorFlowCharmHadronsReduced, processSameEventCharmHadWMultMix, "Process Same Event for Charm-Had with multiplicity pools", false);

  void processMixedEventCharmHadWCentMix(aod::HfcRedFlowColls const& collisions,
                                         soa::Join<aod::HfcRedCharmTrigs, aod::HfcRedCharmMls> const& candidates,
                                         AssocTracks const& tracks)
  {
    ColumnBinningPolicy<aod::hf_collisions_reduced::PosZ, aod::hf_collisions_reduced::Centrality> corrBinningCent{{zPoolBins, centPoolBins}, true};
    fillMixedEvent(collisions, candidates, tracks, corrBinningCent);
  }
  PROCESS_SWITCH(HfCorrelatorFlowCharmHadronsReduced, processMixedEventCharmHadWCentMix, "Process Mixed Event for Charm-Had with centrality pools", false);

  void processMixedEventCharmHadWMultMix(aod::HfcRedFlowColls const& collisions,
                                         soa::Join<aod::HfcRedCharmTrigs, aod::HfcRedCharmMls> const& candidates,
                                         AssocTracks const& tracks)
  {
    ColumnBinningPolicy<aod::hf_collisions_reduced::PosZ, aod::hf_collisions_reduced::Multiplicity> corrBinningMult{{zPoolBins, multPoolBins}, true};
    fillMixedEvent(collisions, candidates, tracks, corrBinningMult);
  }
  PROCESS_SWITCH(HfCorrelatorFlowCharmHadronsReduced, processMixedEventCharmHadWMultMix, "Process Mixed Event for Charm-Had with multiplicity pools", false);

  void processSameEventHadHadWCentMix(aod::HfcRedFlowColls const& collisions,
                                      AssocTracks const& tracks)
  {
    ColumnBinningPolicy<aod::hf_collisions_reduced::PosZ, aod::hf_collisions_reduced::Centrality> corrBinningCent{{zPoolBins, centPoolBins}, true};
    fillSameEvent(collisions, tracks, tracks, corrBinningCent);
  }
  PROCESS_SWITCH(HfCorrelatorFlowCharmHadronsReduced, processSameEventHadHadWCentMix, "Process Same Event for Had-Had with centrality pools", false);

  void processSameEventHadHadWMultMix(aod::HfcRedFlowColls const& collisions,
                                      AssocTracks const& tracks)
  {
    ColumnBinningPolicy<aod::hf_collisions_reduced::PosZ, aod::hf_collisions_reduced::Multiplicity> corrBinningMult{{zPoolBins, multPoolBins}, true};
    fillSameEvent(collisions, tracks, tracks, corrBinningMult);
  }
  PROCESS_SWITCH(HfCorrelatorFlowCharmHadronsReduced, processSameEventHadHadWMultMix, "Process Same Event for Had-Had with multiplicity pools", false);

  void processMixedEventHadHadWCentMix(aod::HfcRedFlowColls const& collisions,
                                       AssocTracks const& tracks)
  {
    ColumnBinningPolicy<aod::hf_collisions_reduced::PosZ, aod::hf_collisions_reduced::Centrality> corrBinningCent{{zPoolBins, centPoolBins}, true};
    fillMixedEvent(collisions, tracks, tracks, corrBinningCent);
  }
  PROCESS_SWITCH(HfCorrelatorFlowCharmHadronsReduced, processMixedEventHadHadWCentMix, "Process Mixed Event for Had-Had with centrality pools", false);

  void processMixedEventHadHadWMultMix(aod::HfcRedFlowColls const& collisions,
                                       AssocTracks const& tracks)
  {
    ColumnBinningPolicy<aod::hf_collisions_reduced::PosZ, aod::hf_collisions_reduced::Multiplicity> corrBinningMult{{zPoolBins, multPoolBins}, true};
    fillMixedEvent(collisions, tracks, tracks, corrBinningMult);
  }
  PROCESS_SWITCH(HfCorrelatorFlowCharmHadronsReduced, processMixedEventHadHadWMultMix, "Process Mixed Event for Had-Had with multiplicity pools", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCorrelatorFlowCharmHadronsReduced>(cfgc)};
}
