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

/// \file correlatorDsHadronsReduced.cxx
/// \brief Ds-Hadrons correlator task for offline analysis
/// \author Samuele Cattaruzzi <samuele.cattaruzzi@cern.ch>

#include <vector>

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "PWGHF/HFC/DataModel/CorrelationTables.h"
#include "PWGHF/HFC/DataModel/DerivedDataCorrelationTables.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::constants::physics;
using namespace o2::constants::math;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// Returns deltaPhi value in range [-pi/2., 3.*pi/2], typically used for correlation studies
double getDeltaPhi(double phiHadron, double phiD)
{
  return RecoDecay::constrainAngle(phiHadron - phiD, -PIHalf);
}

// binning type
using BinningTypeDerived = ColumnBinningPolicy<aod::hf_collisions_reduced::PosZ, aod::hf_collisions_reduced::Multiplicity>;

/// Ds-Hadron correlation pair builder - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
struct HfCorrelatorDsHadronsReduced {
  Produces<aod::DsHadronPair> entryDsHadronPair;
  Produces<aod::DsHadronRecoInfo> entryDsHadronRecoInfo;
  Produces<aod::DsHadronMlInfo> entryDsHadronMlInfo;
  Produces<aod::DsCandRecoInfo> entryDsCandRecoInfo;
  Produces<aod::TrackRecoInfo> entryTrackRecoInfo;
  // Produces<aod::DsHadronGenInfo> entryDsHadronGenInfo;

  Configurable<bool> fillHistoData{"fillHistoData", true, "Flag for filling histograms in data processes"};
  Configurable<int> numberEventsMixed{"numberEventsMixed", 5, "Number of events mixed in ME process"};
  Configurable<std::vector<double>> binsPtD{"binsPtD", std::vector<double>{1., 3., 5., 8., 16., 36.}, "pT bin limits for candidate mass plots"};
  Configurable<std::vector<double>> binsPtHadron{"binsPtHadron", std::vector<double>{0.3, 1., 2., 50.}, "pT bin limits for assoc particle"};

  SliceCache cache;

  // Preslice<aod::AssocTrackReds> tracksPerCol = aod::hf_assoc_track_reduced::hfcRedCollisionId;
  Preslice<aod::AssocTrackReds> tracksPerCol = aod::hf_candidate_reduced::hfcRedCollisionId;
  Preslice<aod::DsCandReduceds> candPerCol = aod::hf_candidate_reduced::hfcRedCollisionId;

  ConfigurableAxis zPoolBins{"zPoolBins", {VARIABLE_WIDTH, -10.0, -2.5, 2.5, 10.0}, "z vertex position pools"};
  ConfigurableAxis multPoolBins{"multPoolBins", {VARIABLE_WIDTH, 0., 900., 1800., 6000.}, "event multiplicity pools (FT0M)"};
  ConfigurableAxis binsMultFT0M{"binsMultFT0M", {100, 0., 10000.}, "Multiplicity as FT0M signal amplitude"};
  ConfigurableAxis binsPosZ{"binsPosZ", {100, -10., 10.}, "primary vertex z coordinate"};
  ConfigurableAxis binsPoolBin{"binsPoolBin", {9, 0., 9.}, "PoolBin"};
  ConfigurableAxis binsEta{"binsEta", {50, -2., 2.}, "#it{#eta}"};
  ConfigurableAxis binsPhi{"binsPhi", {64, -PIHalf, 3. * PIHalf}, "#it{#varphi}"};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    AxisSpec axisMultFT0M = {binsMultFT0M, "MultiplicityFT0M"};
    AxisSpec axisPosZ = {binsPosZ, "PosZ"};
    AxisSpec axisPoolBin = {binsPoolBin, "PoolBin"};
    AxisSpec axisEta = {binsEta, "#it{#eta}"};
    AxisSpec axisPhi = {binsPhi, "#it{#varphi}"};
    AxisSpec axisPtD = {(std::vector<double>)binsPtD, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec axisPtHadron = {(std::vector<double>)binsPtHadron, "#it{p}_{T} Hadron (GeV/#it{c})"};

    // Histograms for data analysis
    if (fillHistoData) {
      registry.add("hMultFT0M", "Multiplicity FT0M", {HistType::kTH1F, {axisMultFT0M}});
      registry.add("hZVtx", "z vertex", {HistType::kTH1F, {axisPosZ}});
      registry.add("hCollisionPoolBin", "Collision pool bin", {HistType::kTH1F, {axisPoolBin}});
      registry.add("hDsPoolBin", "Ds candidates pool bin", {HistType::kTH1F, {axisPoolBin}});
      registry.add("hPhiVsPtCand", "Ds candidates phiVsPt", {HistType::kTH2F, {{axisPhi}, {axisPtD}}});
      registry.add("hPhiVsPtPartAssoc", "Particles associated phiVsPt", {HistType::kTH3F, {{axisPhi}, {axisPtD}, {axisPtHadron}}});
      registry.add("hEtaVsPtCand", "Ds candidates etaVsPt", {HistType::kTH2F, {{axisEta}, {axisPtD}}});
      registry.add("hEtaVsPtPartAssoc", "Particles associated etaVsPt", {HistType::kTH3F, {{axisEta}, {axisPtD}, {axisPtHadron}}});
      registry.add("hTracksPoolBin", "Particles associated pool bin", {HistType::kTH1F, {axisPoolBin}});
    }
  }

  void processDerivedData(aod::HfcRedCollisions const& collisions,
                          soa::Join<aod::DsCandReduceds, aod::DsCandSelInfos> const& candidates,
                          soa::Join<aod::AssocTrackReds, aod::AssocTrackSels> const& tracks)
  {

    BinningTypeDerived corrBinning{{zPoolBins, multPoolBins}, true};

    for (const auto& collision : collisions) {
      int poolBin = corrBinning.getBin(std::make_tuple(collision.posZ(), collision.multiplicity()));
      registry.fill(HIST("hCollisionPoolBin"), poolBin);
      registry.fill(HIST("hMultFT0M"), collision.multiplicity());
      registry.fill(HIST("hZVtx"), collision.posZ());

      auto thisCollId = collision.globalIndex();
      auto candsDsThisColl = candidates.sliceBy(candPerCol, thisCollId);
      auto tracksThisColl = tracks.sliceBy(tracksPerCol, thisCollId);

      for (const auto& candidate : candsDsThisColl) {
        registry.fill(HIST("hDsPoolBin"), poolBin);
        registry.fill(HIST("hPhiVsPtCand"), RecoDecay::constrainAngle(candidate.phiCand(), -PIHalf), candidate.ptCand());
        registry.fill(HIST("hEtaVsPtCand"), candidate.etaCand(), candidate.ptCand());
        entryDsCandRecoInfo(candidate.invMassDs(), candidate.ptCand(), candidate.bdtScorePrompt(), candidate.bdtScoreBkg());
        for (const auto& track : tracksThisColl) {
          registry.fill(HIST("hTracksPoolBin"), poolBin);
          registry.fill(HIST("hPhiVsPtPartAssoc"), RecoDecay::constrainAngle(track.phiAssocTrack(), -PIHalf), candidate.ptCand(), track.ptAssocTrack());
          registry.fill(HIST("hEtaVsPtPartAssoc"), track.etaAssocTrack(), candidate.ptCand(), track.ptAssocTrack());

          entryDsHadronPair(getDeltaPhi(track.phiAssocTrack(), candidate.phiCand()),
                            track.etaAssocTrack() - candidate.etaCand(),
                            candidate.ptCand(),
                            track.ptAssocTrack(),
                            poolBin);
          entryDsHadronRecoInfo(candidate.invMassDs(), false, false);
          entryDsHadronMlInfo(candidate.bdtScorePrompt(), candidate.bdtScoreBkg());
          entryTrackRecoInfo(track.dcaXY(), track.dcaZ(), track.nTpcCrossedRows());
        }
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorDsHadronsReduced, processDerivedData, "Process Derived Data", true);

  void processDerivedDataME(aod::HfcRedCollisions const& collisions,
                            aod::DsCandReduceds const& candidates,
                            aod::AssocTrackReds const& tracks)
  {

    BinningTypeDerived corrBinning{{zPoolBins, multPoolBins}, true};

    for (const auto& collision : collisions) {
      int poolBin = corrBinning.getBin(std::make_tuple(collision.posZ(), collision.multiplicity()));
      registry.fill(HIST("hCollisionPoolBin"), poolBin);
      registry.fill(HIST("hMultFT0M"), collision.multiplicity());
      registry.fill(HIST("hZVtx"), collision.posZ());

      auto thisCollId = collision.globalIndex();
      auto candsDsThisColl = candidates.sliceBy(candPerCol, thisCollId);
      auto tracksThisColl = tracks.sliceBy(tracksPerCol, thisCollId);

      for (const auto& candidate : candsDsThisColl) {
        registry.fill(HIST("hDsPoolBin"), poolBin);
        registry.fill(HIST("hPhiVsPtCand"), RecoDecay::constrainAngle(candidate.phiCand(), -PIHalf), candidate.ptCand());
        registry.fill(HIST("hEtaVsPtCand"), candidate.etaCand(), candidate.ptCand());
        for (const auto& track : tracksThisColl) {
          registry.fill(HIST("hTracksPoolBin"), poolBin);
          registry.fill(HIST("hPhiVsPtPartAssoc"), RecoDecay::constrainAngle(track.phiAssocTrack(), -PIHalf), candidate.ptCand(), track.ptAssocTrack());
          registry.fill(HIST("hEtaVsPtPartAssoc"), track.etaAssocTrack(), candidate.ptCand(), track.ptAssocTrack());
        }
      }
    }

    auto tracksTuple = std::make_tuple(candidates, tracks);

    Pair<aod::HfcRedCollisions, aod::DsCandReduceds, aod::AssocTrackReds, BinningTypeDerived> pairData{corrBinning, numberEventsMixed, -1, collisions, tracksTuple, &cache};

    for (const auto& [c1, tracks1, c2, tracks2] : pairData) {
      if (tracks1.size() == 0) {
        continue;
      }

      int poolBin = corrBinning.getBin(std::make_tuple(c2.posZ(), c2.multiplicity()));
      int poolBinDs = corrBinning.getBin(std::make_tuple(c1.posZ(), c1.multiplicity()));

      if (poolBin != poolBinDs) {
        LOGF(info, "Error, poolBins are diffrent");
      }

      for (const auto& [cand, pAssoc] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        // LOGF(info, "Mixed event tracks pair: (%d, %d) from events (%d, %d), track event: (%d, %d)", cand.index(), pAssoc.index(), c1.index(), c2.index(), cand.hfcRedCollisionId(), pAssoc.hfcRedCollisionId());

        entryDsHadronPair(getDeltaPhi(pAssoc.phiAssocTrack(), cand.phiCand()),
                          pAssoc.etaAssocTrack() - cand.etaCand(),
                          cand.ptCand(),
                          pAssoc.ptAssocTrack(),
                          poolBin);
        entryDsHadronRecoInfo(cand.invMassDs(), false, false);
        // entryDsHadronGenInfo(false, false, 0);
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorDsHadronsReduced, processDerivedDataME, "Process Mixed Event Derived Data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCorrelatorDsHadronsReduced>(cfgc)};
}
