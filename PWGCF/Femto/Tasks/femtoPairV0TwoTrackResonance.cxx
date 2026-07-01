// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file femtoPairV0TwoTrackResonance.cxx
/// \brief Tasks that computes correlation between v0s and resonances decaying into two tracks
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#include "PWGCF/Femto/Core/closePairRejection.h"
#include "PWGCF/Femto/Core/collisionBuilder.h"
#include "PWGCF/Femto/Core/collisionHistManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/pairBuilder.h"
#include "PWGCF/Femto/Core/pairHistManager.h"
#include "PWGCF/Femto/Core/particleCleaner.h"
#include "PWGCF/Femto/Core/partitions.h"
#include "PWGCF/Femto/Core/trackHistManager.h"
#include "PWGCF/Femto/Core/twoTrackResonanceBuilder.h"
#include "PWGCF/Femto/Core/twoTrackResonanceHistManager.h"
#include "PWGCF/Femto/Core/v0Builder.h"
#include "PWGCF/Femto/Core/v0HistManager.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <string>
#include <vector>

using namespace o2::analysis::femto;

struct FemtoPairV0TwoTrackResonance {

  // setup tables
  using Collisions = o2::soa::Join<o2::aod::FCols, o2::aod::FColMasks>;
  using Collision = Collisions::iterator;

  using FilteredCollisions = o2::soa::Filtered<Collisions>;
  using FilteredCollision = FilteredCollisions::iterator;

  using FemtoTracks = o2::soa::Join<o2::aod::FTracks, o2::aod::FTrackMasks>;

  using FemtoLambdas = o2::soa::Join<o2::aod::FLambdas, o2::aod::FLambdaMasks>;
  using FemtoK0shorts = o2::soa::Join<o2::aod::FK0shorts, o2::aod::FK0shortMasks>;

  using FemtoPhis = o2::soa::Join<o2::aod::FPhis, o2::aod::FPhiMasks>;
  using FemtoKstar0s = o2::soa::Join<o2::aod::FKstar0s, o2::aod::FKstar0Masks>;
  using FemtoRho0s = o2::soa::Join<o2::aod::FRho0s, o2::aod::FRho0Masks>;

  o2::framework::SliceCache cache;

  // setup collisions
  collisionbuilder::ConfCollisionSelection collisionSelection;
  o2::framework::expressions::Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);
  colhistmanager::ConfCollisionBinning confCollisionBinning;

  // setup for daughters for v0s
  trackhistmanager::ConfV0PosDauBinning confV0PosDauBinning;
  trackhistmanager::ConfV0NegDauBinning confV0NegDauBinning;

  // setup lambdas
  v0builder::ConfLambdaSelection1 lambdaSelection;
  v0histmanager::ConfLambdaBinning1 confLambdaBinning;
  particlecleaner::ConfLambdaCleaner1 confLambdaCleaner;

  o2::framework::Partition<FemtoLambdas> lambdaPartition = MAKE_LAMBDA_PARTITION(lambdaSelection);
  o2::framework::Preslice<FemtoLambdas> perColLambdas = o2::aod::femtobase::stored::fColId;

  // setup k0shorts
  v0builder::ConfK0shortSelection1 k0shortSelection;
  v0histmanager::ConfK0shortBinning1 confK0shortBinning;
  particlecleaner::ConfK0shortCleaner1 confK0shortCleaner;

  o2::framework::Partition<FemtoK0shorts> k0shortPartition = MAKE_K0SHORT_PARTITION(k0shortSelection);
  o2::framework::Preslice<FemtoK0shorts> perColK0shorts = o2::aod::femtobase::stored::fColId;

  // setup for daughters for resonances
  trackhistmanager::ConfResonancePosDauBinning confResoPosDauBinning;
  trackhistmanager::ConfResonanceNegDauBinning confResoNegDauBinning;

  // setup phis
  twotrackresonancebuilder::ConfPhiSelection phiSelection;
  twotrackresonancehistmanager::ConfPhiBinning confPhiBinning;
  o2::framework::Partition<FemtoPhis> phiPartition = MAKE_RESONANCE_0_PARTITON(phiSelection);
  o2::framework::Preslice<FemtoPhis> perColPhis = o2::aod::femtobase::stored::fColId;

  // setup kstar0
  twotrackresonancebuilder::ConfKstar0Selection kstar0Selection;
  twotrackresonancehistmanager::ConfKstar0Binning confKstar0Binning;
  o2::framework::Partition<FemtoKstar0s> kstar0Partition = MAKE_RESONANCE_1_PARTITON(kstar0Selection);
  o2::framework::Preslice<FemtoKstar0s> perColKstar0s = o2::aod::femtobase::stored::fColId;

  // setup rho0s
  twotrackresonancebuilder::ConfRho0Selection rho0Selection;
  twotrackresonancehistmanager::ConfRho0Binning confRho0Binning;
  o2::framework::Partition<FemtoRho0s> rho0Partition = MAKE_RESONANCE_0_PARTITON(rho0Selection);
  o2::framework::Preslice<FemtoRho0s> perColRho0s = o2::aod::femtobase::stored::fColId;

  // setup pairs
  pairhistmanager::ConfPairBinning confPairBinning;
  pairhistmanager::ConfPairCuts confPairCuts;

  // setup for lambda-phi pairs
  pairbuilder::PairV0TwoTrackResonanceBuilder<
    v0histmanager::PrefixLambda1,
    trackhistmanager::PrefixV01PosDaughter,
    trackhistmanager::PrefixV01NegDaughter,
    twotrackresonancehistmanager::PrefixPhi,
    trackhistmanager::PrefixResonancePosDaughter,
    trackhistmanager::PrefixResonanceNegDaughter,
    pairhistmanager::PrefixV0ResonanceSe,
    pairhistmanager::PrefixV0ResonanceMe,
    closepairrejection::PrefixV0TwoTrackResonancePosSe,
    closepairrejection::PrefixV0TwoTrackResonanceNegSe,
    closepairrejection::PrefixV0TwoTrackResonancePosMe,
    closepairrejection::PrefixV0TwoTrackResonanceNegMe,
    modes::V0::kLambda,
    modes::TwoTrackResonance::kPhi>
    pairLambdaPhiBuilder;

  // setup for lambda-kstar0 pairs
  pairbuilder::PairV0TwoTrackResonanceBuilder<
    v0histmanager::PrefixLambda1,
    trackhistmanager::PrefixV01PosDaughter,
    trackhistmanager::PrefixV01NegDaughter,
    twotrackresonancehistmanager::PrefixKstar,
    trackhistmanager::PrefixResonancePosDaughter,
    trackhistmanager::PrefixResonanceNegDaughter,
    pairhistmanager::PrefixV0ResonanceSe,
    pairhistmanager::PrefixV0ResonanceMe,
    closepairrejection::PrefixV0TwoTrackResonancePosSe,
    closepairrejection::PrefixV0TwoTrackResonanceNegSe,
    closepairrejection::PrefixV0TwoTrackResonancePosMe,
    closepairrejection::PrefixV0TwoTrackResonanceNegMe,
    modes::V0::kLambda,
    modes::TwoTrackResonance::kKstar0>
    pairLambdaKstar0Builder;

  // setup for lambda-rho0 pairs
  pairbuilder::PairV0TwoTrackResonanceBuilder<
    v0histmanager::PrefixLambda1,
    trackhistmanager::PrefixV01PosDaughter,
    trackhistmanager::PrefixV01NegDaughter,
    twotrackresonancehistmanager::PrefixRho,
    trackhistmanager::PrefixResonancePosDaughter,
    trackhistmanager::PrefixResonanceNegDaughter,
    pairhistmanager::PrefixV0ResonanceSe,
    pairhistmanager::PrefixV0ResonanceMe,
    closepairrejection::PrefixV0TwoTrackResonancePosSe,
    closepairrejection::PrefixV0TwoTrackResonanceNegSe,
    closepairrejection::PrefixV0TwoTrackResonancePosMe,
    closepairrejection::PrefixV0TwoTrackResonanceNegMe,
    modes::V0::kLambda,
    modes::TwoTrackResonance::kRho0>
    pairLambdaRho0Builder;

  // setup for k0short-phi pairs
  pairbuilder::PairV0TwoTrackResonanceBuilder<
    v0histmanager::PrefixK0short1,
    trackhistmanager::PrefixV01PosDaughter,
    trackhistmanager::PrefixV01NegDaughter,
    twotrackresonancehistmanager::PrefixPhi,
    trackhistmanager::PrefixResonancePosDaughter,
    trackhistmanager::PrefixResonanceNegDaughter,
    pairhistmanager::PrefixV0ResonanceSe,
    pairhistmanager::PrefixV0ResonanceMe,
    closepairrejection::PrefixV0TwoTrackResonancePosSe,
    closepairrejection::PrefixV0TwoTrackResonanceNegSe,
    closepairrejection::PrefixV0TwoTrackResonancePosMe,
    closepairrejection::PrefixV0TwoTrackResonanceNegMe,
    modes::V0::kK0short,
    modes::TwoTrackResonance::kPhi>
    pairK0shortPhiBuilder;

  // setup for k0short-kstar0 pairs
  pairbuilder::PairV0TwoTrackResonanceBuilder<
    v0histmanager::PrefixK0short1,
    trackhistmanager::PrefixV01PosDaughter,
    trackhistmanager::PrefixV01NegDaughter,
    twotrackresonancehistmanager::PrefixKstar,
    trackhistmanager::PrefixResonancePosDaughter,
    trackhistmanager::PrefixResonanceNegDaughter,
    pairhistmanager::PrefixV0ResonanceSe,
    pairhistmanager::PrefixV0ResonanceMe,
    closepairrejection::PrefixV0TwoTrackResonancePosSe,
    closepairrejection::PrefixV0TwoTrackResonanceNegSe,
    closepairrejection::PrefixV0TwoTrackResonancePosMe,
    closepairrejection::PrefixV0TwoTrackResonanceNegMe,
    modes::V0::kK0short,
    modes::TwoTrackResonance::kKstar0>
    pairK0shortKstar0Builder;

  // setup for k0short-rho0 pairs
  pairbuilder::PairV0TwoTrackResonanceBuilder<
    v0histmanager::PrefixK0short1,
    trackhistmanager::PrefixV01PosDaughter,
    trackhistmanager::PrefixV01NegDaughter,
    twotrackresonancehistmanager::PrefixRho,
    trackhistmanager::PrefixResonancePosDaughter,
    trackhistmanager::PrefixResonanceNegDaughter,
    pairhistmanager::PrefixV0ResonanceSe,
    pairhistmanager::PrefixV0ResonanceMe,
    closepairrejection::PrefixV0TwoTrackResonancePosSe,
    closepairrejection::PrefixV0TwoTrackResonanceNegSe,
    closepairrejection::PrefixV0TwoTrackResonancePosMe,
    closepairrejection::PrefixV0TwoTrackResonanceNegMe,
    modes::V0::kK0short,
    modes::TwoTrackResonance::kRho0>
    pairK0shortRho0Builder;

  // setup mixing
  std::vector<double> defaultVtxBins{10, -10, 10};
  std::vector<double> defaultMultBins{50, 0, 200};
  std::vector<double> defaultCentBins{10, 0, 100};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Mult> mixBinsVtxMult{{defaultVtxBins, defaultMultBins}, true};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Cent> mixBinsVtxCent{{defaultVtxBins, defaultCentBins}, true};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Mult, o2::aod::femtocollisions::Cent> mixBinsVtxMultCent{{defaultVtxBins, defaultMultBins, defaultCentBins}, true};

  pairhistmanager::ConfMixing confMixing;

  // setup cpr
  closepairrejection::ConfCprV0DaughterResoDaughterPos confCprPos;
  closepairrejection::ConfCprV0DaughterResoDaughterNeg confCprNeg;

  o2::framework::HistogramRegistry hRegistry{"FemtoV0TwoTrackResonance", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {

    // all pair combinations share the same histogram prefixes, so only one combination can be processed per workflow instance
    if ((doprocessLambdaPhiSameEvent || doprocessLambdaPhiMixedEvent) +
          (doprocessLambdaKstar0SameEvent || doprocessLambdaKstar0MixedEvent) +
          (doprocessLambdaRho0SameEvent || doprocessLambdaRho0MixedEvent) +
          (doprocessK0shortPhiSameEvent || doprocessK0shortPhiMixedEvent) +
          (doprocessK0shortKstar0SameEvent || doprocessK0shortKstar0MixedEvent) +
          (doprocessK0shortRho0SameEvent || doprocessK0shortRho0MixedEvent) >
        1) {
      LOG(fatal) << "Can only process one v0-resonance combination (lambda/k0short x phi/kstar0/rho0)";
    }

    // setup columnpolicy for binning
    // default values are used during instantiation, so we need to explicity update them here
    mixBinsVtxMult = {{confMixing.vtxBins.value, confMixing.multBins.value}, true};
    mixBinsVtxCent = {{confMixing.vtxBins.value, confMixing.centBins.value}, true};
    mixBinsVtxMultCent = {{confMixing.vtxBins.value, confMixing.multBins.value, confMixing.centBins.value}, true};

    // setup histograms
    auto colHistSpec = colhistmanager::makeColHistSpecMap(confCollisionBinning);

    auto v0PosDauSpec = trackhistmanager::makeTrackHistSpecMap(confV0PosDauBinning);
    auto v0NegDauSpec = trackhistmanager::makeTrackHistSpecMap(confV0NegDauBinning);

    auto resoPosDauSpec = trackhistmanager::makeTrackHistSpecMap(confResoPosDauBinning);
    auto resoNegDauSpec = trackhistmanager::makeTrackHistSpecMap(confResoNegDauBinning);

    auto cprHistSpecPos = closepairrejection::makeCprHistSpecMap(confCprPos);
    auto cprHistSpecNeg = closepairrejection::makeCprHistSpecMap(confCprNeg);

    auto pairV0TwoTrackResonanceHistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning, confMixing);

    // setup for lambda-phi
    if (doprocessLambdaPhiSameEvent || doprocessLambdaPhiMixedEvent) {
      auto lambdaHistSpec = v0histmanager::makeV0HistSpecMap(confLambdaBinning);
      auto phiHistSpec = twotrackresonancehistmanager::makeTwoTrackResonanceHistSpecMap(confPhiBinning);
      pairLambdaPhiBuilder.init<modes::Mode::kSe_Reco, modes::Mode::kMe_Reco>(&hRegistry, confCollisionBinning, lambdaSelection, phiSelection, confCprPos, confCprNeg, confMixing, confPairBinning, confPairCuts, colHistSpec, lambdaHistSpec, v0PosDauSpec, v0NegDauSpec, phiHistSpec, resoPosDauSpec, resoNegDauSpec, pairV0TwoTrackResonanceHistSpec, cprHistSpecPos, cprHistSpecNeg);
    }

    // setup for lambda-kstar0
    if (doprocessLambdaKstar0SameEvent || doprocessLambdaKstar0MixedEvent) {
      auto lambdaHistSpec = v0histmanager::makeV0HistSpecMap(confLambdaBinning);
      auto kstar0HistSpec = twotrackresonancehistmanager::makeTwoTrackResonanceHistSpecMap(confKstar0Binning);
      pairLambdaKstar0Builder.init<modes::Mode::kSe_Reco, modes::Mode::kMe_Reco>(&hRegistry, confCollisionBinning, lambdaSelection, kstar0Selection, confCprPos, confCprNeg, confMixing, confPairBinning, confPairCuts, colHistSpec, lambdaHistSpec, v0PosDauSpec, v0NegDauSpec, kstar0HistSpec, resoPosDauSpec, resoNegDauSpec, pairV0TwoTrackResonanceHistSpec, cprHistSpecPos, cprHistSpecNeg);
    }

    // setup for lambda-rho0
    if (doprocessLambdaRho0SameEvent || doprocessLambdaRho0MixedEvent) {
      auto lambdaHistSpec = v0histmanager::makeV0HistSpecMap(confLambdaBinning);
      auto rho0HistSpec = twotrackresonancehistmanager::makeTwoTrackResonanceHistSpecMap(confRho0Binning);
      pairLambdaRho0Builder.init<modes::Mode::kSe_Reco, modes::Mode::kMe_Reco>(&hRegistry, confCollisionBinning, lambdaSelection, rho0Selection, confCprPos, confCprNeg, confMixing, confPairBinning, confPairCuts, colHistSpec, lambdaHistSpec, v0PosDauSpec, v0NegDauSpec, rho0HistSpec, resoPosDauSpec, resoNegDauSpec, pairV0TwoTrackResonanceHistSpec, cprHistSpecPos, cprHistSpecNeg);
    }

    // setup for k0short-phi
    if (doprocessK0shortPhiSameEvent || doprocessK0shortPhiMixedEvent) {
      auto k0shortHistSpec = v0histmanager::makeV0HistSpecMap(confK0shortBinning);
      auto phiHistSpec = twotrackresonancehistmanager::makeTwoTrackResonanceHistSpecMap(confPhiBinning);
      pairK0shortPhiBuilder.init<modes::Mode::kSe_Reco, modes::Mode::kMe_Reco>(&hRegistry, confCollisionBinning, k0shortSelection, phiSelection, confCprPos, confCprNeg, confMixing, confPairBinning, confPairCuts, colHistSpec, k0shortHistSpec, v0PosDauSpec, v0NegDauSpec, phiHistSpec, resoPosDauSpec, resoNegDauSpec, pairV0TwoTrackResonanceHistSpec, cprHistSpecPos, cprHistSpecNeg);
    }

    // setup for k0short-kstar0
    if (doprocessK0shortKstar0SameEvent || doprocessK0shortKstar0MixedEvent) {
      auto k0shortHistSpec = v0histmanager::makeV0HistSpecMap(confK0shortBinning);
      auto kstar0HistSpec = twotrackresonancehistmanager::makeTwoTrackResonanceHistSpecMap(confKstar0Binning);
      pairK0shortKstar0Builder.init<modes::Mode::kSe_Reco, modes::Mode::kMe_Reco>(&hRegistry, confCollisionBinning, k0shortSelection, kstar0Selection, confCprPos, confCprNeg, confMixing, confPairBinning, confPairCuts, colHistSpec, k0shortHistSpec, v0PosDauSpec, v0NegDauSpec, kstar0HistSpec, resoPosDauSpec, resoNegDauSpec, pairV0TwoTrackResonanceHistSpec, cprHistSpecPos, cprHistSpecNeg);
    }

    // setup for k0short-rho0
    if (doprocessK0shortRho0SameEvent || doprocessK0shortRho0MixedEvent) {
      auto k0shortHistSpec = v0histmanager::makeV0HistSpecMap(confK0shortBinning);
      auto rho0HistSpec = twotrackresonancehistmanager::makeTwoTrackResonanceHistSpecMap(confRho0Binning);
      pairK0shortRho0Builder.init<modes::Mode::kSe_Reco, modes::Mode::kMe_Reco>(&hRegistry, confCollisionBinning, k0shortSelection, rho0Selection, confCprPos, confCprNeg, confMixing, confPairBinning, confPairCuts, colHistSpec, k0shortHistSpec, v0PosDauSpec, v0NegDauSpec, rho0HistSpec, resoPosDauSpec, resoNegDauSpec, pairV0TwoTrackResonanceHistSpec, cprHistSpecPos, cprHistSpecNeg);
    }
  }

  // lambda-phi
  void processLambdaPhiSameEvent(FilteredCollision const& col, FemtoTracks const& tracks, FemtoLambdas const& /*lambdas*/, FemtoPhis const& /*phis*/)
  {
    pairLambdaPhiBuilder.processSameEvent<modes::Mode::kSe_Reco>(col, tracks, lambdaPartition, phiPartition, cache);
  }
  PROCESS_SWITCH(FemtoPairV0TwoTrackResonance, processLambdaPhiSameEvent, "Enable processing same event processing for lambdas and phis", true);

  void processLambdaPhiMixedEvent(FilteredCollisions const& cols, FemtoTracks const& tracks, FemtoLambdas const& /*lambdas*/, FemtoPhis const& /*phis*/)
  {
    pairLambdaPhiBuilder.processMixedEvent<modes::Mode::kMe_Reco>(cols, tracks, lambdaPartition, phiPartition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairV0TwoTrackResonance, processLambdaPhiMixedEvent, "Enable processing mixed event processing for lambdas and phis", true);

  // lambda-kstar0
  void processLambdaKstar0SameEvent(FilteredCollision const& col, FemtoTracks const& tracks, FemtoLambdas const& /*lambdas*/, FemtoKstar0s const& /*kstar0s*/)
  {
    pairLambdaKstar0Builder.processSameEvent<modes::Mode::kSe_Reco>(col, tracks, lambdaPartition, kstar0Partition, cache);
  }
  PROCESS_SWITCH(FemtoPairV0TwoTrackResonance, processLambdaKstar0SameEvent, "Enable processing same event processing for lambdas and kstar0s", false);

  void processLambdaKstar0MixedEvent(FilteredCollisions const& cols, FemtoTracks const& tracks, FemtoLambdas const& /*lambdas*/, FemtoKstar0s const& /*kstar0s*/)
  {
    pairLambdaKstar0Builder.processMixedEvent<modes::Mode::kMe_Reco>(cols, tracks, lambdaPartition, kstar0Partition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairV0TwoTrackResonance, processLambdaKstar0MixedEvent, "Enable processing mixed event processing for lambdas and kstar0s", false);

  // lambda-rho0
  void processLambdaRho0SameEvent(FilteredCollision const& col, FemtoTracks const& tracks, FemtoLambdas const& /*lambdas*/, FemtoRho0s const& /*rho0s*/)
  {
    pairLambdaRho0Builder.processSameEvent<modes::Mode::kSe_Reco>(col, tracks, lambdaPartition, rho0Partition, cache);
  }
  PROCESS_SWITCH(FemtoPairV0TwoTrackResonance, processLambdaRho0SameEvent, "Enable processing same event processing for lambdas and rho0s", false);

  void processLambdaRho0MixedEvent(FilteredCollisions const& cols, FemtoTracks const& tracks, FemtoLambdas const& /*lambdas*/, FemtoRho0s const& /*rho0s*/)
  {
    pairLambdaRho0Builder.processMixedEvent<modes::Mode::kMe_Reco>(cols, tracks, lambdaPartition, rho0Partition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairV0TwoTrackResonance, processLambdaRho0MixedEvent, "Enable processing mixed event processing for lambdas and rho0s", false);

  // k0short-phi
  void processK0shortPhiSameEvent(FilteredCollision const& col, FemtoTracks const& tracks, FemtoK0shorts const& /*k0shorts*/, FemtoPhis const& /*phis*/)
  {
    pairK0shortPhiBuilder.processSameEvent<modes::Mode::kSe_Reco>(col, tracks, k0shortPartition, phiPartition, cache);
  }
  PROCESS_SWITCH(FemtoPairV0TwoTrackResonance, processK0shortPhiSameEvent, "Enable processing same event processing for k0shorts and phis", false);

  void processK0shortPhiMixedEvent(FilteredCollisions const& cols, FemtoTracks const& tracks, FemtoK0shorts const& /*k0shorts*/, FemtoPhis const& /*phis*/)
  {
    pairK0shortPhiBuilder.processMixedEvent<modes::Mode::kMe_Reco>(cols, tracks, k0shortPartition, phiPartition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairV0TwoTrackResonance, processK0shortPhiMixedEvent, "Enable processing mixed event processing for k0shorts and phis", false);

  // k0short-kstar0
  void processK0shortKstar0SameEvent(FilteredCollision const& col, FemtoTracks const& tracks, FemtoK0shorts const& /*k0shorts*/, FemtoKstar0s const& /*kstar0s*/)
  {
    pairK0shortKstar0Builder.processSameEvent<modes::Mode::kSe_Reco>(col, tracks, k0shortPartition, kstar0Partition, cache);
  }
  PROCESS_SWITCH(FemtoPairV0TwoTrackResonance, processK0shortKstar0SameEvent, "Enable processing same event processing for k0shorts and kstar0s", false);

  void processK0shortKstar0MixedEvent(FilteredCollisions const& cols, FemtoTracks const& tracks, FemtoK0shorts const& /*k0shorts*/, FemtoKstar0s const& /*kstar0s*/)
  {
    pairK0shortKstar0Builder.processMixedEvent<modes::Mode::kMe_Reco>(cols, tracks, k0shortPartition, kstar0Partition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairV0TwoTrackResonance, processK0shortKstar0MixedEvent, "Enable processing mixed event processing for k0shorts and kstar0s", false);

  // k0short-rho0
  void processK0shortRho0SameEvent(FilteredCollision const& col, FemtoTracks const& tracks, FemtoK0shorts const& /*k0shorts*/, FemtoRho0s const& /*rho0s*/)
  {
    pairK0shortRho0Builder.processSameEvent<modes::Mode::kSe_Reco>(col, tracks, k0shortPartition, rho0Partition, cache);
  }
  PROCESS_SWITCH(FemtoPairV0TwoTrackResonance, processK0shortRho0SameEvent, "Enable processing same event processing for k0shorts and rho0s", false);

  void processK0shortRho0MixedEvent(FilteredCollisions const& cols, FemtoTracks const& tracks, FemtoK0shorts const& /*k0shorts*/, FemtoRho0s const& /*rho0s*/)
  {
    pairK0shortRho0Builder.processMixedEvent<modes::Mode::kMe_Reco>(cols, tracks, k0shortPartition, rho0Partition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairV0TwoTrackResonance, processK0shortRho0MixedEvent, "Enable processing mixed event processing for k0shorts and rho0s", false);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  o2::framework::WorkflowSpec workflow{
    adaptAnalysisTask<FemtoPairV0TwoTrackResonance>(cfgc),
  };
  return workflow;
}
