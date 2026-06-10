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
/// \brief Tasks that computes correlation between tracks and resonances decaying into two tracks
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#include "PWGCF/Femto/Core/closePairRejection.h"
#include "PWGCF/Femto/Core/collisionBuilder.h"
#include "PWGCF/Femto/Core/collisionHistManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/pairBuilder.h"
#include "PWGCF/Femto/Core/pairHistManager.h"
#include "PWGCF/Femto/Core/partitions.h"
#include "PWGCF/Femto/Core/trackBuilder.h"
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

  o2::framework::Partition<FemtoKstar0s> k0shortPartition = MAKE_K0SHORT_PARTITION(k0shortSelection);
  o2::framework::Preslice<FemtoKstar0s> perColk0shorts = o2::aod::femtobase::stored::fColId;

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

  // rho0s
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

    // setup columnpolicy for binning
    // default values are used during instantiation, so we need to explicity update them here
    mixBinsVtxMult = {{confMixing.vtxBins, confMixing.multBins.value}, true};
    mixBinsVtxCent = {{confMixing.vtxBins.value, confMixing.centBins.value}, true};
    mixBinsVtxMultCent = {{confMixing.vtxBins.value, confMixing.multBins.value, confMixing.centBins.value}, true};

    // setup histograms
    auto colHistSpec = colhistmanager::makeColHistSpecMap(confCollisionBinning);

    auto lambdaHistSpec = v0histmanager::makeV0HistSpecMap(confLambdaBinning);
    auto v0PosDauSpec = trackhistmanager::makeTrackHistSpecMap(confV0PosDauBinning);
    auto v0NegDauSpec = trackhistmanager::makeTrackHistSpecMap(confV0NegDauBinning);

    auto phiHistSpec = twotrackresonancehistmanager::makeTwoTrackResonanceHistSpecMap(confPhiBinning);
    auto resoPosDauSpec = trackhistmanager::makeTrackHistSpecMap(confResoPosDauBinning);
    auto resov0NegDauSpec = trackhistmanager::makeTrackHistSpecMap(confResoNegDauBinning);

    auto cprHistSpecPos = closepairrejection::makeCprHistSpecMap(confCprPos);
    auto cprHistSpecNeg = closepairrejection::makeCprHistSpecMap(confCprNeg);

    auto pairV0TwoTrackResonanceHistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning, confMixing);

    pairLambdaPhiBuilder.init<modes::Mode::kSe_Analysis, modes::Mode::kMe_Analysis>(&hRegistry, confCollisionBinning, lambdaSelection, phiSelection, confCprPos, confCprNeg, confMixing, confPairBinning, confPairCuts, colHistSpec, lambdaHistSpec, v0PosDauSpec, v0NegDauSpec, phiHistSpec, resoPosDauSpec, resov0NegDauSpec, pairV0TwoTrackResonanceHistSpec, cprHistSpecPos, cprHistSpecNeg);
  };

  void processLambdaPhiSameEvent(FilteredCollision const& col, FemtoTracks const& tracks, FemtoLambdas const& /*lambdas*/, FemtoPhis const& /*phis*/)
  {
    pairLambdaPhiBuilder.processSameEvent<modes::Mode::kSe_Analysis>(col, tracks, lambdaPartition, phiPartition, cache);
  }
  PROCESS_SWITCH(FemtoPairV0TwoTrackResonance, processLambdaPhiSameEvent, "Enable processing same event processing for lambdas and phis", true);

  void processLambdaPhiMixedEvent(FilteredCollisions const& cols, FemtoTracks const& tracks, FemtoLambdas const& /*lambda*/, FemtoPhis const& /*phis*/)
  {
    pairLambdaPhiBuilder.processMixedEvent<modes::Mode::kMe_Analysis>(cols, tracks, lambdaPartition, phiPartition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairV0TwoTrackResonance, processLambdaPhiMixedEvent, "Enable processing mixed event processing for lambdas and phis", true);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  o2::framework::WorkflowSpec workflow{
    adaptAnalysisTask<FemtoPairV0TwoTrackResonance>(cfgc),
  };
  return workflow;
}
