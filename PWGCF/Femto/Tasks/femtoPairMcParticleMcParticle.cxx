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

/// \file femtoPairMcParticleMcParticle.cxx
/// \brief Tasks that computes correlation between two tracks
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#include "PWGCF/Femto/Core/closePairRejection.h"
#include "PWGCF/Femto/Core/collisionHistManager.h"
#include "PWGCF/Femto/Core/mcBuilder.h"
#include "PWGCF/Femto/Core/mcParticleHistManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/pairBuilder.h"
#include "PWGCF/Femto/Core/pairHistManager.h"
#include "PWGCF/Femto/Core/particleCleaner.h"
#include "PWGCF/Femto/Core/partitions.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <map>
#include <vector>

using namespace o2::analysis::femto;

struct FemtoPairMcParticleMcParticle {

  // setup tables
  using FemtoMcCollisions = o2::aod::FMcCols;
  using FilteredFemtoMcCollisions = o2::soa::Filtered<FemtoMcCollisions>;
  using FilteredFemtoMcCollision = FilteredFemtoMcCollisions::iterator;

  using FemtoMcParticles = o2::soa::Join<o2::aod::FMcParticles, o2::aod::FMcMotherLabels>;

  o2::framework::SliceCache cache;

  // setup collisions
  mcbuilder::ConfMcCollisionFilters collisionSelection;
  o2::framework::expressions::Filter collisionFilter = MAKE_MC_COLLISION_FILTER(collisionSelection);
  colhistmanager::ConfCollisionBinning confCollisionBinning;
  colhistmanager::CollisionHistManager colHistManager;

  // setup mc particles
  mcbuilder::ConfMcParticleSelection1 confMcParticleSelection1;
  mcparticlehistmanager::ConfMcParticleBinning1 confMcParticleBinning1;
  mcparticlehistmanager::McParticleHistManager<mcparticlehistmanager::PrefixMcParticle1> mcParticleHistManager1;
  particlecleaner::ConfMcParticleCleaner1 confMcParticleCleaner1;
  particlecleaner::ParticleCleaner mcParticleCleaner1;

  o2::framework::Partition<FemtoMcParticles> mcParticlesPartition1 = MAKE_MC_PARTICLE_PARTITION(confMcParticleSelection1);

  mcbuilder::ConfMcParticleSelection2 confMcParticleSelection2;
  mcparticlehistmanager::ConfMcParticleBinning2 confMcParticleBinning2;
  mcparticlehistmanager::McParticleHistManager<mcparticlehistmanager::PrefixMcParticle2> mcParticleHistManager2;
  particlecleaner::ConfMcParticleCleaner2 confMcParticleCleaner2;
  particlecleaner::ParticleCleaner mcParticleCleaner2;

  o2::framework::Partition<FemtoMcParticles> mcParticlesPartition2 = MAKE_MC_PARTICLE_PARTITION(confMcParticleSelection1);

  o2::framework::Preslice<FemtoMcParticles> perColParticles = o2::aod::femtomcparticle::fMcColId;

  // setup pairs
  pairhistmanager::ConfPairBinning confPairBinning;
  pairhistmanager::ConfPairCuts confPairCuts;

  closepairrejection::ConfCprTrackTrack confCpr;

  pairbuilder::PairMcParticleMcParticleBuilder<
    mcparticlehistmanager::PrefixMcParticle1,
    mcparticlehistmanager::PrefixMcParticle2,
    pairhistmanager::PrefixMcParticleMcParticleSe,
    pairhistmanager::PrefixMcParticleMcParticleMe,
    closepairrejection::PrefixMcParticleMcParticleSe,
    closepairrejection::PrefixMcParticleMcParticleMe>
    pairMcParticleMcParticleBuilder;

  // setup mixing
  std::vector<double> defaultVtxBins{10, -10, 10};
  std::vector<double> defaultMultBins{50, 0, 200};
  std::vector<double> defaultCentBins{10, 0, 100};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Mult> mixBinsVtxMult{{defaultVtxBins, defaultMultBins}, true};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Cent> mixBinsVtxCent{{defaultVtxBins, defaultCentBins}, true};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Mult, o2::aod::femtocollisions::Cent> mixBinsVtxMultCent{{defaultVtxBins, defaultMultBins, defaultCentBins}, true};
  pairhistmanager::ConfMixing confMixing;

  o2::framework::HistogramRegistry hRegistry{"FemtoMcParticleMcParticle", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    // setup columnpolicy for binning
    // default values are used during instantiation, so we need to explicity update them here
    mixBinsVtxMult = {{confMixing.vtxBins, confMixing.multBins.value}, true};
    mixBinsVtxCent = {{confMixing.vtxBins.value, confMixing.centBins.value}, true};
    mixBinsVtxMultCent = {{confMixing.vtxBins.value, confMixing.multBins.value, confMixing.centBins.value}, true};

    // setup histogram specs
    std::map<colhistmanager::ColHist, std::vector<o2::framework::AxisSpec>> colHistSpec;
    std::map<mcparticlehistmanager::McParticleHist, std::vector<o2::framework::AxisSpec>> mcParticleHistSpec1;
    std::map<mcparticlehistmanager::McParticleHist, std::vector<o2::framework::AxisSpec>> mcParticleHistSpec2;
    std::map<pairhistmanager::PairHist, std::vector<o2::framework::AxisSpec>> pairHistSpec;
    std::map<closepairrejection::CprHist, std::vector<o2::framework::AxisSpec>> cprHistSpec = closepairrejection::makeCprHistSpecMap(confCpr);

    colHistSpec = colhistmanager::makeColMcHistSpecMap(confCollisionBinning);
    mcParticleHistSpec1 = mcparticlehistmanager::makeMcParticleHistSpecMap(confMcParticleBinning1);
    mcParticleHistSpec2 = mcparticlehistmanager::makeMcParticleHistSpecMap(confMcParticleBinning2);
    pairHistSpec = pairhistmanager::makePairMcTruthHistSpecMap(confPairBinning, confMixing);
    pairMcParticleMcParticleBuilder.init<modes::Mode::kSe_Mc, modes::Mode::kMe_Mc>(&hRegistry, confCollisionBinning, confMcParticleSelection1, confMcParticleSelection2, confMcParticleBinning1, confMcParticleBinning2, confMcParticleCleaner1, confMcParticleCleaner2, confCpr, confMixing, confPairBinning, confPairCuts, colHistSpec, mcParticleHistSpec1, mcParticleHistSpec2, pairHistSpec, cprHistSpec);

    hRegistry.print();
  };

  void processSameEvent(FilteredFemtoMcCollision const& col, FemtoMcParticles const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairMcParticleMcParticleBuilder.processSameEvent<modes::Mode::kSe_Mc>(col, mcParticles, mcMothers, mcPartonicMothers, mcParticlesPartition1, mcParticlesPartition2, cache);
  }
  PROCESS_SWITCH(FemtoPairMcParticleMcParticle, processSameEvent, "Enable processing same event processing", true);

  void processMixedEvent(FilteredFemtoMcCollisions const& cols, FemtoMcParticles const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairMcParticleMcParticleBuilder.processMixedEvent<modes::Mode::kMe_Mc>(cols, mcParticles, mcMothers, mcPartonicMothers, mcParticlesPartition1, mcParticlesPartition2, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairMcParticleMcParticle, processMixedEvent, "Enable processing mixed event processing", true);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& context)
{
  o2::framework::WorkflowSpec workflow{
    adaptAnalysisTask<FemtoPairMcParticleMcParticle>(context),
  };
  return workflow;
}
