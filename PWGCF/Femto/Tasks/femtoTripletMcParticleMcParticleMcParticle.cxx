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

/// \file femtoTripletMcParticleMcParticleMcParticle.cxx
/// \brief Task for triplets of generated particles (mc truth only)
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#include "PWGCF/Femto/Core/closePairRejection.h"
#include "PWGCF/Femto/Core/closeTripletRejection.h"
#include "PWGCF/Femto/Core/collisionHistManager.h"
#include "PWGCF/Femto/Core/mcBuilder.h"
#include "PWGCF/Femto/Core/mcParticleHistManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/pairCleaner.h"
#include "PWGCF/Femto/Core/particleCleaner.h"
#include "PWGCF/Femto/Core/partitions.h"
#include "PWGCF/Femto/Core/tripletBuilder.h"
#include "PWGCF/Femto/Core/tripletCleaner.h"
#include "PWGCF/Femto/Core/tripletHistManager.h"
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

struct FemtoTripletMcParticleMcParticleMcParticle {

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

  // setup mc particles
  mcbuilder::ConfMcParticleSelection1 confMcParticleSelection1;
  mcparticlehistmanager::ConfMcParticleBinning1 confMcParticleBinning1;
  particlecleaner::ConfMcParticleCleaner1 confMcParticleCleaner1;

  mcbuilder::ConfMcParticleSelection2 confMcParticleSelection2;
  mcparticlehistmanager::ConfMcParticleBinning2 confMcParticleBinning2;
  particlecleaner::ConfMcParticleCleaner2 confMcParticleCleaner2;

  mcbuilder::ConfMcParticleSelection3 confMcParticleSelection3;
  mcparticlehistmanager::ConfMcParticleBinning3 confMcParticleBinning3;
  particlecleaner::ConfMcParticleCleaner3 confMcParticleCleaner3;

  o2::framework::Partition<FemtoMcParticles> mcParticlesPartition1 = MAKE_MC_PARTICLE_PARTITION(confMcParticleSelection1);
  o2::framework::Partition<FemtoMcParticles> mcParticlesPartition2 = MAKE_MC_PARTICLE_PARTITION(confMcParticleSelection2);
  o2::framework::Partition<FemtoMcParticles> mcParticlesPartition3 = MAKE_MC_PARTICLE_PARTITION(confMcParticleSelection3);

  // setup triplets
  triplethistmanager::ConfTripletBinning confTripletBinning;
  triplethistmanager::ConfTripletCuts confTripletCuts;
  paircleaner::ConfPairCleanerBinning confTripletCleanerBinning;

  closetripletrejection::ConfCtrMcParticleMcParticleMcParticle confCtr;

  tripletbuilder::TripletMcParticleMcParticleMcParticleBuilder<
    mcparticlehistmanager::PrefixMcParticle1,
    mcparticlehistmanager::PrefixMcParticle2,
    mcparticlehistmanager::PrefixMcParticle3,
    triplethistmanager::PrefixMcParticleMcParticleMcParticleSe,
    triplethistmanager::PrefixMcParticleMcParticleMcParticleMe,
    closetripletrejection::PrefixMcParticle1McParticle2Se,
    closetripletrejection::PrefixMcParticle2McParticle3Se,
    closetripletrejection::PrefixMcParticle1McParticle3Se,
    closetripletrejection::PrefixMcParticle1McParticle2Me,
    closetripletrejection::PrefixMcParticle2McParticle3Me,
    closetripletrejection::PrefixMcParticle1McParticle3Me>
    tripletMcParticleMcParticleMcParticleBuilder;

  // setup mixing
  std::vector<double> defaultVtxBins{10, -10, 10};
  std::vector<double> defaultMultBins{50, 0, 200};
  std::vector<double> defaultCentBins{10, 0, 100};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Mult> mixBinsVtxMult{{defaultVtxBins, defaultMultBins}, true};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Cent> mixBinsVtxCent{{defaultVtxBins, defaultCentBins}, true};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Mult, o2::aod::femtocollisions::Cent> mixBinsVtxMultCent{{defaultVtxBins, defaultMultBins, defaultCentBins}, true};
  triplethistmanager::ConfMixing confMixing;

  o2::framework::HistogramRegistry hRegistry{"FemtoMcParticleMcParticleMcParticle", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    // setup columnpolicy for binning
    // default values are used during instantiation, so we need to explicity update them here
    mixBinsVtxMult = {{confMixing.vtxBins.value, confMixing.multBins.value}, true};
    mixBinsVtxCent = {{confMixing.vtxBins.value, confMixing.centBins.value}, true};
    mixBinsVtxMultCent = {{confMixing.vtxBins.value, confMixing.multBins.value, confMixing.centBins.value}, true};

    // setup histogram specs
    std::map<colhistmanager::ColHist, std::vector<o2::framework::AxisSpec>> colHistSpec = colhistmanager::makeColMcHistSpecMap(confCollisionBinning);
    std::map<mcparticlehistmanager::McParticleHist, std::vector<o2::framework::AxisSpec>> mcParticleHistSpec1 = mcparticlehistmanager::makeMcParticleHistSpecMap(confMcParticleBinning1);
    std::map<mcparticlehistmanager::McParticleHist, std::vector<o2::framework::AxisSpec>> mcParticleHistSpec2 = mcparticlehistmanager::makeMcParticleHistSpecMap(confMcParticleBinning2);
    std::map<mcparticlehistmanager::McParticleHist, std::vector<o2::framework::AxisSpec>> mcParticleHistSpec3 = mcparticlehistmanager::makeMcParticleHistSpecMap(confMcParticleBinning3);
    std::map<triplethistmanager::TripletHist, std::vector<o2::framework::AxisSpec>> tripletHistSpec = triplethistmanager::makeTripletMcTruthHistSpecMap(confTripletBinning, confMixing);
    std::map<closepairrejection::CprHist, std::vector<o2::framework::AxisSpec>> ctrHistSpec = closepairrejection::makeCprHistSpecMap(confCtr);
    std::map<paircleaner::PairCleanerHist, std::vector<o2::framework::AxisSpec>> tripletCleanerHistSpec = paircleaner::makePairCleanerHistSpecMap(confTripletCleanerBinning);

    tripletMcParticleMcParticleMcParticleBuilder.init<modes::Mode::kSe_Mc, modes::Mode::kMe_Mc>(&hRegistry, confCollisionBinning, confMcParticleSelection1, confMcParticleSelection2, confMcParticleSelection3, confMcParticleBinning1, confMcParticleBinning2, confMcParticleBinning3, confMcParticleCleaner1, confMcParticleCleaner2, confMcParticleCleaner3, confCtr, confMixing, confTripletBinning, confTripletCuts, colHistSpec, mcParticleHistSpec1, mcParticleHistSpec2, mcParticleHistSpec3, tripletHistSpec, ctrHistSpec, tripletCleanerHistSpec);

    hRegistry.print();
  };

  void processSameEvent(FilteredFemtoMcCollision const& col, FemtoMcParticles const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    tripletMcParticleMcParticleMcParticleBuilder.processSameEvent<modes::Mode::kSe_Mc>(col, mcParticles, mcMothers, mcPartonicMothers, mcParticlesPartition1, mcParticlesPartition2, mcParticlesPartition3, cache);
  }
  PROCESS_SWITCH(FemtoTripletMcParticleMcParticleMcParticle, processSameEvent, "Enable processing same event processing", true);

  void processMixedEvent(FilteredFemtoMcCollisions const& cols, FemtoMcParticles const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    tripletMcParticleMcParticleMcParticleBuilder.processMixedEvent<modes::Mode::kMe_Mc>(cols, mcParticles, mcMothers, mcPartonicMothers, mcParticlesPartition1, mcParticlesPartition2, mcParticlesPartition3, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoTripletMcParticleMcParticleMcParticle, processMixedEvent, "Enable processing mixed event processing", true);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& context)
{
  o2::framework::WorkflowSpec workflow{
    adaptAnalysisTask<FemtoTripletMcParticleMcParticleMcParticle>(context),
  };
  return workflow;
}
