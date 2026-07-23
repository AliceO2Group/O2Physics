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

/// \file femtoMcParticleQa.cxx
/// \brief QA task for tracks
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#include "PWGCF/Femto/Core/collisionHistManager.h"
#include "PWGCF/Femto/Core/mcBuilder.h"
#include "PWGCF/Femto/Core/mcParticleHistManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/particleCleaner.h"
#include "PWGCF/Femto/Core/partitions.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <map>
#include <vector>

using namespace o2::analysis::femto;

struct FemtoMcParticleQa {

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

  o2::framework::Partition<FemtoMcParticles> mcParticlesPartition1 = MAKE_MC_PARTICLE_PARTITION(confMcParticleSelection1);
  o2::framework::Preslice<FemtoMcParticles> perColReco = o2::aod::femtomcparticle::fMcColId;

  particlecleaner::ConfMcParticleCleaner1 confMcParticleCleaner1;
  particlecleaner::ParticleCleaner mcParticleCleaner;

  o2::framework::HistogramRegistry hRegistry{"FemtoMcParticleQA", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    mcParticleCleaner.init(confMcParticleCleaner1);

    std::map<colhistmanager::ColHist, std::vector<o2::framework::AxisSpec>> colHistSpec;
    std::map<mcparticlehistmanager::McParticleHist, std::vector<o2::framework::AxisSpec>> mcParticleHistSpec;

    colHistSpec = colhistmanager::makeColMcHistSpecMap(confCollisionBinning);
    colHistManager.init<modes::Mode::kMc>(&hRegistry, colHistSpec, confCollisionBinning);
    mcParticleHistSpec = mcparticlehistmanager::makeMcParticleHistSpecMap(confMcParticleBinning1);
    mcParticleHistManager1.init(&hRegistry, mcParticleHistSpec, confMcParticleBinning1);

    hRegistry.print();
  };

  void process(FilteredFemtoMcCollision const& col, FemtoMcParticles const& /*mcParticles*/, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    auto mcParticleSlice = mcParticlesPartition1->sliceByCached(o2::aod::femtomcparticle::fMcColId, col.globalIndex(), cache);
    if (mcParticleSlice.size() == 0) {
      return;
    }
    colHistManager.fill<modes::Mode::kMc>(col);
    for (auto const& mcParticle : mcParticleSlice) {
      if (!mcParticleCleaner.isClean(mcParticle, mcMothers, mcPartonicMothers)) {
        continue;
      }
      mcParticleHistManager1.fill(mcParticle, mcMothers, mcPartonicMothers);
    }
  }
};

o2::framework::WorkflowSpec
  defineDataProcessing(o2::framework::ConfigContext const& context)
{
  o2::framework::WorkflowSpec workflow{
    adaptAnalysisTask<FemtoMcParticleQa>(context),
  };
  return workflow;
}
