// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file femtoD0Qa.cxx
/// \brief QA task for D0 mesons
/// \author Igor Ptak, WUT, igor.tomasz.ptak@cern.ch

#include "PWGCF/Femto/Core/charmHadronBuilder.h"
#include "PWGCF/Femto/Core/charmHadronHistManager.h"
#include "PWGCF/Femto/Core/collisionBuilder.h"
#include "PWGCF/Femto/Core/collisionHistManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/partitions.h"
#include "PWGCF/Femto/Core/trackHistManager.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
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

struct FemtoD0Qa {

  // setup collisions
  using FemtoCollisions = o2::soa::Join<o2::aod::FCols, o2::aod::FColMasks, o2::aod::FColPos, o2::aod::FColSphericities, o2::aod::FColMults, o2::aod::FColCents>;
  using FilteredFemtoCollisions = o2::soa::Filtered<FemtoCollisions>;
  using FilteredFemtoCollision = FilteredFemtoCollisions::iterator;

  // setup D0s, joined with the mask for the partition and the QA columns
  using FemtoD0s = o2::soa::Join<o2::aod::FD0s, o2::aod::FD0Masks, o2::aod::FD0Extras>;
  // setup tracks with full pid information for the daughter QA
  using FemtoTracks = o2::soa::Join<o2::aod::FTracks, o2::aod::FTrackMass, o2::aod::FTrackDcas, o2::aod::FTrackExtras, o2::aod::FTrackPids>;

  // setup monte carlo, joining the labels that link reco to generated particles
  using FemtoCollisionsWithLabel = o2::soa::Join<FemtoCollisions, o2::aod::FColLabels>;
  using FilteredFemtoCollisionsWithLabel = o2::soa::Filtered<FemtoCollisionsWithLabel>;
  using FilteredFemtoCollisionWithLabel = FilteredFemtoCollisionsWithLabel::iterator;

  using FemtoD0sWithLabel = o2::soa::Join<FemtoD0s, o2::aod::FD0Labels>;
  using FemtoTracksWithLabel = o2::soa::Join<FemtoTracks, o2::aod::FTrackLabels>;
  using FemtoMcParticlesWithLabel = o2::soa::Join<o2::aod::FMcParticles, o2::aod::FMcMotherLabels>;

  o2::framework::SliceCache cache;

  // setup for collisions
  collisionbuilder::ConfCollisionSelection collisionSelection;
  o2::framework::expressions::Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);
  colhistmanager::CollisionHistManager colHistManager;
  colhistmanager::ConfCollisionBinning confCollisionBinning;
  colhistmanager::ConfCollisionQaBinning confCollisionQaBinning;

  // setup for D0s
  charmhadronbuilder::ConfD0Selection1 confD0Selection;

  o2::framework::Partition<FemtoD0s> d0Partition = MAKE_D0_PARTITION(confD0Selection);
  o2::framework::Preslice<FemtoD0s> perColD0s = o2::aod::femtobase::stored::fColId;
  o2::framework::Partition<FemtoD0sWithLabel> d0WithLabelPartition = MAKE_D0_PARTITION(confD0Selection);
  o2::framework::Preslice<FemtoD0sWithLabel> perColD0sWithLabel = o2::aod::femtobase::stored::fColId;

  charmhadronhistmanager::ConfD0Binning1 confD0Binning;
  charmhadronhistmanager::ConfD0QaBinning1 confD0QaBinning;
  charmhadronhistmanager::CharmHadronHistManager<
    charmhadronhistmanager::PrefixD0Qa,
    trackhistmanager::PrefixD01PosDaughterQa,
    trackhistmanager::PrefixD01NegDaughterQa,
    modes::CharmHadron::kD0>
    d0HistManager;

  // setup for daughters
  trackhistmanager::ConfD01PosDauBinning confD01PosDaughterBinning;
  trackhistmanager::ConfD01PosDauQaBinning confD01PosDaughterQaBinning;

  trackhistmanager::ConfD01NegDauBinning confD01NegDaughterBinning;
  trackhistmanager::ConfD01NegDauQaBinning confD01NegDaughterQaBinning;

  o2::framework::HistogramRegistry hRegistry{"FemtoD0Qa", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    bool processData = doprocessD0;

    std::map<colhistmanager::ColHist, std::vector<o2::framework::AxisSpec>> colHistSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> posDaughterHistSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> negDaughterHistSpec;
    std::map<charmhadronhistmanager::CharmHadronHist, std::vector<o2::framework::AxisSpec>> d0HistSpec;
    std::map<charmhadronhistmanager::CharmHadronHist, std::vector<o2::framework::AxisSpec>> d0QaHistSpec;

    if (processData) {
      colHistSpec = colhistmanager::makeColQaHistSpecMap(confCollisionBinning, confCollisionQaBinning);
      colHistManager.init<modes::Mode::kReco_Qa>(&hRegistry, colHistSpec, confCollisionBinning, confCollisionQaBinning);
      posDaughterHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confD01PosDaughterBinning, confD01PosDaughterQaBinning);
      negDaughterHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confD01NegDaughterBinning, confD01NegDaughterQaBinning);
      d0HistSpec = charmhadronhistmanager::makeD0HistSpecMap(confD0Binning);
      d0QaHistSpec = charmhadronhistmanager::makeD0QaHistSpecMap(confD0QaBinning);
      d0HistManager.init<modes::Mode::kReco_Qa>(&hRegistry, d0HistSpec, d0QaHistSpec, confD0Selection, confD0QaBinning, posDaughterHistSpec, confD01PosDaughterQaBinning, negDaughterHistSpec, confD01NegDaughterQaBinning);
    } else {
      colHistSpec = colhistmanager::makeColMcQaHistSpecMap(confCollisionBinning, confCollisionQaBinning);
      colHistManager.init<modes::Mode::kReco_Qa_Mc>(&hRegistry, colHistSpec, confCollisionBinning, confCollisionQaBinning);
      posDaughterHistSpec = trackhistmanager::makeTrackMcQaHistSpecMap(confD01PosDaughterBinning, confD01PosDaughterQaBinning);
      negDaughterHistSpec = trackhistmanager::makeTrackMcQaHistSpecMap(confD01NegDaughterBinning, confD01NegDaughterQaBinning);
      d0HistSpec = charmhadronhistmanager::makeD0McQaHistSpecMap(confD0Binning, confD0QaBinning);
      d0QaHistSpec = charmhadronhistmanager::makeD0QaHistSpecMap(confD0QaBinning);
      d0HistManager.init<modes::Mode::kReco_Qa_Mc>(&hRegistry, d0HistSpec, d0QaHistSpec, confD0Selection, confD0QaBinning, posDaughterHistSpec, confD01PosDaughterQaBinning, negDaughterHistSpec, confD01NegDaughterQaBinning);
    }

    hRegistry.print();
  };

  void processD0(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoD0s const& /*d0s*/)
  {
    auto d0Slice = d0Partition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    if (d0Slice.size() == 0) {
      return;
    }
    colHistManager.fill<modes::Mode::kReco_Qa>(col);
    for (auto const& d0 : d0Slice) {
      d0HistManager.fill<modes::Mode::kReco_Qa>(d0, tracks);
    }
  }
  PROCESS_SWITCH(FemtoD0Qa, processD0, "Process D0s", true);

  void processD0Mc(FilteredFemtoCollisionWithLabel const& col, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoD0sWithLabel const& /*d0s*/, FemtoMcParticlesWithLabel const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    auto d0Slice = d0WithLabelPartition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    if (d0Slice.size() == 0) {
      return;
    }
    colHistManager.fill<modes::Mode::kReco_Qa_Mc>(col, mcCols);
    for (auto const& d0 : d0Slice) {
      d0HistManager.fill<modes::Mode::kReco_Qa_Mc>(d0, tracks, col, mcParticles, mcMothers, mcPartonicMothers);
    }
  }
  PROCESS_SWITCH(FemtoD0Qa, processD0Mc, "Process D0s with MC information", false);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& context)
{
  o2::framework::WorkflowSpec workflow{
    adaptAnalysisTask<FemtoD0Qa>(context),
  };
  return workflow;
}
