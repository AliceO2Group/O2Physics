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

/// \file femtoCharmHadronQa.cxx
/// \brief QA task for charm hadrons
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

struct FemtoCharmHadronQa {
  // setup collisions
  using FemtoCollisions = o2::soa::Join<o2::aod::FCols, o2::aod::FColMasks, o2::aod::FColPos, o2::aod::FColSphericities, o2::aod::FColMults, o2::aod::FColCents>;
  using FilteredFemtoCollisions = o2::soa::Filtered<FemtoCollisions>;
  using FilteredFemtoCollision = FilteredFemtoCollisions::iterator;

  // setup D0s, joined with the mask for the partition and the QA columns
  using FemtoD0s = o2::soa::Join<o2::aod::FD0s, o2::aod::FD0Masks, o2::aod::FD0Extras>;
  // setup Lcs, joined with the mask for the partition and the QA columns
  using FemtoLcs = o2::soa::Join<o2::aod::FLcs, o2::aod::FLcMasks, o2::aod::FLcExtras>;
  // setup tracks with full pid information for the daughter QA
  using FemtoTracks = o2::soa::Join<o2::aod::FTracks, o2::aod::FTrackMass, o2::aod::FTrackDcas, o2::aod::FTrackExtras, o2::aod::FTrackPids>;

  // setup monte carlo, joining the labels that link reco to generated particles
  using FemtoCollisionsWithLabel = o2::soa::Join<FemtoCollisions, o2::aod::FColLabels>;
  using FilteredFemtoCollisionsWithLabel = o2::soa::Filtered<FemtoCollisionsWithLabel>;
  using FilteredFemtoCollisionWithLabel = FilteredFemtoCollisionsWithLabel::iterator;

  using FemtoD0sWithLabel = o2::soa::Join<FemtoD0s, o2::aod::FD0Labels>;
  using FemtoLcsWithLabel = o2::soa::Join<FemtoLcs, o2::aod::FLcLabels>;
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

  // setup for Lcs
  charmhadronbuilder::ConfLcSelection1 confLcSelection;

  o2::framework::Partition<FemtoLcs> lcPartition = MAKE_CHARM3PRONG_PARTITION(confLcSelection);
  o2::framework::Preslice<FemtoLcs> perColLcs = o2::aod::femtobase::stored::fColId;
  o2::framework::Partition<FemtoLcsWithLabel> lcWithLabelPartition = MAKE_CHARM3PRONG_PARTITION(confLcSelection);
  o2::framework::Preslice<FemtoLcsWithLabel> perColLcsWithLabel = o2::aod::femtobase::stored::fColId;

  charmhadronhistmanager::ConfLcBinning1 confLcBinning;
  charmhadronhistmanager::ConfLcQaBinning1 confLcQaBinning;

  // D0
  charmhadronhistmanager::CharmHadronHistManager<
    charmhadronhistmanager::PrefixD0Qa,
    trackhistmanager::PrefixD01PosDaughterQa,
    trackhistmanager::PrefixD01NegDaughterQa,
    trackhistmanager::PrefixD01NegDaughterQa, // unused for 2-prong candidates
    modes::CharmHadron::kD0>
    d0HistManager;
  // Lc
  charmhadronhistmanager::CharmHadronHistManager<
    charmhadronhistmanager::PrefixLcQa,
    trackhistmanager::PrefixLc1ProtonDaughterQa,
    trackhistmanager::PrefixLc1KaonDaughterQa,
    trackhistmanager::PrefixLc1PionDaughterQa,
    modes::CharmHadron::kLc>
    lcHistManager;

  // setup for daughters
  trackhistmanager::ConfD01PosDauBinning confD01PosDaughterBinning;
  trackhistmanager::ConfD01PosDauQaBinning confD01PosDaughterQaBinning;

  trackhistmanager::ConfD01NegDauBinning confD01NegDaughterBinning;
  trackhistmanager::ConfD01NegDauQaBinning confD01NegDaughterQaBinning;

  trackhistmanager::ConfLc1ProtonDauBinning confLc1ProtonDaughterBinning;
  trackhistmanager::ConfLc1ProtonDauQaBinning confLc1ProtonDaughterQaBinning;

  trackhistmanager::ConfLc1KaonDauBinning confLc1KaonDaughterBinning;
  trackhistmanager::ConfLc1KaonDauQaBinning confLc1KaonDaughterQaBinning;

  trackhistmanager::ConfLc1PionDauBinning confLc1PionDaughterBinning;
  trackhistmanager::ConfLc1PionDauQaBinning confLc1PionDaughterQaBinning;

  o2::framework::HistogramRegistry hRegistry{"FemtoCharmHadronQa", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    if ((static_cast<int>(doprocessD0) + static_cast<int>(doprocessD0Mc) + static_cast<int>(doprocessLc) + static_cast<int>(doprocessLcMc)) != 1) {
      LOG(fatal) << "Only one process can be activated";
    }
    bool processData = doprocessD0 || doprocessLc;

    std::map<colhistmanager::ColHist, std::vector<o2::framework::AxisSpec>> colHistSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> posDaughterHistSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> negDaughterHistSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> protonHistSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> kaonHistSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> pionHistSpec;
    std::map<charmhadronhistmanager::CharmHadronHist, std::vector<o2::framework::AxisSpec>> d0HistSpec;
    std::map<charmhadronhistmanager::CharmHadronHist, std::vector<o2::framework::AxisSpec>> d0QaHistSpec;
    std::map<charmhadronhistmanager::CharmHadronHist, std::vector<o2::framework::AxisSpec>> lcHistSpec;
    std::map<charmhadronhistmanager::CharmHadronHist, std::vector<o2::framework::AxisSpec>> lcQaHistSpec;

    if (processData) {
      colHistSpec = colhistmanager::makeColQaHistSpecMap(confCollisionBinning, confCollisionQaBinning);
      colHistManager.init<modes::Mode::kReco_Qa>(&hRegistry, colHistSpec, confCollisionBinning, confCollisionQaBinning);

      if (doprocessD0) {
        posDaughterHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confD01PosDaughterBinning, confD01PosDaughterQaBinning);
        negDaughterHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confD01NegDaughterBinning, confD01NegDaughterQaBinning);
        d0HistSpec = charmhadronhistmanager::makeCharmHadronHistSpecMap(confD0Binning);
        d0QaHistSpec = charmhadronhistmanager::makeD0QaHistSpecMap(confD0QaBinning);
        d0HistManager.init<modes::Mode::kReco_Qa>(&hRegistry, d0HistSpec, d0QaHistSpec, confD0Selection, confD0QaBinning, posDaughterHistSpec, confD01PosDaughterQaBinning, negDaughterHistSpec, confD01NegDaughterQaBinning);
      }
      if (doprocessLc) {
        protonHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confLc1ProtonDaughterBinning, confLc1ProtonDaughterQaBinning);
        kaonHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confLc1KaonDaughterBinning, confLc1KaonDaughterQaBinning);
        pionHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confLc1PionDaughterBinning, confLc1PionDaughterQaBinning);
        lcHistSpec = charmhadronhistmanager::makeCharmHadronHistSpecMap(confLcBinning);
        lcQaHistSpec = charmhadronhistmanager::makeLcQaHistSpecMap(confLcQaBinning);
        lcHistManager.init<modes::Mode::kReco_Qa>(&hRegistry, lcHistSpec, lcQaHistSpec, confLcSelection, confLcQaBinning, protonHistSpec, confLc1ProtonDaughterQaBinning, kaonHistSpec, confLc1KaonDaughterQaBinning, pionHistSpec, confLc1PionDaughterQaBinning);
      }
    } else {
      colHistSpec = colhistmanager::makeColMcQaHistSpecMap(confCollisionBinning, confCollisionQaBinning);
      colHistManager.init<modes::Mode::kReco_Qa_Mc>(&hRegistry, colHistSpec, confCollisionBinning, confCollisionQaBinning);

      if (doprocessD0Mc) {
        posDaughterHistSpec = trackhistmanager::makeTrackMcQaHistSpecMap(confD01PosDaughterBinning, confD01PosDaughterQaBinning);
        negDaughterHistSpec = trackhistmanager::makeTrackMcQaHistSpecMap(confD01NegDaughterBinning, confD01NegDaughterQaBinning);
        d0HistSpec = charmhadronhistmanager::makeD0McQaHistSpecMap(confD0Binning, confD0QaBinning);
        d0QaHistSpec = charmhadronhistmanager::makeD0QaHistSpecMap(confD0QaBinning);
        d0HistManager.init<modes::Mode::kReco_Qa_Mc>(&hRegistry, d0HistSpec, d0QaHistSpec, confD0Selection, confD0QaBinning, posDaughterHistSpec, confD01PosDaughterQaBinning, negDaughterHistSpec, confD01NegDaughterQaBinning);
      }
      if (doprocessLcMc) {
        protonHistSpec = trackhistmanager::makeTrackMcQaHistSpecMap(confLc1ProtonDaughterBinning, confLc1ProtonDaughterQaBinning);
        kaonHistSpec = trackhistmanager::makeTrackMcQaHistSpecMap(confLc1KaonDaughterBinning, confLc1KaonDaughterQaBinning);
        pionHistSpec = trackhistmanager::makeTrackMcQaHistSpecMap(confLc1PionDaughterBinning, confLc1PionDaughterQaBinning);
        lcHistSpec = charmhadronhistmanager::makeLcMcQaHistSpecMap(confLcBinning, confLcQaBinning);
        lcQaHistSpec = charmhadronhistmanager::makeLcQaHistSpecMap(confLcQaBinning);
        lcHistManager.init<modes::Mode::kReco_Qa_Mc>(&hRegistry, lcHistSpec, lcQaHistSpec, confLcSelection, confLcQaBinning, protonHistSpec, confLc1ProtonDaughterQaBinning, kaonHistSpec, confLc1KaonDaughterQaBinning, pionHistSpec, confLc1PionDaughterQaBinning);
      }
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
  PROCESS_SWITCH(FemtoCharmHadronQa, processD0, "Process D0s", true);

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
  PROCESS_SWITCH(FemtoCharmHadronQa, processD0Mc, "Process D0s with MC information", false);

  void processLc(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoLcs const& /*lcs*/)
  {
    auto lcSlice = lcPartition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    if (lcSlice.size() == 0) {
      return;
    }
    colHistManager.fill<modes::Mode::kReco_Qa>(col);
    for (auto const& lc : lcSlice) {
      lcHistManager.fill<modes::Mode::kReco_Qa>(lc, tracks);
    }
  }
  PROCESS_SWITCH(FemtoCharmHadronQa, processLc, "Process Lcs", false);

  void processLcMc(FilteredFemtoCollisionWithLabel const& col, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoLcsWithLabel const& /*lcs*/, FemtoMcParticlesWithLabel const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    auto lcSlice = lcWithLabelPartition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    if (lcSlice.size() == 0) {
      return;
    }
    colHistManager.fill<modes::Mode::kReco_Qa_Mc>(col, mcCols);
    for (auto const& lc : lcSlice) {
      lcHistManager.fill<modes::Mode::kReco_Qa_Mc>(lc, tracks, col, mcParticles, mcMothers, mcPartonicMothers);
    }
  }
  PROCESS_SWITCH(FemtoCharmHadronQa, processLcMc, "Process Lcs with MC information", false);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& context)
{
  o2::framework::WorkflowSpec workflow{
    adaptAnalysisTask<FemtoCharmHadronQa>(context),
  };
  return workflow;
}
