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

/// \file femtoCascadeQa.cxx
/// \brief Tasks for Qa of cascades
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#include "PWGCF/Femto/Core/cascadeBuilder.h"
#include "PWGCF/Femto/Core/cascadeHistManager.h"
#include "PWGCF/Femto/Core/collisionBuilder.h"
#include "PWGCF/Femto/Core/collisionHistManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/particleCleaner.h"
#include "PWGCF/Femto/Core/partitions.h"
#include "PWGCF/Femto/Core/trackHistManager.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/Expressions.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/InitContext.h"
#include "Framework/OutputObjHeader.h"
#include "Framework/runDataProcessing.h"

using namespace o2::analysis::femto;

struct FemtoCascadeQa {

  using FemtoCollisions = o2::soa::Join<o2::aod::FCols, o2::aod::FColMasks, o2::aod::FColPos, o2::aod::FColSphericities, o2::aod::FColMults>;
  using FilteredFemtoCollisions = o2::soa::Filtered<FemtoCollisions>;
  using FilteredFemtoCollision = FilteredFemtoCollisions::iterator;

  using FemtoCollisionsWithLabel = o2::soa::Join<FemtoCollisions, o2::aod::FColLabels>;
  using FilteredFemtoCollisionsWithLabel = o2::soa::Filtered<FemtoCollisionsWithLabel>;
  using FilteredFemtoCollisionWithLabel = FilteredFemtoCollisionsWithLabel::iterator;

  using FemtoXis = o2::soa::Join<o2::aod::FXis, o2::aod::FXiMasks, o2::aod::FXiExtras>;
  using FemtoOmegas = o2::soa::Join<o2::aod::FOmegas, o2::aod::FOmegaMasks, o2::aod::FOmegaExtras>;
  using FemtoTracks = o2::soa::Join<o2::aod::FTracks, o2::aod::FTrackDcas, o2::aod::FTrackExtras, o2::aod::FTrackPids>;

  using FemtoXisWithLabel = o2::soa::Join<FemtoXis, o2::aod::FXiLabels>;
  using FemtoOmegasWithLabel = o2::soa::Join<FemtoOmegas, o2::aod::FOmegaLabels>;
  using FemtoTracksWithLabel = o2::soa::Join<FemtoTracks, o2::aod::FTrackLabels>;

  o2::framework::SliceCache cache;

  // setup collisions
  collisionbuilder::ConfCollisionSelection collisionSelection;
  o2::framework::expressions::Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);
  colhistmanager::CollisionHistManager colHistManager;
  colhistmanager::ConfCollisionBinning confCollisionBinning;
  colhistmanager::ConfCollisionQaBinning confCollisionQaBinning;

  // setup for xis
  cascadebuilder::ConfXiSelection confXiSelection;

  particlecleaner::ConfXiCleaner1 confXiCleaner;
  particlecleaner::ParticleCleaner xiCleaner;

  o2::framework::Partition<FemtoXis> xiPartition = MAKE_CASCADE_PARTITION(confXiSelection);
  o2::framework::Preslice<FemtoXis> preColXis = o2::aod::femtobase::stored::fColId;

  o2::framework::Partition<FemtoXisWithLabel> xiWithLabelPartition = MAKE_CASCADE_PARTITION(confXiSelection);
  o2::framework::Preslice<FemtoXisWithLabel> perColXisWithLabel = o2::aod::femtobase::stored::fColId;

  cascadehistmanager::ConfXiBinning confXiBinning;
  cascadehistmanager::ConfXiQaBinning confXiQaBinning;
  cascadehistmanager::CascadeHistManager<
    cascadehistmanager::PrefixXiQa,
    trackhistmanager::PrefixCascadeBachelorQa,
    trackhistmanager::PrefixV0PosDaughterQa,
    trackhistmanager::PrefixV0NegDaughterQa,
    modes::Cascade::kXi>
    xiHistManager;

  // setup for omegas
  cascadebuilder::ConfOmegaSelection confOmegaSelection;

  particlecleaner::ConfOmegaCleaner1 confOmegaCleaner;
  particlecleaner::ParticleCleaner omegaCleaner;

  o2::framework::Partition<FemtoOmegas> omegaPartition = MAKE_CASCADE_PARTITION(confOmegaSelection);
  o2::framework::Preslice<FemtoOmegas> preColOmegas = o2::aod::femtobase::stored::fColId;

  o2::framework::Partition<FemtoOmegasWithLabel> omegaWithLabelPartition = MAKE_CASCADE_PARTITION(confOmegaSelection);
  o2::framework::Preslice<FemtoOmegasWithLabel> perColOmegasWithLabel = o2::aod::femtobase::stored::fColId;

  cascadehistmanager::ConfOmegaBinning confOmegaBinning;
  cascadehistmanager::ConfOmegaQaBinning confOmegaQaBinning;
  cascadehistmanager::CascadeHistManager<
    cascadehistmanager::PrefixOmegaQa,
    trackhistmanager::PrefixCascadeBachelorQa,
    trackhistmanager::PrefixV0PosDaughterQa,
    trackhistmanager::PrefixV0NegDaughterQa,
    modes::Cascade::kOmega>
    omegaHistManager;

  // setup for daughters/bachelor
  trackhistmanager::ConfCascadePosDauBinning confPosDaughterBinning;
  trackhistmanager::ConfCascadePosDauQaBinning confPosDaughterQaBinning;
  trackhistmanager::ConfCascadeNegDauBinning confNegDaughterBinning;
  trackhistmanager::ConfCascadeNegDauQaBinning confNegDaughterQaBinning;
  trackhistmanager::ConfCascadeBachelorBinning confBachelorBinning;
  trackhistmanager::ConfCascadeBachelorQaBinning confBachelorQaBinning;

  o2::framework::HistogramRegistry hRegistry{"FemtoCascadeQA", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    if ((doprocessXi + doprocessXiMc + doprocessOmega + doprocessOmegaMc) > 1) {
      LOG(fatal) << "Only one process can be activated";
    }
    bool processData = doprocessXi || doprocessOmega;

    xiCleaner.init(confXiCleaner);
    omegaCleaner.init(confOmegaCleaner);

    std::map<colhistmanager::ColHist, std::vector<o2::framework::AxisSpec>> colHistSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> bachelorHistSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> posDaughterHistSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> negDaughterHistSpec;
    std::map<cascadehistmanager::CascadeHist, std::vector<o2::framework::AxisSpec>> xiHistSpec;
    std::map<cascadehistmanager::CascadeHist, std::vector<o2::framework::AxisSpec>> omegaHistSpec;

    if (processData) {
      colHistSpec = colhistmanager::makeColQaHistSpecMap(confCollisionBinning, confCollisionQaBinning);
      colHistManager.init<modes::Mode::kAnalysis_Qa>(&hRegistry, colHistSpec, confCollisionQaBinning);
      bachelorHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confBachelorBinning, confBachelorQaBinning);
      posDaughterHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confPosDaughterBinning, confPosDaughterQaBinning);
      negDaughterHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confNegDaughterBinning, confNegDaughterQaBinning);
      if (doprocessXi) {
        xiHistSpec = cascadehistmanager::makeCascadeQaHistSpecMap(confXiBinning, confXiQaBinning);
        xiHistManager.init<modes::Mode::kAnalysis_Qa>(&hRegistry, xiHistSpec, confXiSelection, confXiQaBinning, bachelorHistSpec, confBachelorQaBinning, posDaughterHistSpec, confPosDaughterQaBinning, negDaughterHistSpec, confNegDaughterQaBinning);
      }
      if (doprocessOmega) {
        omegaHistSpec = cascadehistmanager::makeCascadeQaHistSpecMap(confOmegaBinning, confOmegaQaBinning);
        omegaHistManager.init<modes::Mode::kAnalysis_Qa>(&hRegistry, omegaHistSpec, confOmegaSelection, confOmegaQaBinning, bachelorHistSpec, confBachelorQaBinning, posDaughterHistSpec, confPosDaughterQaBinning, negDaughterHistSpec, confNegDaughterQaBinning);
      }
    } else {
      colHistSpec = colhistmanager::makeColMcQaHistSpecMap(confCollisionBinning, confCollisionQaBinning);
      colHistManager.init<modes::Mode::kAnalysis_Qa_Mc>(&hRegistry, colHistSpec, confCollisionQaBinning);
      bachelorHistSpec = trackhistmanager::makeTrackMcQaHistSpecMap(confBachelorBinning, confBachelorQaBinning);
      posDaughterHistSpec = trackhistmanager::makeTrackMcQaHistSpecMap(confPosDaughterBinning, confPosDaughterQaBinning);
      negDaughterHistSpec = trackhistmanager::makeTrackMcQaHistSpecMap(confNegDaughterBinning, confNegDaughterQaBinning);
      if (doprocessXiMc) {
        xiHistSpec = cascadehistmanager::makeCascadeMcQaHistSpecMap(confXiBinning, confXiQaBinning);
        xiHistManager.init<modes::Mode::kAnalysis_Qa_Mc>(&hRegistry, xiHistSpec, confXiSelection, confXiQaBinning, bachelorHistSpec, confBachelorQaBinning, posDaughterHistSpec, confPosDaughterQaBinning, negDaughterHistSpec, confNegDaughterQaBinning);
      }
      if (doprocessOmegaMc) {
        omegaHistSpec = cascadehistmanager::makeCascadeMcQaHistSpecMap(confOmegaBinning, confOmegaQaBinning);
        omegaHistManager.init<modes::Mode::kAnalysis_Qa_Mc>(&hRegistry, omegaHistSpec, confOmegaSelection, confOmegaQaBinning, bachelorHistSpec, confBachelorQaBinning, posDaughterHistSpec, confPosDaughterQaBinning, negDaughterHistSpec, confNegDaughterQaBinning);
      }
    }
    hRegistry.print();
  };

  void processXi(FilteredFemtoCollision const& col, FemtoXis const& /*xis*/, FemtoTracks const& tracks)
  {
    colHistManager.fill<modes::Mode::kAnalysis_Qa>(col);
    auto xiSlice = xiPartition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    for (auto const& xi : xiSlice) {
      xiHistManager.fill<modes::Mode::kAnalysis_Qa>(xi, tracks);
    }
  }
  PROCESS_SWITCH(FemtoCascadeQa, processXi, "Process Xis", true);

  void processXiMc(FilteredFemtoCollisionWithLabel const& col, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoXisWithLabel const& /*xis*/, o2::aod::FMcParticles const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    colHistManager.fill<modes::Mode::kAnalysis_Qa_Mc>(col, mcCols);
    auto xiSlice = xiWithLabelPartition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    for (auto const& xi : xiSlice) {
      if (!xiCleaner.isClean(xi, mcParticles, mcMothers, mcPartonicMothers)) {
        continue;
      }
      xiHistManager.fill<modes::Mode::kAnalysis_Qa_Mc>(xi, tracks, mcParticles, mcMothers, mcPartonicMothers);
    }
  }
  PROCESS_SWITCH(FemtoCascadeQa, processXiMc, "Process Xis with MC information", false);

  void processOmega(FilteredFemtoCollision const& col, FemtoOmegas const& /*omegas*/, FemtoTracks const& tracks)
  {
    colHistManager.fill<modes::Mode::kAnalysis_Qa>(col);
    auto omegaSlice = omegaPartition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    for (auto const& omega : omegaSlice) {
      omegaHistManager.fill<modes::Mode::kAnalysis_Qa>(omega, tracks);
    }
  }
  PROCESS_SWITCH(FemtoCascadeQa, processOmega, "Process Omegas", false);

  void processOmegaMc(FilteredFemtoCollisionWithLabel const& col, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoOmegasWithLabel const& /*omegas*/, o2::aod::FMcParticles const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    colHistManager.fill<modes::Mode::kAnalysis_Qa_Mc>(col, mcCols);
    auto omegaSlice = omegaWithLabelPartition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    for (auto const& omega : omegaSlice) {
      if (!omegaCleaner.isClean(omega, mcParticles, mcMothers, mcPartonicMothers)) {
        continue;
      }
      omegaHistManager.fill<modes::Mode::kAnalysis_Qa_Mc>(omega, tracks, mcParticles, mcMothers, mcPartonicMothers);
    }
  }
  PROCESS_SWITCH(FemtoCascadeQa, processOmegaMc, "Process Omegas with MC information", false);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  o2::framework::WorkflowSpec workflow{
    adaptAnalysisTask<FemtoCascadeQa>(cfgc),
  };
  return workflow;
}
