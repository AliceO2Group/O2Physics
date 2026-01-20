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

/// \file femtoKinkQa.cxx
/// \brief QA task for kinks
/// \author Anton Riedel, TU München, anton.riedel@cern.ch
/// \author Henrik Fribert, TU München, henrik.fribert@cern.ch

#include "PWGCF/Femto/Core/collisionBuilder.h"
#include "PWGCF/Femto/Core/collisionHistManager.h"
#include "PWGCF/Femto/Core/kinkBuilder.h"
#include "PWGCF/Femto/Core/kinkHistManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/particleCleaner.h"
#include "PWGCF/Femto/Core/partitions.h"
#include "PWGCF/Femto/Core/trackHistManager.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"
#include "PWGLF/DataModel/LFKinkDecayTables.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/Expressions.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/InitContext.h"
#include "Framework/OutputObjHeader.h"
#include "Framework/runDataProcessing.h"

#include <string>
#include <vector>

using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femto;

struct FemtoKinkQa {

  // setup tables
  using FemtoCollisions = o2::soa::Join<FCols, FColMasks, FColPos, FColSphericities, FColMults>;
  using FilteredFemtoCollisions = o2::soa::Filtered<FemtoCollisions>;
  using FilteredFemtoCollision = FilteredFemtoCollisions::iterator;

  using FemtoCollisionsWithLabel = o2::soa::Join<FemtoCollisions, FColLabels>;
  using FilteredFemtoCollisionsWithLabel = o2::soa::Filtered<FemtoCollisionsWithLabel>;
  using FilteredFemtoCollisionWithLabel = FilteredFemtoCollisionsWithLabel::iterator;

  // Define kink/sigma tables (joining tables for comprehensive information)
  using FemtoSigmas = o2::soa::Join<FSigmas, FSigmaMasks, FSigmaExtras>;
  using FemtoSigmaPlus = o2::soa::Join<FSigmaPlus, FSigmaPlusMasks, FSigmaPlusExtras>;
  using FemtoTracks = o2::soa::Join<FTracks, FTrackDcas, FTrackExtras, FTrackPids>;

  using FemtoSigmasWithLabel = o2::soa::Join<FemtoSigmas, FSigmaLabels>;
  using FemtoSigmaPlusWithLabel = o2::soa::Join<FemtoSigmaPlus, FK0shortLabels>;
  using FemtoTracksWithLabel = o2::soa::Join<FemtoTracks, FTrackLabels>;

  SliceCache cache;

  // setup for collisions
  collisionbuilder::ConfCollisionSelection collisionSelection;
  Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);

  colhistmanager::CollisionHistManager colHistManager;
  colhistmanager::ConfCollisionBinning confCollisionBinning;
  colhistmanager::ConfCollisionQaBinning confCollisionQaBinning;

  // setup for sigmas
  kinkbuilder::ConfSigmaSelection1 confSigmaSelection;

  particlecleaner::ConfSigmaCleaner1 confSigmaCleaner;
  particlecleaner::ParticleCleaner sigmaCleaner;

  Partition<FemtoSigmas> sigmaPartition = MAKE_SIGMA_PARTITION(confSigmaSelection);
  Preslice<FemtoSigmas> perColSigmas = femtobase::stored::fColId;

  Partition<FemtoSigmasWithLabel> sigmaWithLabelPartition = MAKE_SIGMA_PARTITION(confSigmaSelection);
  Preslice<FemtoSigmasWithLabel> perColSigmasWithLabel = femtobase::stored::fColId;

  kinkhistmanager::ConfSigmaBinning1 confSigmaBinning;
  kinkhistmanager::ConfSigmaQaBinning1 confSigmaQaBinning;
  kinkhistmanager::KinkHistManager<
    kinkhistmanager::PrefixSigmaQa,
    trackhistmanager::PrefixKinkChaDaughterQa,
    modes::Kink::kSigma>
    sigmaHistManager;

  // setup for sigma plus
  kinkbuilder::ConfSigmaPlusSelection1 confSigmaPlusSelection;

  particlecleaner::ConfSigmaPlusCleaner1 confSigmaPlusCleaner;
  particlecleaner::ParticleCleaner sigmaPlusCleaner;

  Partition<FemtoSigmaPlus> sigmaPlusPartition = MAKE_SIGMAPLUS_PARTITION(confSigmaPlusSelection);
  Preslice<FemtoSigmaPlus> perColSigmaPlus = femtobase::stored::fColId;

  Partition<FemtoSigmaPlusWithLabel> sigmaPlusWithLabelPartition = MAKE_SIGMAPLUS_PARTITION(confSigmaSelection);
  Preslice<FemtoSigmaPlusWithLabel> perColSigmaPlussWithLabel = femtobase::stored::fColId;

  kinkhistmanager::ConfSigmaPlusBinning1 confSigmaPlusBinning;
  kinkhistmanager::ConfSigmaPlusQaBinning1 confSigmaPlusQaBinning;
  kinkhistmanager::KinkHistManager<
    kinkhistmanager::PrefixSigmaPlusQa,
    trackhistmanager::PrefixKinkChaDaughterQa,
    modes::Kink::kSigmaPlus>
    sigmaPlusHistManager;

  // setup for daughters
  trackhistmanager::ConfKinkChaDauBinning confKinkChaDaughterBinning;
  trackhistmanager::ConfKinkChaDauQaBinning confKinkChaDaughterQaBinning;

  HistogramRegistry hRegistry{"FemtoKinkQa", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    if ((doprocessSigma + doprocessSigmaMc + doprocessSigmaPlus + doprocessSigmaPlusMc) > 1) {
      LOG(fatal) << "Only one process can be activated";
    }

    bool processData = doprocessSigma || doprocessSigmaPlus;
    sigmaCleaner.init(confSigmaCleaner);
    sigmaPlusCleaner.init(confSigmaPlusCleaner);
    if (processData) {
      auto colHistSpec = colhistmanager::makeColQaHistSpecMap(confCollisionBinning, confCollisionQaBinning);
      colHistManager.init<modes::Mode::kAnalysis_Qa>(&hRegistry, colHistSpec, confCollisionQaBinning);
      auto chaDauHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confKinkChaDaughterBinning, confKinkChaDaughterQaBinning);
      if (doprocessSigma) {
        auto sigmaHistSpec = kinkhistmanager::makeKinkQaHistSpecMap(confSigmaBinning, confSigmaQaBinning);
        sigmaHistManager.init<modes::Mode::kAnalysis_Qa>(&hRegistry, sigmaHistSpec, confSigmaSelection, confSigmaQaBinning, chaDauHistSpec, confKinkChaDaughterQaBinning);
      }
      if (doprocessSigmaPlus) {
        auto sigmaPlusHistSpec = kinkhistmanager::makeKinkQaHistSpecMap(confSigmaPlusBinning, confSigmaPlusQaBinning);
        sigmaPlusHistManager.init<modes::Mode::kAnalysis_Qa>(&hRegistry, sigmaPlusHistSpec, confSigmaPlusSelection, confSigmaPlusQaBinning, chaDauHistSpec, confKinkChaDaughterQaBinning);
      }
    } else {
      auto colHistSpec = colhistmanager::makeColMcQaHistSpecMap(confCollisionBinning, confCollisionQaBinning);
      colHistManager.init<modes::Mode::kAnalysis_Qa_Mc>(&hRegistry, colHistSpec, confCollisionQaBinning);
      auto chaDauHistSpec = trackhistmanager::makeTrackMcQaHistSpecMap(confKinkChaDaughterBinning, confKinkChaDaughterQaBinning);
      if (doprocessSigmaMc) {
        auto sigmaHistSpec = kinkhistmanager::makeKinkMcQaHistSpecMap(confSigmaBinning, confSigmaQaBinning);
        sigmaHistManager.init<modes::Mode::kAnalysis_Qa_Mc>(&hRegistry, sigmaHistSpec, confSigmaSelection, confSigmaQaBinning, chaDauHistSpec, confKinkChaDaughterQaBinning);
      }
      if (doprocessSigmaPlusMc) {
        auto sigmaPlusHistSpec = kinkhistmanager::makeKinkMcQaHistSpecMap(confSigmaPlusBinning, confSigmaPlusQaBinning);
        sigmaPlusHistManager.init<modes::Mode::kAnalysis_Qa_Mc>(&hRegistry, sigmaPlusHistSpec, confSigmaPlusSelection, confSigmaPlusQaBinning, chaDauHistSpec, confKinkChaDaughterQaBinning);
      }
    }
    hRegistry.print();
  };

  void processSigma(FilteredFemtoCollision const& col, FemtoSigmas const& /*sigmas*/, FemtoTracks const& tracks)
  {
    colHistManager.fill<modes::Mode::kAnalysis_Qa>(col);
    auto sigmaSlice = sigmaPartition->sliceByCached(femtobase::stored::fColId, col.globalIndex(), cache);
    for (auto const& sigma : sigmaSlice) {
      sigmaHistManager.fill<modes::Mode::kAnalysis_Qa>(sigma, tracks);
    }
  }
  PROCESS_SWITCH(FemtoKinkQa, processSigma, "Process sigmas", true);

  void processSigmaMc(FilteredFemtoCollisionWithLabel const& col, FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoSigmasWithLabel const& /*sigmas*/, FMcParticles const& mcParticles, FMcMothers const& mcMothers, FMcPartMoths const& mcPartonicMothers)
  {
    colHistManager.fill<modes::Mode::kAnalysis_Qa_Mc>(col, mcCols);
    auto sigmaSlice = sigmaWithLabelPartition->sliceByCached(femtobase::stored::fColId, col.globalIndex(), cache);
    for (auto const& sigma : sigmaSlice) {
      if (!sigmaCleaner.isClean(sigma, mcParticles, mcMothers, mcPartonicMothers)) {
        continue;
      }
      sigmaHistManager.fill<modes::Mode::kAnalysis_Qa_Mc>(sigma, tracks, mcParticles, mcMothers, mcPartonicMothers);
    }
  }
  PROCESS_SWITCH(FemtoKinkQa, processSigmaMc, "Process sigmas", false);

  void processSigmaPlus(FilteredFemtoCollision const& col, FemtoSigmaPlus const& /*sigmaplus*/, FemtoTracks const& tracks)
  {
    colHistManager.fill<modes::Mode::kAnalysis_Qa>(col);

    auto sigmaplusSlice = sigmaPlusPartition->sliceByCached(femtobase::stored::fColId, col.globalIndex(), cache);

    for (auto const& sp : sigmaplusSlice) {
      sigmaPlusHistManager.fill<modes::Mode::kAnalysis_Qa>(sp, tracks);
    }
  }
  PROCESS_SWITCH(FemtoKinkQa, processSigmaPlus, "Process sigma plus", false);

  void processSigmaPlusMc(FilteredFemtoCollisionWithLabel const& col, FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoSigmaPlusWithLabel const& /*sigmaPlus*/, FMcParticles const& mcParticles, FMcMothers const& mcMothers, FMcPartMoths const& mcPartonicMothers)
  {
    colHistManager.fill<modes::Mode::kAnalysis_Qa_Mc>(col, mcCols);
    auto sigmaPlusSlice = sigmaPlusWithLabelPartition->sliceByCached(femtobase::stored::fColId, col.globalIndex(), cache);
    for (auto const& sigmaPlus : sigmaPlusSlice) {
      if (!sigmaPlusCleaner.isClean(sigmaPlus, mcParticles, mcMothers, mcPartonicMothers)) {
        continue;
      }
      sigmaPlusHistManager.fill<modes::Mode::kAnalysis_Qa_Mc>(sigmaPlus, tracks, mcParticles, mcMothers, mcPartonicMothers);
    }
  }
  PROCESS_SWITCH(FemtoKinkQa, processSigmaPlusMc, "Process sigmas", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtoKinkQa>(cfgc),
  };
  return workflow;
}
