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

#include <map>
#include <vector>

using namespace o2::analysis::femto;

struct FemtoKinkQa {

  // setup tables
  using FemtoCollisions = o2::soa::Join<o2::aod::FCols, o2::aod::FColMasks, o2::aod::FColPos, o2::aod::FColSphericities, o2::aod::FColMults>;
  using FilteredFemtoCollisions = o2::soa::Filtered<FemtoCollisions>;
  using FilteredFemtoCollision = FilteredFemtoCollisions::iterator;

  using FemtoCollisionsWithLabel = o2::soa::Join<FemtoCollisions, o2::aod::FColLabels>;
  using FilteredFemtoCollisionsWithLabel = o2::soa::Filtered<FemtoCollisionsWithLabel>;
  using FilteredFemtoCollisionWithLabel = FilteredFemtoCollisionsWithLabel::iterator;

  // Define kink/sigma tables (joining tables for comprehensive information)
  using FemtoSigmas = o2::soa::Join<o2::aod::FSigmas, o2::aod::FSigmaMasks, o2::aod::FSigmaExtras>;
  using FemtoSigmaPlus = o2::soa::Join<o2::aod::FSigmaPlus, o2::aod::FSigmaPlusMasks, o2::aod::FSigmaPlusExtras>;
  using FemtoTracks = o2::soa::Join<o2::aod::FTracks, o2::aod::FTrackDcas, o2::aod::FTrackExtras, o2::aod::FTrackPids>;

  using FemtoSigmasWithLabel = o2::soa::Join<FemtoSigmas, o2::aod::FSigmaLabels>;
  using FemtoSigmaPlusWithLabel = o2::soa::Join<FemtoSigmaPlus, o2::aod::FSigmaPlusLabels>;
  using FemtoTracksWithLabel = o2::soa::Join<FemtoTracks, o2::aod::FTrackLabels>;

  o2::framework::SliceCache cache;

  // setup for collisions
  collisionbuilder::ConfCollisionSelection collisionSelection;
  o2::framework::expressions::Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);

  colhistmanager::CollisionHistManager colHistManager;
  colhistmanager::ConfCollisionBinning confCollisionBinning;
  colhistmanager::ConfCollisionQaBinning confCollisionQaBinning;

  // setup for sigmas
  kinkbuilder::ConfSigmaSelection1 confSigmaSelection;

  particlecleaner::ConfSigmaCleaner1 confSigmaCleaner;
  particlecleaner::ParticleCleaner sigmaCleaner;

  o2::framework::Partition<FemtoSigmas> sigmaPartition = MAKE_SIGMA_PARTITION(confSigmaSelection);
  o2::framework::Preslice<FemtoSigmas> perColSigmas = o2::aod::femtobase::stored::fColId;

  o2::framework::Partition<FemtoSigmasWithLabel> sigmaWithLabelPartition = MAKE_SIGMA_PARTITION(confSigmaSelection);
  o2::framework::Preslice<FemtoSigmasWithLabel> perColSigmasWithLabel = o2::aod::femtobase::stored::fColId;

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

  o2::framework::Partition<FemtoSigmaPlus> sigmaPlusPartition = MAKE_SIGMAPLUS_PARTITION(confSigmaPlusSelection);
  o2::framework::Preslice<FemtoSigmaPlus> perColSigmaPlus = o2::aod::femtobase::stored::fColId;

  o2::framework::Partition<FemtoSigmaPlusWithLabel> sigmaPlusWithLabelPartition = MAKE_SIGMAPLUS_PARTITION(confSigmaSelection);
  o2::framework::Preslice<FemtoSigmaPlusWithLabel> perColSigmaPlussWithLabel = o2::aod::femtobase::stored::fColId;

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

  o2::framework::HistogramRegistry hRegistry{"FemtoKinkQa", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    if ((doprocessSigma + doprocessSigmaMc + doprocessSigmaPlus + doprocessSigmaPlusMc) > 1) {
      LOG(fatal) << "Only one process can be activated";
    }

    bool processData = doprocessSigma || doprocessSigmaPlus;

    sigmaCleaner.init(confSigmaCleaner);
    sigmaPlusCleaner.init(confSigmaPlusCleaner);

    std::map<colhistmanager::ColHist, std::vector<o2::framework::AxisSpec>> colHistSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> chaDauHistSpec;
    std::map<kinkhistmanager::KinkHist, std::vector<o2::framework::AxisSpec>> sigmaHistSpec;
    std::map<kinkhistmanager::KinkHist, std::vector<o2::framework::AxisSpec>> sigmaPlusHistSpec;

    if (processData) {
      colHistSpec = colhistmanager::makeColQaHistSpecMap(confCollisionBinning, confCollisionQaBinning);
      colHistManager.init<modes::Mode::kAnalysis_Qa>(&hRegistry, colHistSpec, confCollisionQaBinning);
      chaDauHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confKinkChaDaughterBinning, confKinkChaDaughterQaBinning);
      if (doprocessSigma) {
        sigmaHistSpec = kinkhistmanager::makeKinkQaHistSpecMap(confSigmaBinning, confSigmaQaBinning);
        sigmaHistManager.init<modes::Mode::kAnalysis_Qa>(&hRegistry, sigmaHistSpec, confSigmaSelection, confSigmaQaBinning, chaDauHistSpec, confKinkChaDaughterQaBinning);
      }
      if (doprocessSigmaPlus) {
        sigmaPlusHistSpec = kinkhistmanager::makeKinkQaHistSpecMap(confSigmaPlusBinning, confSigmaPlusQaBinning);
        sigmaPlusHistManager.init<modes::Mode::kAnalysis_Qa>(&hRegistry, sigmaPlusHistSpec, confSigmaPlusSelection, confSigmaPlusQaBinning, chaDauHistSpec, confKinkChaDaughterQaBinning);
      }
    } else {
      colHistSpec = colhistmanager::makeColMcQaHistSpecMap(confCollisionBinning, confCollisionQaBinning);
      colHistManager.init<modes::Mode::kAnalysis_Qa_Mc>(&hRegistry, colHistSpec, confCollisionQaBinning);
      chaDauHistSpec = trackhistmanager::makeTrackMcQaHistSpecMap(confKinkChaDaughterBinning, confKinkChaDaughterQaBinning);
      if (doprocessSigmaMc) {
        sigmaHistSpec = kinkhistmanager::makeKinkMcQaHistSpecMap(confSigmaBinning, confSigmaQaBinning);
        sigmaHistManager.init<modes::Mode::kAnalysis_Qa_Mc>(&hRegistry, sigmaHistSpec, confSigmaSelection, confSigmaQaBinning, chaDauHistSpec, confKinkChaDaughterQaBinning);
      }
      if (doprocessSigmaPlusMc) {
        sigmaPlusHistSpec = kinkhistmanager::makeKinkMcQaHistSpecMap(confSigmaPlusBinning, confSigmaPlusQaBinning);
        sigmaPlusHistManager.init<modes::Mode::kAnalysis_Qa_Mc>(&hRegistry, sigmaPlusHistSpec, confSigmaPlusSelection, confSigmaPlusQaBinning, chaDauHistSpec, confKinkChaDaughterQaBinning);
      }
    }
    hRegistry.print();
  };

  void processSigma(FilteredFemtoCollision const& col, FemtoSigmas const& /*sigmas*/, FemtoTracks const& tracks)
  {
    colHistManager.fill<modes::Mode::kAnalysis_Qa>(col);
    auto sigmaSlice = sigmaPartition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    for (auto const& sigma : sigmaSlice) {
      sigmaHistManager.fill<modes::Mode::kAnalysis_Qa>(sigma, tracks);
    }
  }
  PROCESS_SWITCH(FemtoKinkQa, processSigma, "Process sigmas", true);

  void processSigmaMc(FilteredFemtoCollisionWithLabel const& col, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoSigmasWithLabel const& /*sigmas*/, o2::aod::FMcParticles const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    colHistManager.fill<modes::Mode::kAnalysis_Qa_Mc>(col, mcCols);
    auto sigmaSlice = sigmaWithLabelPartition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
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

    auto sigmaplusSlice = sigmaPlusPartition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);

    for (auto const& sp : sigmaplusSlice) {
      sigmaPlusHistManager.fill<modes::Mode::kAnalysis_Qa>(sp, tracks);
    }
  }
  PROCESS_SWITCH(FemtoKinkQa, processSigmaPlus, "Process sigma plus", false);

  void processSigmaPlusMc(FilteredFemtoCollisionWithLabel const& col, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoSigmaPlusWithLabel const& /*sigmaPlus*/, o2::aod::FMcParticles const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    colHistManager.fill<modes::Mode::kAnalysis_Qa_Mc>(col, mcCols);
    auto sigmaPlusSlice = sigmaPlusWithLabelPartition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    for (auto const& sigmaPlus : sigmaPlusSlice) {
      if (!sigmaPlusCleaner.isClean(sigmaPlus, mcParticles, mcMothers, mcPartonicMothers)) {
        continue;
      }
      sigmaPlusHistManager.fill<modes::Mode::kAnalysis_Qa_Mc>(sigmaPlus, tracks, mcParticles, mcMothers, mcPartonicMothers);
    }
  }
  PROCESS_SWITCH(FemtoKinkQa, processSigmaPlusMc, "Process sigmas", false);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  o2::framework::WorkflowSpec workflow{
    adaptAnalysisTask<FemtoKinkQa>(cfgc),
  };
  return workflow;
}
