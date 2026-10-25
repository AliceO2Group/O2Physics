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

/// \file femtoTrackTrackTwoPhotonResonanceQa.cxx
/// \brief QA task for resonances built from two tracks and a two-photon resonance (omega(782) -> pi+ pi- pi0)
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#include "PWGCF/Femto/Core/collisionBuilder.h"
#include "PWGCF/Femto/Core/collisionHistManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/partitions.h"
#include "PWGCF/Femto/Core/photonHistManager.h"
#include "PWGCF/Femto/Core/trackHistManager.h"
#include "PWGCF/Femto/Core/trackTrackTwoPhotonResonanceBuilder.h"
#include "PWGCF/Femto/Core/trackTrackTwoPhotonResonanceHistManager.h"
#include "PWGCF/Femto/Core/twoPhotonResonanceHistManager.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

using namespace o2::analysis::femto;

struct FemtoTrackTrackTwoPhotonResonanceQa {

  // setup tables
  using FemtoCollisions = o2::soa::Join<o2::aod::FCols, o2::aod::FColMasks, o2::aod::FColPos, o2::aod::FColSphericities, o2::aod::FColMults, o2::aod::FColCents>;
  using FilteredFemtoCollisions = o2::soa::Filtered<FemtoCollisions>;
  using FilteredFemtoCollision = FilteredFemtoCollisions::iterator;

  using FemtoOmegaMesons = o2::soa::Join<o2::aod::FOmegaMesons, o2::aod::FOmegaMesonMasks>;
  using FemtoTracks = o2::soa::Join<o2::aod::FTracks, o2::aod::FTrackMasks, o2::aod::FTrackMass, o2::aod::FTrackDcas, o2::aod::FTrackExtras, o2::aod::FTrackPids>;
  using FemtoPi0s = o2::soa::Join<o2::aod::FPi0s, o2::aod::FPi0Masks>;
  using FemtoPhotons = o2::soa::Join<o2::aod::FPhotons, o2::aod::FPhotonMasks, o2::aod::FPhotonExtras>;

  o2::framework::SliceCache cache;

  // setup for collisions
  collisionbuilder::ConfCollisionSelection collisionSelection;
  o2::framework::expressions::Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);
  colhistmanager::CollisionHistManager colHistManager;
  colhistmanager::ConfCollisionBinning confCollisionBinning;
  colhistmanager::ConfCollisionQaBinning confCollisionQaBinning;

  // setup for omega mesons
  trackstwophotonresonancebuilder::ConfOmegaMesonSelection1 confOmegaMesonSelection;
  o2::framework::Partition<FemtoOmegaMesons> omegaMesonPartition = MAKE_TRACKTRACKTWOPHOTONRESONANCE_PARTITION(confOmegaMesonSelection);
  o2::framework::Preslice<FemtoOmegaMesons> perColOmegaMesons = o2::aod::femtobase::stored::fColId;

  trackstwophotonresonancehistmanager::ConfOmegaMesonBinning confOmegaMesonBinning;

  // setup for the two charged (pion) daughters
  trackhistmanager::ConfResonancePosDauBinning confPosDauBinning;
  trackhistmanager::ConfResonancePosDauQaBinning confPosDauQaBinning;
  trackhistmanager::ConfResonanceNegDauBinning confNegDauBinning;
  trackhistmanager::ConfResonanceNegDauQaBinning confNegDauQaBinning;

  // setup for the neutral (pi0) daughter -- reuses the regular Pi0 binning/QA config, since the
  // embedded candidate gets the same TwoPhotonResonanceHistManager treatment (including its own
  // photon-leg QA) as the dedicated Pi0/Eta QA task
  twophotonresonancehistmanager::ConfPi0Binning confNeutralDauBinning;
  photonhistmanager::ConfPhotonBinning confNeutralDauPhotonBinning;
  photonhistmanager::ConfPhotonQaBinning confNeutralDauPhotonQaBinning;

  trackstwophotonresonancehistmanager::TrackTrackTwoPhotonResonanceHistManager<
    trackstwophotonresonancehistmanager::PrefixOmegaMeson,
    trackhistmanager::PrefixResonancePosDaughterSe,
    trackhistmanager::PrefixResonanceNegDaughterSe,
    trackstwophotonresonancehistmanager::PrefixOmegaMesonNeutralDau,
    trackstwophotonresonancehistmanager::PrefixOmegaMesonNeutralDauPhotonDau1,
    trackstwophotonresonancehistmanager::PrefixOmegaMesonNeutralDauPhotonDau2,
    modes::TrackTrackTwoPhotonResonance::kOmegaMeson,
    modes::TwoPhotonResonance::kPi0>
    omegaMesonHistManager;

  o2::framework::HistogramRegistry hRegistry{"FemtoTrackTrackTwoPhotonResonanceQa", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    auto colHistSpec = colhistmanager::makeColQaHistSpecMap(confCollisionBinning, confCollisionQaBinning);
    colHistManager.init<modes::Mode::kReco_Qa>(&hRegistry, colHistSpec, confCollisionBinning, confCollisionQaBinning);

    auto omegaMesonHistSpec = trackstwophotonresonancehistmanager::makeTrackTrackTwoPhotonResonanceQaHistSpecMap(confOmegaMesonBinning);
    auto posDauSpec = trackhistmanager::makeTrackQaHistSpecMap(confPosDauBinning, confPosDauQaBinning);
    auto negDauSpec = trackhistmanager::makeTrackQaHistSpecMap(confNegDauBinning, confNegDauQaBinning);
    auto neutralDauSpec = twophotonresonancehistmanager::makeTwoPhotonResonanceQaHistSpecMap(confNeutralDauBinning);
    auto neutralDauPhotonSpec = photonhistmanager::makePhotonQaHistSpecMap(confNeutralDauPhotonBinning, confNeutralDauPhotonQaBinning);

    omegaMesonHistManager.init<modes::Mode::kReco_Qa>(&hRegistry, omegaMesonHistSpec, posDauSpec, negDauSpec, neutralDauSpec, neutralDauPhotonSpec, confPosDauQaBinning, confNegDauQaBinning, confNeutralDauPhotonQaBinning);
  }

  void processOmegaMesons(FilteredFemtoCollision const& col, FemtoOmegaMesons const& /*omegaMesons*/, FemtoTracks const& tracks, FemtoPi0s const& pi0s, FemtoPhotons const& photons)
  {
    auto omegaMesonSlice = omegaMesonPartition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    if (omegaMesonSlice.size() == 0) {
      return;
    }
    colHistManager.fill<modes::Mode::kReco_Qa>(col);
    for (auto const& omegaMeson : omegaMesonSlice) {
      omegaMesonHistManager.fill<modes::Mode::kReco_Qa>(omegaMeson, tracks, pi0s, photons);
    }
  }
  PROCESS_SWITCH(FemtoTrackTrackTwoPhotonResonanceQa, processOmegaMesons, "Process OmegaMesons", true);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& context)
{
  o2::framework::WorkflowSpec workflow{
    adaptAnalysisTask<FemtoTrackTrackTwoPhotonResonanceQa>(context),
  };
  return workflow;
}
