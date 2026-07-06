// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file treeCreatorPidTpcQa.cxx
/// \brief Task to produce table with clean selections for TPC PID calibration
///
/// \author Ana Marin <ana.marin@cern.ch>
/// \author Oleksii Lubynets <oleksii.lubynets@cern.ch>

#include "treeCreatorPidTpcQa.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "DPG/Tasks/TPC/tpcSkimsTableCreator.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/StaticFor.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/PID.h>

#include <cmath>

using namespace o2;
using namespace o2::framework;
using namespace o2::track;
using namespace o2::dpg_pidtpcqa;

struct treeCreatorPidTpcQa {
  Produces<o2::aod::QaPidTpc> rowPidTpcQa;

  Configurable<int> applyEvSel{"applyEvSel", 2, "Flag to apply event selection: 0 -> no event selection, 1 -> Run 2 event selection, 2 -> Run 3 event selection"};
  Configurable<float> cutVtxZ{"cutVtxZ", 10.f, "Cut on vertex Z position [cm]. Set negative value to switch this cut off"};
  Configurable<int> trackSelection{"trackSelection", 1, "Track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> kQualityTracks, 5 -> kInAcceptanceTracks"};
  Configurable<bool> requireGlobalTrack{"requireGlobalTrack", true, "Skip non-global tracks"};
  Configurable<bool> requireIts{"requireIts", true, "Skip tracks without ITS"};
  Configurable<int16_t> cutMinTPCNcls{"cutMinTPCNcls", 0, "Minimum number or TPC Clusters for tracks"};
  Configurable<float> cutRapidity{"cutRapidity", 0.5, "Rapidity cut. Set negative value to switch this cut off"};
  Configurable<float> nClNorm{"nClNorm", 152., "Number of cluster normalization. Run 2: 159, Run 3 152"};

  using CollisionsExtra = soa::Join<aod::Collisions, aod::Mults, aod::EvSels>;
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection>;

  Preslice<TrackCandidates> perCollisionTracks = aod::track::collisionId;

  int mEnabledParticles{0};

  // Track selection
  template <typename TrackType>
  bool isTrackSelected(const TrackType& track)
  {
    bool isSelected{false};
    isSelected |= trackSelection == TrackSelectionNoCut;
    isSelected |= (trackSelection == TrackSelectionGlobalTrack) && track.isGlobalTrack();
    isSelected |= (trackSelection == TrackSelectionTrackWoPtEta) && track.isGlobalTrackWoPtEta();
    isSelected |= (trackSelection == TrackSelectionGlobalTrackWoDCA) && track.isGlobalTrackWoDCA();
    isSelected |= (trackSelection == TrackSelectionQualityTracks) && track.isQualityTrack();
    isSelected |= (trackSelection == TrackSelectionInAcceptanceTracks) && track.isInAcceptanceTrack();
    isSelected &= (!requireGlobalTrack || track.isGlobalTrack());
    isSelected &= (!requireIts || track.hasITS());
    isSelected &= track.hasTPC();
    isSelected &= (track.tpcNClsFound() >= cutMinTPCNcls);

    return isSelected;
  }

  template <o2::track::PID::ID Id>
  bool initPerParticle()
  {
    static_assert(Id >= 0 && Id <= PID::Alpha && "Particle index outside limits");
    int enabledProcesses{0};

    switch (Id) {
#define PARTICLE_CASE(ParticleId)                                                                    \
  case PID::ParticleId:                                                                              \
    if (!doprocess##ParticleId && !doprocessFull##ParticleId && !doprocessFullWithTOF##ParticleId) { \
      return false;                                                                                  \
    }                                                                                                \
    if (doprocess##ParticleId) {                                                                     \
      ++enabledProcesses;                                                                            \
    }                                                                                                \
    if (doprocessFull##ParticleId) {                                                                 \
      ++enabledProcesses;                                                                            \
    }                                                                                                \
    if (doprocessFullWithTOF##ParticleId) {                                                          \
      ++enabledProcesses;                                                                            \
    }                                                                                                \
    LOG(info) << "Enabled TPC QA for " << #ParticleId;                                               \
    break;

      PARTICLE_CASE(Electron);
      PARTICLE_CASE(Muon);
      PARTICLE_CASE(Pion);
      PARTICLE_CASE(Kaon);
      PARTICLE_CASE(Proton);
      PARTICLE_CASE(Deuteron);
      PARTICLE_CASE(Triton);
      PARTICLE_CASE(Helium3);
      PARTICLE_CASE(Alpha);
#undef PARTICLE_CASE
    }
    if (enabledProcesses != 1) {
      LOG(fatal) << "Cannot enable more than one process function per particle, check and retry!";
    }
    return true;
  }

  void init(o2::framework::InitContext&)
  {
    static_for<0, PID::Alpha>([&](auto Id) {
      mEnabledParticles += static_cast<int>(initPerParticle<Id>());
    });
  }

  template <o2::track::PID::ID Id, typename TrackType>
  void processSingleParticle(CollisionsExtra const& collisions,
                             TrackType const& tracks)
  {
    rowPidTpcQa.reserve(tracks.size() * mEnabledParticles);

    for (const auto& collision : collisions) {
      if (!isEventSelected(collision, applyEvSel) || ((cutVtxZ > 0.f) && std::abs(collision.posZ()) > cutVtxZ)) {
        continue;
      }

      const float ft0Occ = collision.ft0cOccupancyInTimeRange();
      const float multTPC = collision.multTPC() / dpg_tpcskimstablecreator::MultiplicityNorm;

      const auto tracksFromCollision = tracks.sliceBy(perCollisionTracks, static_cast<int>(collision.globalIndex()));

      for (const auto& track : tracksFromCollision) {
        if (!isTrackSelected(track)) {
          continue;
        }
        const float rapidity = track.rapidity(PID::getMass(Id));
        if (cutRapidity > 0.f && std::fabs(rapidity) > cutRapidity) {
          continue;
        }

        const int nClNormalized = std::sqrt(nClNorm / track.tpcNClsFound());
        const float phi = track.phi();
        const float tgl = track.tgl();
        const float tpcInnerParam = track.tpcInnerParam();

        rowPidTpcQa(Id, ft0Occ, multTPC, nClNormalized, phi, tgl, tpcInnerParam, rapidity);
      } // tracks
    } // collisions
  }

  // QA of nsigma only tables
#define MAKE_PROCESS_FUNCTION(PidTableTPC, ParticleId)                            \
  void process##ParticleId(CollisionsExtra const& collisions,                     \
                           soa::Join<TrackCandidates, PidTableTPC> const& tracks) \
  {                                                                               \
    processSingleParticle<PID::ParticleId>(collisions, tracks);                   \
  }                                                                               \
  PROCESS_SWITCH(treeCreatorPidTpcQa, process##ParticleId, Form("Process for the %s hypothesis for TPC NSigma QA", #ParticleId), false);

  MAKE_PROCESS_FUNCTION(aod::pidTPCEl, Electron);
  MAKE_PROCESS_FUNCTION(aod::pidTPCMu, Muon);
  MAKE_PROCESS_FUNCTION(aod::pidTPCPi, Pion);
  MAKE_PROCESS_FUNCTION(aod::pidTPCKa, Kaon);
  MAKE_PROCESS_FUNCTION(aod::pidTPCPr, Proton);
  MAKE_PROCESS_FUNCTION(aod::pidTPCDe, Deuteron);
  MAKE_PROCESS_FUNCTION(aod::pidTPCTr, Triton);
  MAKE_PROCESS_FUNCTION(aod::pidTPCHe, Helium3);
  MAKE_PROCESS_FUNCTION(aod::pidTPCAl, Alpha);
#undef MAKE_PROCESS_FUNCTION

// QA of full tables
#define MAKE_PROCESS_FUNCTION(PidTableTPC, ParticleId)                                \
  void processFull##ParticleId(CollisionsExtra const& collisions,                     \
                               soa::Join<TrackCandidates, PidTableTPC> const& tracks) \
  {                                                                                   \
    processSingleParticle<PID::ParticleId>(collisions, tracks);                       \
  }                                                                                   \
  PROCESS_SWITCH(treeCreatorPidTpcQa, processFull##ParticleId, Form("Process for the %s hypothesis for full TPC PID QA", #ParticleId), false);

  MAKE_PROCESS_FUNCTION(aod::pidTPCFullEl, Electron);
  MAKE_PROCESS_FUNCTION(aod::pidTPCFullMu, Muon);
  MAKE_PROCESS_FUNCTION(aod::pidTPCFullPi, Pion);
  MAKE_PROCESS_FUNCTION(aod::pidTPCFullKa, Kaon);
  MAKE_PROCESS_FUNCTION(aod::pidTPCFullPr, Proton);
  MAKE_PROCESS_FUNCTION(aod::pidTPCFullDe, Deuteron);
  MAKE_PROCESS_FUNCTION(aod::pidTPCFullTr, Triton);
  MAKE_PROCESS_FUNCTION(aod::pidTPCFullHe, Helium3);
  MAKE_PROCESS_FUNCTION(aod::pidTPCFullAl, Alpha);
#undef MAKE_PROCESS_FUNCTION

  // QA of full tables with TOF information
#define MAKE_PROCESS_FUNCTION(PidTableTPC, PidTableTOF, ParticleId)                                       \
  void processFullWithTOF##ParticleId(CollisionsExtra const& collisions,                                  \
                                      soa::Join<TrackCandidates, PidTableTPC, PidTableTOF> const& tracks) \
  {                                                                                                       \
    processSingleParticle<PID::ParticleId>(collisions, tracks);                                           \
  }                                                                                                       \
  PROCESS_SWITCH(treeCreatorPidTpcQa, processFullWithTOF##ParticleId, Form("Process for the %s hypothesis for full TPC PID QA with the TOF info added", #ParticleId), false);

  MAKE_PROCESS_FUNCTION(aod::pidTPCFullEl, aod::pidTOFFullEl, Electron);
  MAKE_PROCESS_FUNCTION(aod::pidTPCFullMu, aod::pidTOFFullMu, Muon);
  MAKE_PROCESS_FUNCTION(aod::pidTPCFullPi, aod::pidTOFFullPi, Pion);
  MAKE_PROCESS_FUNCTION(aod::pidTPCFullKa, aod::pidTOFFullKa, Kaon);
  MAKE_PROCESS_FUNCTION(aod::pidTPCFullPr, aod::pidTOFFullPr, Proton);
  MAKE_PROCESS_FUNCTION(aod::pidTPCFullDe, aod::pidTOFFullDe, Deuteron);
  MAKE_PROCESS_FUNCTION(aod::pidTPCFullTr, aod::pidTOFFullTr, Triton);
  MAKE_PROCESS_FUNCTION(aod::pidTPCFullHe, aod::pidTOFFullHe, Helium3);
  MAKE_PROCESS_FUNCTION(aod::pidTPCFullAl, aod::pidTOFFullAl, Alpha);
#undef MAKE_PROCESS_FUNCTION
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<treeCreatorPidTpcQa>(cfgc)};
}
