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
//
/// \author Alexander Tiekoetter <alexander.tiekoetter@cern.ch>, Muenster
/// \brief Global Track estimator with track-to-collision association
/// \file globalTrackEstimator.cxx

#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"

#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/McCollisionExtra.h"
#include "Common/DataModel/Multiplicity.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/runDataProcessing.h>

#include <string>
#include <unordered_set>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;
using namespace o2::soa;

using MyEventsBasic = Join<ReducedEvents, ReducedEventsExtended>;
using MyBarrelTracks = Join<ReducedTracks, ReducedTracksBarrel, ReducedTracksBarrelCov, ReducedTracksBarrelPID>;

using MyEventsDirect = Join<Collisions, EvSels>;
using MyTracksDirect = Join<Tracks, TracksExtra, TracksCov, TracksDCA>;
using MyEventsMcWithMults = Join<McCollisions, McCollsExtra, MultMCExtras_001>;

constexpr static uint32_t GkEventFillMap = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended;
constexpr static uint32_t GkTrackFillMap = VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelCov;

constexpr static uint32_t GkEventDirectFillMap = VarManager::ObjTypes::Collision;
constexpr static uint32_t GkTrackDirectFillMap = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackCov | VarManager::ObjTypes::TrackDCA;

constexpr float MaxChargeTolerated = 1e-3;

struct GlobalTrackEstimator {
  Produces<ReducedEventsGlobalTrk> reducedEventsGlobalTrk;
  Produces<GlobalTrackMultiplicity> globalTrackMultiplicity;

  Produces<ReducedMCGlobalMultiplicities> reducedMCGlobalMultiplicities;
  Produces<McGlobalMultiplicities> mcGlobalMultiplicities;

  Service<o2::framework::O2DatabasePDG> pdg;

  PresliceUnsorted<ReducedMCTracks> mcTracksPerReducedEvent = reducedtrackMC::reducedMCeventId;

  Configurable<float> cfgShiftDCAxy{"cfgShiftDCAxy", 0.0105, "Constant shift of DCAxy paramterization"};
  Configurable<float> cfgPtMagnitudeDCAxy{"cfgPtMagnitudeDCAxy", 0.035, "Magnitude of the pT dependency of the DCAxy paramterization"};
  Configurable<float> cfgExponentDCAxy{"cfgExponentDCAxy", 1.1, "Exponent of the pT dependency of the DCAxy paramterization"};
  Configurable<float> cfgMaxDCAz{"cfgMaxDCAz", 2.0, "Maximum DCAz"};
  Configurable<float> cfgMinPt{"cfgMinPt", 0.1, "Minimum required pT"};
  Configurable<float> cfgMinMCPt{"cfgMinMCPt", 0.0, "Minimum required pT for MC particles"};
  Configurable<float> cfgEtaRange{"cfgEtaRange", 0.8, "Absolute pseudorapidity range"};
  Configurable<float> cfgMinTPCCluster{"cfgMinTPCCluster", 70, "Minimum number of crossed TPC rows"};
  Configurable<float> cfgRatioTPCClusterOverFindableCluster{"cfgRatioTPCClusterOverFindableCluster", 0.8, "Minimum ratio of crossed TPC rows over findable clusters"};
  Configurable<float> cfgMaxChi2TPC{"cfgMaxChi2TPC", 4.0, "Maximum chi2 per cluster in TPC"};
  Configurable<float> cfgMaxChi2ITS{"cfgMaxChi2ITS", 36.0, "Maximum chi2 per cluster in ITS"};
  Configurable<bool> cfgRequireITSibAny{"cfgRequireITSibAny", true, "Require a hit in the ITS inner barrel"};
  Configurable<bool> cfgRequirePVContributor{"cfgRequirePVContributor", false, "Require the track to be a PV contributor"};
  Configurable<bool> cfgRecalculateDCA{"cfgRecalculateDCA", true, "Recalculate DCA vector with associations"};

  o2::framework::Preslice<ReducedTracksAssoc> preslice = reducedtrack_association::reducedeventId;
  o2::framework::Preslice<TrackAssoc> presliceDirect = track_association::collisionId;

  void init(InitContext& context)
  {
    bool isSkimmedEnabled = context.mOptions.get<bool>("processGlobalTracksSkimmed") || context.mOptions.get<bool>("processMCSkimmed");
    bool isDirectEnabled = context.mOptions.get<bool>("processGlobalTracksDirect") || context.mOptions.get<bool>("processMCDirect");
    bool isDummyEnabled = context.mOptions.get<bool>("processDummy");

    if (isDummyEnabled) {
      if (isSkimmedEnabled || isDirectEnabled) {
        LOG(fatal) << "You have enabled both the dummy process and at least one normal process function.";
      }
      return;
    }

    if (!isSkimmedEnabled && !isDirectEnabled && !isDummyEnabled) {
      LOG(fatal) << "No process function activated.";
      return;
    }

    if (isSkimmedEnabled && isDirectEnabled) {
      LOG(fatal) << "Process functions for skimmed data model and direct data model can not run simultaneously.";
      return;
    }
  }

  template <uint32_t fillMap, typename T, typename C>
  bool isTrackSelected(T& track, C const& event)
  {
    VarManager::ResetValues(VarManager::kNBarrelTrackVariables, VarManager::kNBarrelTrackVariables);
    VarManager::FillTrack<fillMap>(track);

    if (cfgRecalculateDCA) {
      VarManager::FillTrackCollision<fillMap>(track, event);
    }

    if (track.pt() < cfgMinPt || std::abs(track.eta()) > cfgEtaRange) {
      return false;
    }

    if (cfgRequirePVContributor && !track.isPVContributor()) {
      return false;
    }

    if (track.tpcChi2NCl() > cfgMaxChi2TPC || track.itsChi2NCl() > cfgMaxChi2ITS) {
      return false;
    }

    const int nFindable = track.tpcNClsFindable();
    const int nCrossedRows = track.tpcNClsCrossedRows();

    if (nCrossedRows < static_cast<int>(cfgMinTPCCluster.value)) {
      return false;
    }

    if (nFindable <= 0) {
      return false;
    }

    const float crossedOverFindable = static_cast<float>(nCrossedRows) / static_cast<float>(nFindable);

    if (crossedOverFindable < cfgRatioTPCClusterOverFindableCluster.value) {
      return false;
    }

    if (cfgRequireITSibAny.value && ((track.itsClusterMap() & 0x7) == 0)) {
      return false;
    }

    if (std::abs(VarManager::fgValues[VarManager::kTrackDCAz]) > cfgMaxDCAz.value) {
      return false;
    }

    if (std::abs(VarManager::fgValues[VarManager::kTrackDCAxy]) > cfgShiftDCAxy + cfgPtMagnitudeDCAxy / std::pow(track.pt(), cfgExponentDCAxy)) {
      return false;
    }

    return true;
  }

  template <typename TParticles>
  int countMCGlobalParticles(TParticles const& particles)
  {
    int count = 0;

    for (auto const& particle : particles) {
      if (!particle.isPhysicalPrimary() || particle.pt() < cfgMinMCPt.value || std::abs(particle.eta()) > cfgEtaRange.value) {
        continue;
      }

      const auto* pdgParticle = pdg->GetParticle(particle.pdgCode());
      if (pdgParticle == nullptr || std::abs(pdgParticle->Charge()) < MaxChargeTolerated) {
        continue;
      }

      ++count;
    }

    return count;
  }

  void processGlobalTracksSkimmed(MyEventsBasic const& events, ReducedTracksAssoc const& assocs, MyBarrelTracks const& tracks)
  {
    std::unordered_set<int64_t> selectedTrackIds;
    selectedTrackIds.reserve(tracks.size());

    for (auto const& track : tracks) {
      selectedTrackIds.insert(track.globalIndex());
    }

    for (const auto& event : events) {
      VarManager::ResetValues(0, VarManager::kNBarrelTrackVariables);
      VarManager::FillEvent<GkEventFillMap>(event);

      auto groupedAssocs = assocs.sliceBy(preslice, event.globalIndex());

      int nTrkGlobal = 0;

      for (auto const& assoc : groupedAssocs) {
        if (selectedTrackIds.find(assoc.reducedtrackId()) == selectedTrackIds.end()) {
          continue;
        }

        auto track = assoc.template reducedtrack_as<MyBarrelTracks>();

        if (!isTrackSelected<GkTrackFillMap>(track, event)) {
          continue;
        }

        ++nTrkGlobal;
      }

      reducedEventsGlobalTrk(nTrkGlobal);
    }
  }

  PROCESS_SWITCH(GlobalTrackEstimator, processGlobalTracksSkimmed, "Count global tracks from skimmed data model", false);

  void processGlobalTracksDirect(MyEventsDirect const& events, TrackAssoc const& assocs, MyTracksDirect const& tracks)
  {
    std::unordered_set<int64_t> selectedTrackIds;
    selectedTrackIds.reserve(tracks.size());

    for (auto const& track : tracks) {
      selectedTrackIds.insert(track.globalIndex());
    }

    for (const auto& event : events) {
      VarManager::ResetValues(0, VarManager::kNBarrelTrackVariables);
      VarManager::FillEvent<GkEventDirectFillMap>(event);

      auto groupedAssocs = assocs.sliceBy(presliceDirect, event.globalIndex());

      int nTrkGlobal = 0;

      for (auto const& assoc : groupedAssocs) {
        if (selectedTrackIds.find(assoc.trackId()) == selectedTrackIds.end()) {
          continue;
        }

        auto track = assoc.template track_as<MyTracksDirect>();

        if (!isTrackSelected<GkTrackDirectFillMap>(track, event)) {
          continue;
        }

        ++nTrkGlobal;
      }

      globalTrackMultiplicity(nTrkGlobal);
    }
  }

  PROCESS_SWITCH(GlobalTrackEstimator, processGlobalTracksDirect, "Count global tracks from unskimmed data", false);

  void processMCSkimmed(ReducedMCEvents const& events, ReducedMCTracks const& particles)
  {
    for (auto const& event : events) {
      auto eventParticles = particles.sliceBy(mcTracksPerReducedEvent, event.globalIndex());
      reducedMCGlobalMultiplicities(countMCGlobalParticles(eventParticles));
    }
  }

  PROCESS_SWITCH(GlobalTrackEstimator, processMCSkimmed, "Count primary charged particles in DQ skimmed MC", false);

  void processMCDirect(McCollision const&, McParticles const& particles)
  {
    mcGlobalMultiplicities(countMCGlobalParticles(particles));
  }

  PROCESS_SWITCH(GlobalTrackEstimator, processMCDirect, "Count primary charged particles in unskimmed MC", false);

  void processDummy(MyEventsBasic const&)
  {
    // do nothing
  }

  PROCESS_SWITCH(GlobalTrackEstimator, processDummy, "Dummy process function", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<GlobalTrackEstimator>(cfgc)};
}
