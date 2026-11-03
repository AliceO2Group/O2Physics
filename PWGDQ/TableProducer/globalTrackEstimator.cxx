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

#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/MCSignal.h"
#include "PWGDQ/Core/MCSignalLibrary.h"
#include "PWGDQ/Core/VarManager.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <string>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

using MyEventsBasic = o2::soa::Join<o2::aod::ReducedEvents, o2::aod::ReducedEventsExtended>;
using MyBarrelTracks = o2::soa::Join<o2::aod::ReducedTracks, o2::aod::ReducedTracksBarrel, o2::aod::ReducedTracksBarrelCov, o2::aod::ReducedTracksBarrelPID>;

constexpr static uint32_t GkEventFillMap = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended;
constexpr static uint32_t GkTrackFillMap = VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelCov;

struct GlobalTrackEstimator {
  Produces<o2::aod::ReducedEventsGlobalTrk> reducedEventsGlobalTrk;
  Produces<o2::aod::ReducedEventsMultAzimuth> reducedEventsMultAzimuth;

  Configurable<float> cfgMaxDCAxy{"cfgMaxDCAxy", 2.0, "Maximum DCAxy"};
  Configurable<float> cfgMaxDCAz{"cfgMaxDCAz", 2.0, "Maximum DCAz"};
  Configurable<float> cfgMinPt{"cfgMinPt", 0.1, "Minimum required pT"};
  Configurable<float> cfgEtaRange{"cfgEtaRange", 0.8, "Pseudorapidity range"};
  Configurable<float> cfgMinTPCCluster{"cfgMinTPCCluster", 70, "Minimum number of crossed TPC rows"};
  Configurable<float> cfgRatioTPCClusterOverFindableCluster{"cfgRatioTPCClusterOverFindableCluster", 0.8, "Minimum ratio of crossed TPC rows over findable clusters"};
  Configurable<float> cfgMaxChi2TPC{"cfgMaxChi2TPC", 4.0, "Maximum chi2 per cluster in TPC"};
  Configurable<float> cfgMaxChi2ITS{"cfgMaxChi2ITS", 36.0, "Maximum chi2 per cluster in ITS"};
  Configurable<bool> cfgRequireITSibAny{"cfgRequireITSibAny", true, "Pseudorapidity range"};

  o2::framework::Preslice<o2::aod::ReducedTracksAssoc> preslice = o2::aod::reducedtrack_association::reducedeventId;

  o2::framework::expressions::Filter trackFilter = (o2::aod::reducedtrack::pt >= cfgMinPt) && (nabs(o2::aod::reducedtrack::eta) <= cfgEtaRange) && (o2::aod::track::tpcChi2NCl <= cfgMaxChi2TPC) && (o2::aod::track::itsChi2NCl <= cfgMaxChi2ITS);

  void init(InitContext& context)
  {
    LOGF(info, "Started global track counting");

    bool isAnyProcessEnabled = context.mOptions.get<bool>("processGlobalTracksSkimmed") || context.mOptions.get<bool>("processGlobalTracksSkimmed");
    bool isDummyEnabled = context.mOptions.get<bool>("processDummy");

    if (isDummyEnabled) {
      if (isAnyProcessEnabled) {
        LOG(fatal) << "You have enabled both the dummy process and at least one normal process function! Check your config!";
      }
      return;
    }

    if (!isAnyProcessEnabled) {
      return;
    }
  }

  void processGlobalTracksSkimmed(MyEventsBasic const& events, o2::aod::ReducedTracksAssoc const& assocs, o2::soa::Filtered<MyBarrelTracks> const& tracks)
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

        auto track = assoc.template reducedtrack_as<o2::soa::Filtered<MyBarrelTracks>>();

        VarManager::FillTrack<GkTrackFillMap>(track);
        VarManager::FillTrackCollision<GkTrackFillMap>(track, event);

        const int nFindable = track.tpcNClsFindable();
        const int nCrossedRows = track.tpcNClsCrossedRows();

        if (nCrossedRows < static_cast<int>(cfgMinTPCCluster.value)) {
          continue;
        }

        if (nFindable <= 0) {
          continue;
        }

        const float crossedOverFindable = static_cast<float>(nCrossedRows) / static_cast<float>(nFindable);

        if (crossedOverFindable < cfgRatioTPCClusterOverFindableCluster.value) {
          continue;
        }

        if (cfgRequireITSibAny.value && ((track.itsClusterMap() & 0x7) == 0)) {
          continue;
        }

        if (std::abs(VarManager::fgValues[VarManager::kTrackDCAxy]) > cfgMaxDCAxy.value || std::abs(VarManager::fgValues[VarManager::kTrackDCAz]) > cfgMaxDCAz.value) {
          continue;
        }

        ++nTrkGlobal;
      }

      reducedEventsGlobalTrk(nTrkGlobal);
    }
  }

  PROCESS_SWITCH(GlobalTrackEstimator, processGlobalTracksSkimmed, "Count global tracks from skimmed data model", false);

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
