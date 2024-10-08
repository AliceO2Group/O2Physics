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

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "PWGLF/Utils/inelGt.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct mcInelGt0 {
  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    registry.add("selected", "selected", HistType::kTH1D, {{10, 0, 10}});
    registry.add("ITScontributors", "ITScontributors", HistType::kTH2D, {{10, 0, 10}, {10, 0, 10}});
  }

  Service<o2::framework::O2DatabasePDG> pdgDB;
  SliceCache cache;
  Preslice<o2::aod::Tracks> perCollision = o2::aod::track::collisionId;
  void process(o2::aod::McCollision const& /*collisionMC*/,
               o2::soa::SmallGroups<o2::soa::Join<o2::aod::Collisions, o2::aod::McCollisionLabels>> const& collisions,
               const soa::Join<o2::aod::TracksIU, o2::aod::TracksExtra, o2::aod::McTrackLabels>& tracks,
               o2::aod::McParticles const& mcParticles)
  {
    if (pwglf::isINELgtNmc(mcParticles, 0, pdgDB)) {
      registry.fill(HIST("selected"), 0.5);
    }
    if (pwglf::isINELgtNmc(mcParticles, 1, pdgDB)) {
      registry.fill(HIST("selected"), 1.5);
    }
    for (auto const& collision : collisions) {
      const auto groupedTracks = tracks.sliceBy(perCollision, collision.globalIndex());
      for (auto const& track : groupedTracks) {
        if (!track.isPVContributor()) {
          continue;
        }
        if (std::abs(track.eta()) > 1) {
          LOG(info) << "Track with eta > 1: " << track.eta()
                    << (track.hasTPC()
                          ? "hasTPC"
                          : "no TPC ")
                    << (track.hasITS()
                          ? "hasITS"
                          : "no ITS");
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<mcInelGt0>(cfgc)}; }
