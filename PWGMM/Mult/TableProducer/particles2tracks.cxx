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

///
/// \file   particles2tracks.cxx
/// \author Anton Alkin
/// \since  2022-04-27
/// \brief  A task to create a reverse index from McParticles to Tracks
///

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"

#include "Index.h"

using namespace o2;
using namespace o2::framework;

struct ParticlesToTracks {
  using LabeledTracks = soa::Join<aod::Tracks, aod::McTrackLabels>;
  Produces<aod::ParticlesToTracks> p2t;

  using LabeledMFTTracks = soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels>;
  Produces<aod::ParticlesToMftTracks> p2tmft;

  std::vector<int> trackIds;

  void init(InitContext&)
  {
  }

  void process(aod::McParticles const& particles)
  {
    if (doprocessIndexingCentral) {
      p2t.reserve(particles.size());
    }
    if (doprocessIndexingFwd) {
      p2tmft.reserve(particles.size());
    }
  }

  void processIndexingCentral(aod::McParticle const&, soa::SmallGroups<LabeledTracks> const& tracks)
  {
    trackIds.clear();
    for (auto& track : tracks) {
      trackIds.push_back(track.globalIndex());
    }
    p2t(trackIds);
  }

  PROCESS_SWITCH(ParticlesToTracks, processIndexingCentral, "Create reverse index from particles to tracks", false);

  void processIndexingFwd(aod::McParticle const&, soa::SmallGroups<LabeledMFTTracks> const& tracks)
  {
    trackIds.clear();
    for (auto& track : tracks) {
      trackIds.push_back(track.globalIndex());
    }
    p2tmft(trackIds);
  }

  PROCESS_SWITCH(ParticlesToTracks, processIndexingFwd, "Create reverse index from particles to tracks", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<ParticlesToTracks>(cfgc)};
}
