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
/// \file   udMcParticles2udTracks.cxx
/// \author Roman Lavička
/// \since  2025-04-15
/// \brief  A task to create a reverse index from UDMcParticles to UDTracks
///

#include <vector>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "PWGUD/DataModel/UDTables.h"
#include "PWGUD/DataModel/UDIndex.h"

using namespace o2;
using namespace o2::framework;

struct UDMcParticlesToUDTracks {
  using LabeledTracks = soa::Join<aod::UDTracks, aod::UDMcTrackLabels>;
  Produces<aod::UDMcParticlesToUDTracks> udp2udt;

  std::vector<int> trackIds;

  void init(InitContext&)
  {
  }

  void process(aod::UDMcParticles const& particles)
  {
    if (doprocessIndexingCentral || doprocessIndexingCentralFast) {
      udp2udt.reserve(particles.size());
    }
  }

  void processIndexingCentralFast(aod::UDMcParticles const& mcParticles, LabeledTracks const& tracks)
  {
    // faster version, but will use more memory due to pre-allocation
    std::vector<std::vector<int>> part2track(mcParticles.size());
    for (const auto& track : tracks) {
      if (track.has_udMcParticle())
        part2track[track.udMcParticleId()].push_back(track.globalIndex());
    }
    for (const auto& mcParticle : mcParticles) {
      udp2udt(part2track[mcParticle.globalIndex()]);
    }
  }
  PROCESS_SWITCH(UDMcParticlesToUDTracks, processIndexingCentralFast, "Create reverse index from particles to tracks: more memory use but potentially faster", true);

  void processIndexingCentral(aod::UDMcParticles const&, soa::SmallGroups<LabeledTracks> const& tracks)
  {
    trackIds.clear();
    for (const auto& track : tracks) {
      trackIds.push_back(track.globalIndex());
    }
    udp2udt(trackIds);
  }
  PROCESS_SWITCH(UDMcParticlesToUDTracks, processIndexingCentral, "Create reverse index from particles to tracks", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<UDMcParticlesToUDTracks>(cfgc)};
}
