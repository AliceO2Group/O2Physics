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
/// \file   ECALqa.cxx
/// \author Nicolo' Jacazio
/// \since  14/09/2021
/// \brief  Task to use the ALICE3 ECAL table
///

// O2 includes
#include "Framework/AnalysisTask.h"
#include "ALICE3/DataModel/ECAL.h"
#include "Common/Core/MC.h"
#include "Common/Core/PID/PIDResponse.h"
#include "ReconstructionDataFormats/PID.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::track;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{

namespace indices
{
DECLARE_SOA_INDEX_COLUMN(Track, track);
DECLARE_SOA_INDEX_COLUMN(ECAL, ecal);
} // namespace indices

DECLARE_SOA_INDEX_TABLE_USER(ECALTracksIndex, Tracks, "ECALTRK", indices::TrackId, indices::ECALId);
} // namespace o2::aod

struct ecalIndexBuilder { // Builder of the ECAL-track index linkage
  Builds<o2::aod::ECALTracksIndex> ind;
  void init(o2::framework::InitContext&)
  {
  }
};

struct ecalQaMc {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::QAObject};

  void init(o2::framework::InitContext&)
  {
    histos.add("energy", ";Energy;Entries", HistType::kTH1F, {{1000, 0, 10000}});
  }

  using Trks = soa::Join<aod::Tracks, aod::ECALTracksIndex, aod::TracksExtra>;
  void process(const aod::McParticles& mcParticles,
               const Trks& tracks,
               const aod::McTrackLabels& labels,
               const aod::ECALs&,
               const aod::Collisions& colls)
  {
    for (auto& track : tracks) {
      if (!track.has_ecal())
        continue;
      histos.fill(HIST("energy"), track.ecal().energy());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  auto workflow = WorkflowSpec{adaptAnalysisTask<ecalIndexBuilder>(cfg)};
  workflow.push_back(adaptAnalysisTask<ecalQaMc>(cfg));
  return workflow;
}
