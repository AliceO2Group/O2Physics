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
/// \brief Load all tables in the data model to check if they can be read correctly
/// \author
/// \since

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/ConfigParamSpec.h>
#include <Framework/Variant.h>

#include <TH1.h>

#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// custom configurable for switching between run2 and run3 selection types
void customize(std::vector<ConfigParamSpec>& workflowOptions)
{
  workflowOptions.push_back(ConfigParamSpec{"selection-run", VariantType::Int, 2, {"selection type: 2 - run 2, 3 - run 3"}});
  // workflowOptions.push_back(ConfigParamSpec{"isMC", VariantType::Bool, false, {"Check also MC tables if set"}});
}

#include <Framework/runDataProcessing.h>

template <typename Table>
struct LoadTable {
  OutputObj<TH1F> counter{TH1F("counter", "counter", 2, 0., 2)};
  void process(Table const& table)
  {
    LOGF(info, "Table has %d entries", table.size());
    counter->Fill(0.5);
    counter->Fill(1.5, table.size());
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  // TODO MC -->  if (cfgc.options().get<int>("isMC") == true) {

  if (cfgc.options().get<int>("selection-run") == 2) {
    return WorkflowSpec{
      adaptAnalysisTask<LoadTable<aod::BCs>>(cfgc, TaskName("BCs")),
      adaptAnalysisTask<LoadTable<aod::Collisions>>(cfgc, TaskName("Collisions")),

      adaptAnalysisTask<LoadTable<aod::Tracks>>(cfgc, TaskName("Tracks")),
      adaptAnalysisTask<LoadTable<aod::TracksCov>>(cfgc, TaskName("TracksCov")),
      adaptAnalysisTask<LoadTable<aod::TracksExtra>>(cfgc, TaskName("TracksExtra")),

      adaptAnalysisTask<LoadTable<aod::FwdTracks>>(cfgc, TaskName("FwdTracks")),
      adaptAnalysisTask<LoadTable<aod::FwdTracksCov>>(cfgc, TaskName("FwdTracksCov")),

      adaptAnalysisTask<LoadTable<aod::Calos>>(cfgc, TaskName("Calos")),
      adaptAnalysisTask<LoadTable<aod::CaloTriggers>>(cfgc, TaskName("CaloTriggers")),

      adaptAnalysisTask<LoadTable<aod::HMPIDs>>(cfgc, TaskName("HMPIDs")),
      adaptAnalysisTask<LoadTable<aod::Zdcs>>(cfgc, TaskName("Zdcs")),
      adaptAnalysisTask<LoadTable<aod::FV0As>>(cfgc, TaskName("FV0As")),
      adaptAnalysisTask<LoadTable<aod::FT0s>>(cfgc, TaskName("FT0s")),
      adaptAnalysisTask<LoadTable<aod::FDDs>>(cfgc, TaskName("FDDs")),

      adaptAnalysisTask<LoadTable<aod::V0s>>(cfgc, TaskName("V0s")),
      adaptAnalysisTask<LoadTable<aod::Cascades>>(cfgc, TaskName("Cascades")),

      adaptAnalysisTask<LoadTable<aod::Run2BCInfos>>(cfgc, TaskName("Run2BCInfos")),
      adaptAnalysisTask<LoadTable<aod::FV0Cs>>(cfgc, TaskName("FV0Cs")),
    };
  } else if (cfgc.options().get<int>("selection-run") == 3) {
    return WorkflowSpec{
      adaptAnalysisTask<LoadTable<aod::BCs>>(cfgc, TaskName("BCs")),
      adaptAnalysisTask<LoadTable<aod::Collisions>>(cfgc, TaskName("Collisions")),

      adaptAnalysisTask<LoadTable<aod::TracksIU>>(cfgc, TaskName("TracksIU")),
      adaptAnalysisTask<LoadTable<aod::TracksCovIU>>(cfgc, TaskName("TracksCovIU")),
      adaptAnalysisTask<LoadTable<aod::TracksExtra>>(cfgc, TaskName("TracksExtra")),

      adaptAnalysisTask<LoadTable<aod::FwdTracks>>(cfgc, TaskName("FwdTracks")),
      adaptAnalysisTask<LoadTable<aod::FwdTracksCov>>(cfgc, TaskName("FwdTracksCov")),

      adaptAnalysisTask<LoadTable<aod::Calos>>(cfgc, TaskName("Calos")),
      adaptAnalysisTask<LoadTable<aod::CaloTriggers>>(cfgc, TaskName("CaloTriggers")),

      adaptAnalysisTask<LoadTable<aod::HMPIDs>>(cfgc, TaskName("HMPIDs")),
      adaptAnalysisTask<LoadTable<aod::Zdcs>>(cfgc, TaskName("Zdcs")),
      adaptAnalysisTask<LoadTable<aod::FV0As>>(cfgc, TaskName("FV0As")),
      adaptAnalysisTask<LoadTable<aod::FT0s>>(cfgc, TaskName("FT0s")),
      adaptAnalysisTask<LoadTable<aod::FDDs>>(cfgc, TaskName("FDDs")),

      adaptAnalysisTask<LoadTable<aod::V0s>>(cfgc, TaskName("V0s")),
      adaptAnalysisTask<LoadTable<aod::Cascades>>(cfgc, TaskName("Cascades")),

      adaptAnalysisTask<LoadTable<aod::AmbiguousTracks>>(cfgc, TaskName("AmbiguousTracks")),
      adaptAnalysisTask<LoadTable<aod::MFTTracks>>(cfgc, TaskName("MFTTracks")),
      adaptAnalysisTask<LoadTable<aod::AmbiguousMFTTracks>>(cfgc, TaskName("AmbiguousMFTTracks")),
      adaptAnalysisTask<LoadTable<aod::AmbiguousFwdTracks>>(cfgc, TaskName("AmbiguousFwdTracks")),
    };
  } else {
    LOGF(fatal, "Invalid setting for run: %d", cfgc.options().get<int>("selection-run"));
    return WorkflowSpec{};
  }
}
