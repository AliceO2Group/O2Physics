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
#include "Framework/ConfigParamSpec.h"

using namespace o2;
using namespace o2::framework;

// custom configurable for switching between run2 and run3 selection types
void customize(std::vector<ConfigParamSpec>& workflowOptions)
{
  workflowOptions.push_back(ConfigParamSpec{"selection-run", VariantType::Int, 2, {"selection type: 2 - run 2, 3 - run 3"}});
}

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "iostream"

struct MultiplicityTableTaskIndexed {
  Produces<aod::Mults> mult;
  Partition<soa::Join<aod::Tracks, aod::TracksExtra>> run2tracklets = (aod::track::trackType == static_cast<uint8_t>(o2::aod::track::TrackTypeEnum::Run2Tracklet));
  Partition<soa::Join<aod::Tracks, aod::TracksExtra>> tracksWithTPC = (aod::track::tpcNClsFindable > (uint8_t)0);
  void processRun2(aod::Run2MatchedSparse::iterator const& collision, soa::Join<aod::Tracks, aod::TracksExtra> const& tracksExtra, aod::BCs const&, aod::Zdcs const&, aod::FV0As const& fv0as, aod::FV0Cs const& fv0cs, aod::FT0s const& ft0s)
  {
    float multV0A = -1.f;
    float multV0C = -1.f;
    float multT0A = -1.f;
    float multT0C = -1.f;
    float multZNA = -1.f;
    float multZNC = -1.f;
    int multTracklets = run2tracklets.size();
    int multTPC = tracksWithTPC.size();

    if (collision.has_fv0a()) {
      for (auto amplitude : collision.fv0a().amplitude()) {
        multV0A += amplitude;
      }
    }
    if (collision.has_fv0c()) {
      for (auto amplitude : collision.fv0c().amplitude()) {
        multV0C += amplitude;
      }
    }
    if (collision.has_ft0()) {
      auto ft0 = collision.ft0();
      for (auto amplitude : ft0.amplitudeA()) {
        multT0A += amplitude;
      }
      for (auto amplitude : ft0.amplitudeC()) {
        multT0C += amplitude;
      }
    }
    if (collision.has_zdc()) {
      auto zdc = collision.zdc();
      multZNA = zdc.energyCommonZNA();
      multZNC = zdc.energyCommonZNC();
    }

    LOGF(debug, "multV0A=%5.0f multV0C=%5.0f multT0A=%5.0f multT0C=%5.0f multZNA=%6.0f multZNC=%6.0f multTracklets=%i multTPC=%i", multV0A, multV0C, multT0A, multT0C, multZNA, multZNC, multTracklets, multTPC);
    mult(multV0A, multV0C, multT0A, multT0C, multZNA, multZNC, multTracklets, multTPC);
  }
  PROCESS_SWITCH(MultiplicityTableTaskIndexed, processRun2, "Produce Run 2 multiplicity tables", true);

  void processRun3(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Join<aod::Tracks, aod::TracksExtra> const& tracksExtra, aod::BCs const& bcs, aod::Zdcs const& zdcs, aod::FV0As const& fv0as, aod::FT0s const& ft0s)
  {
    float multV0A = -1.f;
    float multV0C = -1.f;
    float multT0A = -1.f;
    float multT0C = -1.f;
    float multZNA = -1.f;
    float multZNC = -1.f;
    int multTracklets = -1;
    int multTPC = tracksWithTPC.size();

    // using FT0 row index from event selection task
    if (collision.has_foundFT0()) {
      auto ft0 = collision.foundFT0();
      for (auto amplitude : ft0.amplitudeA()) {
        multT0A += amplitude;
      }
      for (auto amplitude : ft0.amplitudeC()) {
        multT0C += amplitude;
      }
    }

    LOGF(debug, "multV0A=%5.0f multV0C=%5.0f multT0A=%5.0f multT0C=%5.0f multZNA=%6.0f multZNC=%6.0f multTracklets=%i multTPC=%i", multV0A, multV0C, multT0A, multT0C, multZNA, multZNC, multTracklets, multTPC);
    mult(multV0A, multV0C, multT0A, multT0C, multZNA, multZNC, multTracklets, multTPC);
  }
  PROCESS_SWITCH(MultiplicityTableTaskIndexed, processRun3, "Produce Run 3 multiplicity tables", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<MultiplicityTableTaskIndexed>(cfgc, TaskName{"multiplicity-table"})};
}
