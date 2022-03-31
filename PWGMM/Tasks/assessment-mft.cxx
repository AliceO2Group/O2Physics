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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

#include "MathUtils/Utils.h"
#include "CommonConstants/LHCConstants.h"

using namespace o2;
using namespace o2::framework;

struct AssessmentMFT {

  HistogramRegistry registry{
    "registry",
    {

      {"TracksPhiEta", "; #varphi; #eta; tracks", {HistType::kTH2F, {{600, 0, 2 * M_PI}, {35, -4.5, -1.}}}}, //
      {"TracksTime", "; time; #count", {HistType::kTH1D, {{6000000, 0, 60000}}}},                            //
    }                                                                                                        //
  };

  void process(aod::Collisions::iterator const& collision, aod::MFTTracks const& tracks, aod::BCs const& bcs)
  {
    for (auto& track : tracks) {
      float phi = track.phi();
      o2::math_utils::bringTo02Pi(phi);
      registry.fill(HIST("TracksPhiEta"), phi, track.eta());

      auto collisionBCID = collision.bcId();
      auto collisionBC = (bcs.iteratorAt(collisionBCID)).globalBC();
      //printf("collisionBC %d\n", collisionBC);
      double seconds = (collisionBC * o2::constants::lhc::LHCBunchSpacingNS + track.trackTime()) / 1e9;
      //this is in s
      //printf("seconds  %f\n", seconds);
      registry.fill(HIST("TracksTime"), seconds);
    }
    //printf("minTrackTime = %f, maxTrackTime = %f\n", minTrackTime, maxTrackTime);
  }
};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<AssessmentMFT>(cfgc)};
}
