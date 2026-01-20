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
// ========================
//
// This code runs loop over ULS ee pars for virtual photon QC.
//    Please write to: daiki.sekihata@cern.ch

#include "PWGEM/Dilepton/DataModel/dileptonTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct trackConverter1 {
  Produces<aod::EMPrimaryTracks_001> track_001;

  void process(aod::EMPrimaryTracks_000 const& tracks)
  {
    for (const auto& track : tracks) {
      track_001(track.collisionId(),
                track.trackId(),
                track.sign() / track.pt(),
                track.eta(),
                track.phi(),
                track.trackBit());
    } // end of track loop
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<trackConverter1>(cfgc, TaskName{"track-converter1"})};
}
