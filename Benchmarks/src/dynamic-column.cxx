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

#include "tables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct ProduceDynamicExtension {
  Produces<aod::ExtTracksDM> etdm;
  Configurable<float> factor{"factor", 1.0f, "Configurable factor"};

  void process(aod::TracksIU const& tracks)
  {
    auto tracksD = soa::Attach<aod::TracksIU, aod::extensions::Direct<aod::track::X, aod::track::Y, aod::track::Z>>(tracks);
    for (auto& track : tracksD) {
      etdm(track.direct((float)factor));
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return {adaptAnalysisTask<ProduceDynamicExtension>(cfgc)};
}
