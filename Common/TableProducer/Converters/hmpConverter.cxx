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

// HMPID converter to new format
// to be used with Run 2 converted data and older AO2Ds

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/runDataProcessing.h>

using namespace o2;
using namespace o2::framework;

struct hmpConverter {
  Produces<aod::HMPID_001> HMPID_001;

  void process(aod::HMPID_000 const& hmpLegacy, aod::Tracks const&)
  {
    for (auto& hmpData : hmpLegacy) {

      float phots[] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
      auto trackid = hmpData.trackId();
      auto hmpidSignal = hmpData.hmpidSignal();
      auto hmpidXTrack = -999.; // dummy
      auto hmpidYTrack = -999.; // dummy
      auto hmpidXMip = -999.;   // dummy
      auto hmpidYMip = -999.;   // dummy
      auto hmpidNPhotons = hmpData.hmpidNPhotons();
      auto hmpidQMip = hmpData.hmpidQMip();
      auto hmpidClusSize = -999; // dummy
      auto hmpidMom = -999;      // dummy
      auto hmpidPhotsCharge = phots;

      HMPID_001(trackid,
                hmpidSignal,
                hmpidXTrack,
                hmpidYTrack,
                hmpidXMip,
                hmpidYMip,
                hmpidNPhotons,
                hmpidQMip,
                hmpidClusSize,
                hmpidMom,
                hmpidPhotsCharge);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<hmpConverter>(cfgc),
  };
}
