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
#include "Framework/AnalysisDataModel.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

using namespace o2;
using namespace o2::framework;

// Converts V0 version 001 to 002
struct strarawcentsconverter {
  Produces<aod::StraRawCents_001> straRawCents_001;
  Produces<aod::StraRawCents_003> straRawCents_003;

  void process000to001(aod::StraRawCents_000 const& straRawCents_000)
  {
    for (auto& values : straRawCents_000) {
      straRawCents_001(values.multFT0A(), values.multFT0C(), values.multFV0A(), values.multNTracksPVeta1(), 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
    }
  }
  void process002to003(aod::StraRawCents_002 const& straRawCents_002)
  {
    for (auto& values : straRawCents_002) {
      straRawCents_003(values.multFT0A(),
                       values.multFT0C(),
                       values.multFT0A(),
                       values.multNTracksPVeta1(),
                       0, 0,
                       values.multNTracksITSTPC(),
                       values.multAllTracksTPCOnly(),
                       values.multAllTracksITSTPC(),
                       values.multZNA(),
                       values.multZNC(),
                       values.multZEM1(),
                       values.multZEM2(),
                       values.multZPA(),
                       values.multZPC());
    }
  }

  PROCESS_SWITCH(strarawcentsconverter, process000to001, "from raw 000 to 001", false);
  PROCESS_SWITCH(strarawcentsconverter, process002to003, "from raw 002 to 003", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<strarawcentsconverter>(cfgc)};
}
