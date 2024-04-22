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

struct strarawcentsconverter {
  Produces<aod::StraRawCents_001> straRawCents_001;
  Produces<aod::StraRawCents_003> straRawCents_003;

  void processStraRawCents000to001(aod::StraRawCents_000 const& straRawCents_000)
  {
    for (auto& values : straRawCents_000) {
      straRawCents_001(values.multFT0A(), values.multFT0C(), values.multFV0A(), values.multNTracksPVeta1(), 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
    }
  }
  void processStraRawCents002to003(aod::StraRawCents_002 const& straRawCents_002)
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

  PROCESS_SWITCH(strarawcentsconverter, processStraRawCents000to001, "from StraRawCents 000 to 001", false);
  PROCESS_SWITCH(strarawcentsconverter, processStraRawCents002to003, "from StraRawCents 002 to 003", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<strarawcentsconverter>(cfgc)};
}
