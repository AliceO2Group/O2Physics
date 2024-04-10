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
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"

using namespace o2;
using namespace o2::framework;

// Converts V0 version 001 to 002
struct strangeDerivedConverter {
  Produces<aod::StraRawCents_001> straRawCents_001;
  Produces<aod::StraRawCents_003> straRawCents_003;
  Produces<aod::DauTrackTOFPIDs> dautracktofpids;

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

  void processGenerateDauTracksTOFPIDs(soa::Join<aod::V0Cores, aod::V0Extras, aod::V0TOFs> const& v0s, soa::Join<aod::CascCores, aod::CascExtras, aod::CascTOFs> const& cascs, aod::DauTrackExtras const& dauTracks)
  {
    // prepare arrays with the relevant information
    std::vector<float> lLengths, lTOFSignals, lTOFEvTimes;
    lLengths.reserve(dauTracks.size());
    lTOFSignals.reserve(dauTracks.size());
    lTOFEvTimes.reserve(dauTracks.size());
    for (int ii = 0; ii < dauTracks.size(); ii++) {
      lLengths[ii] = 1e+6;
      lTOFSignals[ii] = -1e+3f;
      lTOFEvTimes[ii] = -1e+3f;
    }
    for (auto& v0 : v0s) {
      lLengths[v0.posTrackExtraId()] = v0.posTOFLengthToPV();
      lTOFSignals[v0.posTrackExtraId()] = v0.posTOFSignal();
      lTOFEvTimes[v0.posTrackExtraId()] = v0.posTOFEventTime();
      lLengths[v0.negTrackExtraId()] = v0.negTOFLengthToPV();
      lTOFSignals[v0.negTrackExtraId()] = v0.negTOFSignal();
      lTOFEvTimes[v0.negTrackExtraId()] = v0.negTOFEventTime();
    }
    for (auto& casc : cascs) {
      lLengths[casc.posTrackExtraId()] = casc.posTOFLengthToPV();
      lTOFSignals[casc.posTrackExtraId()] = casc.posTOFSignal();
      lTOFEvTimes[casc.posTrackExtraId()] = casc.posTOFEventTime();
      lLengths[casc.negTrackExtraId()] = casc.negTOFLengthToPV();
      lTOFSignals[casc.negTrackExtraId()] = casc.negTOFSignal();
      lTOFEvTimes[casc.negTrackExtraId()] = casc.negTOFEventTime();
      lLengths[casc.bachTrackExtraId()] = casc.bachTOFLengthToPV();
      lTOFSignals[casc.bachTrackExtraId()] = casc.bachTOFSignal();
      lTOFEvTimes[casc.bachTrackExtraId()] = casc.bachTOFEventTime();
    }
    for (int ii = 0; ii < dauTracks.size(); ii++) {
      dautracktofpids(lTOFSignals[ii], lTOFEvTimes[ii], lLengths[ii]);
    }
  }

  PROCESS_SWITCH(strangeDerivedConverter, processStraRawCents000to001, "from StraRawCents 000 to 001", false);
  PROCESS_SWITCH(strangeDerivedConverter, processStraRawCents002to003, "from StraRawCents 002 to 003", false);
  PROCESS_SWITCH(strangeDerivedConverter, processGenerateDauTracksTOFPIDs, "from DauTrackTOFPIDs to V0TOFs/CascTOFs", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<strangeDerivedConverter>(cfgc)};
}
