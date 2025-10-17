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
struct stradautrackstofpidconverter {
  Produces<aod::DauTrackTOFPIDs> dautracktofpids;

  void process(soa::Join<aod::V0Cores, aod::V0Extras, aod::V0TOFs> const& v0s, soa::Join<aod::CascCores, aod::CascExtras, aod::CascTOFs> const& cascs, aod::DauTrackExtras const& dauTracks)
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
      dautracktofpids(-1, -1, lTOFSignals[ii], lTOFEvTimes[ii], lLengths[ii], 0.0f);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<stradautrackstofpidconverter>(cfgc)};
}