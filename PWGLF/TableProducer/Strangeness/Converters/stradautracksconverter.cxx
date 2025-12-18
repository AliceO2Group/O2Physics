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
struct stradautracksconverter {
  Produces<aod::DauTrackTOFPIDs_000> dautracktofpids;

  void process(soa::Join<aod::V0Cores, aod::V0Extras, aod::V0TOFs> const& v0s, soa::Join<aod::CascCores, aod::CascExtras, aod::CascTOFs> const& cascs, aod::DauTrackExtras const& dauTracks)
  {
    // prepare arrays with the relevant information
    std::vector<float> lLengths(dauTracks.size(), 1.e+6), lTOFSignals(dauTracks.size(), -1e+3f), lTOFEvTimes(dauTracks.size(), -1e+3f);
    for (const auto& v0 : v0s) {
      lLengths[v0.posTrackExtraId()] = v0.posTOFLengthToPV();
      lTOFSignals[v0.posTrackExtraId()] = v0.posTOFSignal();
      lTOFEvTimes[v0.posTrackExtraId()] = v0.posTOFEventTime();
      lLengths[v0.negTrackExtraId()] = v0.negTOFLengthToPV();
      lTOFSignals[v0.negTrackExtraId()] = v0.negTOFSignal();
      lTOFEvTimes[v0.negTrackExtraId()] = v0.negTOFEventTime();
    }
    for (const auto& casc : cascs) {
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
    for (unsigned int ii = 0; ii < dauTracks.size(); ii++) {
      dautracktofpids(lTOFSignals[ii], lTOFEvTimes[ii], lLengths[ii]);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<stradautracksconverter>(cfgc)};
}
