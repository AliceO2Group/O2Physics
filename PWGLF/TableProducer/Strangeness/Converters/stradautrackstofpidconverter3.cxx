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
/// \file stradautrackstofpidconverter3.cxx
/// \brief Converts DauTrackTOFPID_001 and StraEvTimes_000 into DauTrackTOFPID_002 and StraEvTimes_001
///
/// \author David Dobrigkeit Chinellato <david.dobrigkeit.chinellato@cern.ch>, Austrian Academy of Sciences & MBI
/// \author Romain Schotter <romain.schotter@cern.ch>, Austrian Academy of Sciences & MBI

#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;

// converts DauTrackTOFPIDs_001 to _002
struct stradautrackstofpidconverter3 {
  Produces<aod::DauTrackTOFPIDs_002> dautracktofpids;
  Produces<aod::StraEvTimes_001> straEvTimes;

  void process(aod::DauTrackTOFPIDs_001 const& dauTracks, aod::StraEvTimes_000 const& straEvTimes_000)
  {
    // create new TOFPIDs
    dautracktofpids.reserve(dauTracks.size());
    for (const auto& dauTrack : dauTracks) {
      dautracktofpids(
        -1,
        -1,
        dauTrack.tofSignal(),
        dauTrack.tofEvTime(),
        999.0f, /*dummy event time error for TOF*/
        dauTrack.length(),
        0.0f);
    }
    straEvTimes.reserve(straEvTimes_000.size());
    for (const auto& value : straEvTimes_000) {
      straEvTimes(value.eventTime(), 999.0f /*dummy event time error for TOF*/);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<stradautrackstofpidconverter3>(cfgc)};
}
