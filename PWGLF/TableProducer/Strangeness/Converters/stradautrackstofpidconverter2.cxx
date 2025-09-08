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
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;

// converts DauTrackTOFPIDs_000 to _001
struct stradautrackstofpidconverter2 {
  Produces<aod::DauTrackTOFPIDs_001> dautracktofpids;

  void process(aod::DauTrackTOFPIDs_000 const& dauTrackTOFPIDs)
  {
    for (int ii = 0; ii < dauTrackTOFPIDs.size(); ii++) {
      auto dauTrackTOFPID = dauTrackTOFPIDs.rawIteratorAt(ii);
      dautracktofpids(-1, -1, dauTrackTOFPID.tofSignal(), dauTrackTOFPID.tofEvTime(), dauTrackTOFPID.length(), 0.0f);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<stradautrackstofpidconverter2>(cfgc)};
}
 