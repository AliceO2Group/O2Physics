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
struct stradautrackstpcpidconverter {
  Produces<aod::DauTrackTPCPIDs_001> dautrackpcpids;

  void process(aod::DauTrackTPCPIDs_000 const& v000s)
  {
    for (int ii = 0; ii < v000s.size(); ii++) {
      auto dauTrackTPCPID = v000s.rawIteratorAt(ii);
      dautrackpcpids(dauTrackTPCPID.tpcSignal(),
                     aod::dautrack::packing::packInInt8(dauTrackTPCPID.tpcNSigmaEl()),
                     aod::dautrack::packing::packInInt8(dauTrackTPCPID.tpcNSigmaPi()),
                     aod::dautrack::packing::packInInt8(dauTrackTPCPID.tpcNSigmaKa()),
                     aod::dautrack::packing::packInInt8(dauTrackTPCPID.tpcNSigmaPr()));
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<stradautrackstpcpidconverter>(cfgc)};
}
