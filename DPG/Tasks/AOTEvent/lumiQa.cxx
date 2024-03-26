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
#include "CCDB/BasicCCDBManager.h"
#include "Framework/HistogramRegistry.h"
#include "TList.h"
#include "TH1.h"

using namespace o2;
using namespace o2::framework;
using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::Run3MatchedToBCSparse>;

struct LumiQaTask {
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  int lastRunNumber = -1;
  TH1* hCalibT0C = nullptr;

  void init(InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    const AxisSpec axisMultT0C{1000, 0., 70000., "T0C multiplicity"};
    const AxisSpec axisCentT0C{1000, 0., 100., "T0C centrality"};
    histos.add("hMultT0C", "", kTH1F, {axisMultT0C});
    histos.add("hCentT0C", "", kTH1F, {axisCentT0C});
  }

  void process(BCsRun3 const& bcs, aod::Zdcs const& zdcs, aod::FT0s const& ft0s)
  {
    int runNumber = bcs.iteratorAt(0).runNumber();
    LOGP(info, "runNumber={}", runNumber);
    if (runNumber != lastRunNumber) {
      TList* callst = ccdb->getForTimeStamp<TList>("Centrality/Estimators", bcs.iteratorAt(0).timestamp());
      lastRunNumber = runNumber;
      hCalibT0C = reinterpret_cast<TH1*>(callst->FindObject("hCalibZeqFT0C"));
    }
    if (!hCalibT0C) {
      return;
    }

    for (const auto& bc : bcs) {
      if (!bc.has_ft0()) {
        continue;
      }
      float multT0C = bc.ft0().sumAmpC();
      float centT0C = hCalibT0C->GetBinContent(hCalibT0C->FindFixBin(multT0C));
      histos.fill(HIST("hMultT0C"), multT0C);
      histos.fill(HIST("hCentT0C"), centT0C);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<LumiQaTask>(cfgc)};
}
