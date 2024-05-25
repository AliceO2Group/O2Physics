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
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "DataFormatsFT0/Digit.h"
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
  static const int nBCsPerOrbit = o2::constants::lhc::LHCMaxBunches;
  std::bitset<nBCsPerOrbit> bcPatternB;

  void init(InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    const AxisSpec axisMultT0M{1000, 0., 270000., "T0M multiplicity"};
    const AxisSpec axisMultT0A{1000, 0., 200000., "T0A multiplicity"};
    const AxisSpec axisMultT0C{1000, 0., 70000., "T0C multiplicity"};
    const AxisSpec axisCentT0C{100, 0., 100., "T0C centrality"};
    histos.add("hMultT0M", "", kTH1F, {axisMultT0M});
    histos.add("hMultT0A", "", kTH1F, {axisMultT0A});
    histos.add("hMultT0C", "", kTH1F, {axisMultT0C});
    histos.add("hCentT0C", "", kTH1F, {axisCentT0C});
    histos.add("hMultT0MselTSC", "", kTH1F, {axisMultT0M});
    histos.add("hMultT0MselTCE", "", kTH1F, {axisMultT0M});
    histos.add("hMultT0CselTCE", "", kTH1F, {axisMultT0C});
    histos.add("hCentT0CselTCE", "", kTH1F, {axisCentT0C});
    histos.add("hMultT0CselTVXTCE", "", kTH1F, {axisMultT0C});
    histos.add("hCentT0CselTVXTCE", "", kTH1F, {axisCentT0C});
    histos.add("hMultT0CselTVXTCEB", "", kTH1F, {axisMultT0C});
    histos.add("hCentT0CselTVXTCEB", "", kTH1F, {axisCentT0C});

    histos.add("hCounterTCE", "", kTH1D, {{1, 0., 1.}});
    histos.add("hCounterZNA", "", kTH1D, {{1, 0., 1.}});
    histos.add("hCounterZNC", "", kTH1D, {{1, 0., 1.}});
    histos.add("hCounterZEM", "", kTH1D, {{1, 0., 1.}});
  }

  void process(BCsRun3 const& bcs, aod::Zdcs const&, aod::FT0s const&)
  {
    int runNumber = bcs.iteratorAt(0).runNumber();
    LOGP(info, "runNumber={}", runNumber);
    const char* srun = Form("%d", runNumber);

    if (runNumber != lastRunNumber) {
      lastRunNumber = runNumber;
      int64_t ts = bcs.iteratorAt(0).timestamp();

      auto grplhcif = ccdb->getForTimeStamp<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", ts);
      bcPatternB = grplhcif->getBunchFilling().getBCPattern();

      TList* callst = ccdb->getForTimeStamp<TList>("Centrality/Estimators", ts);
      if (callst == nullptr) {
        LOGF(info, "Centrality calibration is not available in CCDB for run=%d at timestamp=%llu", runNumber, ts);
        return;
      }

      hCalibT0C = reinterpret_cast<TH1*>(callst->FindObject("hCalibZeqFT0C"));
      if (hCalibT0C == nullptr) {
        LOGF(info, "hCalibZeqFT0C histogram is not available for run=%d at timestamp=%llu", runNumber, ts);
        return;
      }
    }

    for (const auto& bc : bcs) {
      if (bc.has_zdc()) {
        float timeZNA = bc.zdc().timeZNA();
        float timeZNC = bc.zdc().timeZNC();
        if (fabs(timeZNA) < 2) {
          histos.get<TH1>(HIST("hCounterZNA"))->Fill(srun, 1);
        }
        if (fabs(timeZNC) < 2) {
          histos.get<TH1>(HIST("hCounterZNC"))->Fill(srun, 1);
        }
        if (fabs(timeZNA) < 2 || fabs(timeZNC) < 2) {
          histos.get<TH1>(HIST("hCounterZEM"))->Fill(srun, 1);
        }
      }

      if (!bc.has_ft0()) {
        continue;
      }
      float multT0A = bc.ft0().sumAmpA();
      float multT0C = bc.ft0().sumAmpC();
      float multT0M = multT0A + multT0C;
      float centT0C = hCalibT0C->GetBinContent(hCalibT0C->FindFixBin(multT0C));
      histos.fill(HIST("hMultT0M"), multT0M);
      histos.fill(HIST("hMultT0A"), multT0A);
      histos.fill(HIST("hMultT0C"), multT0C);
      histos.fill(HIST("hCentT0C"), centT0C);

      if (TESTBIT(bc.ft0().triggerMask(), o2::ft0::Triggers::bitSCen)) { // TSC
        histos.fill(HIST("hMultT0MselTSC"), multT0M);
      }

      if (!TESTBIT(bc.ft0().triggerMask(), o2::ft0::Triggers::bitCen)) { // TCE
        continue;
      }
      histos.fill(HIST("hMultT0MselTCE"), multT0M);
      histos.fill(HIST("hMultT0CselTCE"), multT0C);
      histos.fill(HIST("hCentT0CselTCE"), centT0C);

      if (!TESTBIT(bc.ft0().triggerMask(), o2::ft0::Triggers::bitVertex)) { // TVX
        continue;
      }
      histos.fill(HIST("hMultT0CselTVXTCE"), multT0C);
      histos.fill(HIST("hCentT0CselTVXTCE"), centT0C);

      if (!bcPatternB[bc.globalBC() % nBCsPerOrbit]) { // B-mask
        continue;
      }
      histos.fill(HIST("hMultT0CselTVXTCEB"), multT0C);
      histos.fill(HIST("hCentT0CselTVXTCEB"), centT0C);

      histos.get<TH1>(HIST("hCounterTCE"))->Fill(srun, 1);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<LumiQaTask>(cfgc)};
}
