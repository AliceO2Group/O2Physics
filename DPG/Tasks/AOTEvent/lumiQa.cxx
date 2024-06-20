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
  std::bitset<nBCsPerOrbit> bcPatternA;
  std::bitset<nBCsPerOrbit> bcPatternC;
  std::bitset<nBCsPerOrbit> beamPatternA;
  std::bitset<nBCsPerOrbit> beamPatternC;

  void init(InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    const AxisSpec axisMultZNA{2000, 0., 400., "ZNA multiplicity"};
    const AxisSpec axisMultZNC{2000, 0., 400., "ZNC multiplicity"};
    const AxisSpec axisMultT0M{1000, 0., 270000., "T0M multiplicity"};
    const AxisSpec axisMultT0A{1000, 0., 200000., "T0A multiplicity"};
    const AxisSpec axisMultT0C{1000, 0., 70000., "T0C multiplicity"};
    const AxisSpec axisCentT0C{100, 0., 100., "T0C centrality"};
    const AxisSpec axisTime{700, -35., 35., "time (ns)"};
    const AxisSpec axisMultChannelT0A{5000, 0., 5000., "T0A channel multiplicity"};
    const AxisSpec axisMultChannelT0C{5000, 0., 5000., "T0C channel multiplicity"};
    int nChannelsT0A = 96;
    int nChannelsT0C = 112;
    const AxisSpec axisChannelsT0A{nChannelsT0A, 0, static_cast<double>(nChannelsT0A)};
    const AxisSpec axisChannelsT0C{nChannelsT0C, 0, static_cast<double>(nChannelsT0C)};

    histos.add("hMultZNA", "", kTH1F, {axisMultZNA});
    histos.add("hMultZNC", "", kTH1F, {axisMultZNC});
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

    histos.add("hTimeZNA", "", kTH1F, {axisTime});
    histos.add("hTimeZNC", "", kTH1F, {axisTime});
    histos.add("hTimeZNAselTVX", "", kTH1F, {axisTime});
    histos.add("hTimeZNCselTVX", "", kTH1F, {axisTime});
    histos.add("hTimeZNAselB", "", kTH1F, {axisTime});
    histos.add("hTimeZNCselB", "", kTH1F, {axisTime});
    histos.add("hTimeZNAselA", "", kTH1F, {axisTime});
    histos.add("hTimeZNCselA", "", kTH1F, {axisTime});
    histos.add("hTimeZNAselC", "", kTH1F, {axisTime});
    histos.add("hTimeZNCselC", "", kTH1F, {axisTime});

    histos.add("hCounterTCE", "", kTH1D, {{1, 0., 1.}});
    histos.add("hCounterZNA", "", kTH1D, {{1, 0., 1.}});
    histos.add("hCounterZNC", "", kTH1D, {{1, 0., 1.}});
    histos.add("hCounterZEM", "", kTH1D, {{1, 0., 1.}});
    histos.add("hMultT0AperChannel", "", kTH2D, {axisMultChannelT0A, axisChannelsT0A});
    histos.add("hMultT0CperChannel", "", kTH2D, {axisMultChannelT0C, axisChannelsT0C});
  }

  void process(BCsRun3 const& bcs, aod::Zdcs const&, aod::FT0s const&)
  {
    int runNumber = bcs.iteratorAt(0).runNumber();
    const char* srun = Form("%d", runNumber);

    if (runNumber != lastRunNumber) {
      LOGP(info, "runNumber={}", runNumber);
      lastRunNumber = runNumber;
      int64_t ts = bcs.iteratorAt(0).timestamp();

      auto grplhcif = ccdb->getForTimeStamp<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", ts);
      beamPatternA = grplhcif->getBunchFilling().getBeamPattern(0);
      beamPatternC = grplhcif->getBunchFilling().getBeamPattern(1);
      bcPatternB = grplhcif->getBunchFilling().getBCPattern();
      bcPatternA = beamPatternA & ~beamPatternC;
      bcPatternC = ~beamPatternA & beamPatternC;

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
        float multZNA = bc.zdc().energyCommonZNA();
        float multZNC = bc.zdc().energyCommonZNC();

        histos.fill(HIST("hMultZNA"), multZNA);
        histos.fill(HIST("hMultZNC"), multZNC);
        histos.fill(HIST("hTimeZNA"), timeZNA);
        histos.fill(HIST("hTimeZNC"), timeZNC);
        if (bc.has_ft0() && TESTBIT(bc.ft0().triggerMask(), o2::ft0::Triggers::bitVertex)) { // TVX
          histos.fill(HIST("hTimeZNAselTVX"), timeZNA);
          histos.fill(HIST("hTimeZNCselTVX"), timeZNC);
        }

        if (bcPatternB[bc.globalBC() % nBCsPerOrbit]) { // B-mask
          histos.fill(HIST("hTimeZNAselB"), timeZNA);
          histos.fill(HIST("hTimeZNCselB"), timeZNC);
        }
        if (bcPatternA[bc.globalBC() % nBCsPerOrbit]) { // A-mask
          histos.fill(HIST("hTimeZNAselA"), timeZNA);
          histos.fill(HIST("hTimeZNCselA"), timeZNC);
        }
        if (bcPatternC[bc.globalBC() % nBCsPerOrbit]) { // C-mask
          histos.fill(HIST("hTimeZNAselC"), timeZNA);
          histos.fill(HIST("hTimeZNCselC"), timeZNC);
        }

        double meanTimeZNA = 0;
        double meanTimeZNC = 0;
        if (runNumber == 544795) {
          meanTimeZNA = 0.49;
          meanTimeZNC = -5.19;
        } else if (runNumber == 544911) {
          meanTimeZNA = -1.44;
          meanTimeZNC = -11.39;
        }

        if (fabs(timeZNA - meanTimeZNA) < 2) {
          histos.get<TH1>(HIST("hCounterZNA"))->Fill(srun, 1);
        }
        if (fabs(timeZNC - meanTimeZNC) < 2) {
          histos.get<TH1>(HIST("hCounterZNC"))->Fill(srun, 1);
        }
        if (fabs(timeZNA - meanTimeZNA) < 2 || fabs(timeZNC - meanTimeZNC) < 2) {
          histos.get<TH1>(HIST("hCounterZEM"))->Fill(srun, 1);
        }
      }

      if (!bc.has_ft0()) {
        continue;
      }
      for (unsigned int ic = 0; ic < bc.ft0().amplitudeA().size(); ic++) {
        histos.fill(HIST("hMultT0AperChannel"), bc.ft0().amplitudeA()[ic], bc.ft0().channelA()[ic]);
      }
      for (unsigned int ic = 0; ic < bc.ft0().amplitudeC().size(); ic++) {
        histos.fill(HIST("hMultT0CperChannel"), bc.ft0().amplitudeC()[ic], bc.ft0().channelC()[ic]);
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
