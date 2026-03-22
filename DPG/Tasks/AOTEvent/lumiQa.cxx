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

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/LHCConstants.h>
#include <DataFormatsFT0/Digit.h>
#include <DataFormatsParameters/AggregatedRunInfo.h>
#include <DataFormatsParameters/GRPLHCIFData.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TList.h>
#include <TString.h>

#include <Rtypes.h>

#include <bitset>
#include <cmath>
#include <cstdint>

using namespace o2;
using namespace o2::framework;
using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::Run3MatchedToBCSparse>;

struct LumiQaTask {
  Configurable<float> confTimeBinWidthInSec{"TimeBinWidthInSec", 60., "Width of time bins in seconds"}; // o2-linter: disable=name/configurable (temporary fix)
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  int lastRunNumber = -1;
  double maxSec = 1;
  double minSec = 0;
  TH1* hCalibV0A = nullptr;
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
    const AxisSpec axisMultV0A{1000, 0., 200000., "V0A multiplicity"};
    const AxisSpec axisMultT0M{1000, 0., 270000., "T0M multiplicity"};
    const AxisSpec axisMultT0A{1000, 0., 200000., "T0A multiplicity"};
    const AxisSpec axisMultT0C{1000, 0., 70000., "T0C multiplicity"};
    const AxisSpec axisCentV0A{100, 0., 100., "V0A centrality"};
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
    histos.add("hMultV0A", "", kTH1F, {axisMultV0A});
    histos.add("hMultT0M", "", kTH1F, {axisMultT0M});
    histos.add("hMultT0A", "", kTH1F, {axisMultT0A});
    histos.add("hMultT0C", "", kTH1F, {axisMultT0C});
    histos.add("hCentV0A", "", kTH1F, {axisCentV0A});
    histos.add("hCentT0C", "", kTH1F, {axisCentT0C});

    histos.add("hMultV0AselTVXB", "", kTH1F, {axisMultV0A});
    histos.add("hCentV0AselTVXB", "", kTH1F, {axisCentV0A});
    histos.add("hMultVCHselTVXB", "", kTH1F, {axisMultV0A});
    histos.add("hCentVCHselTVXB", "", kTH1F, {axisCentV0A});
    histos.add("hMultV0AselZACTVXB", "", kTH1F, {axisMultV0A});
    histos.add("hCentV0AselZACTVXB", "", kTH1F, {axisCentV0A});
    histos.add("hMultVCHselZACTVXB", "", kTH1F, {axisMultV0A});
    histos.add("hCentVCHselZACTVXB", "", kTH1F, {axisCentV0A});

    histos.add("hMultT0CselTVXB", "", kTH1F, {axisMultT0C});
    histos.add("hCentT0CselTVXB", "", kTH1F, {axisCentT0C});
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

    histos.add("hCounterTCEselB", "", kTH1D, {{1, 0., 1.}});
    histos.add("hCounterZNAselB", "", kTH1D, {{1, 0., 1.}});
    histos.add("hCounterZNCselB", "", kTH1D, {{1, 0., 1.}});
    histos.add("hCounterZEMselB", "", kTH1D, {{1, 0., 1.}});
    histos.add("hCounterVCHselB", "", kTH1D, {{1, 0., 1.}});

    histos.add("hCounterTCEselA", "", kTH1D, {{1, 0., 1.}});
    histos.add("hCounterZNAselA", "", kTH1D, {{1, 0., 1.}});
    histos.add("hCounterZNCselA", "", kTH1D, {{1, 0., 1.}});
    histos.add("hCounterZEMselA", "", kTH1D, {{1, 0., 1.}});
    histos.add("hCounterVCHselA", "", kTH1D, {{1, 0., 1.}});

    histos.add("hCounterTCEselC", "", kTH1D, {{1, 0., 1.}});
    histos.add("hCounterZNAselC", "", kTH1D, {{1, 0., 1.}});
    histos.add("hCounterZNCselC", "", kTH1D, {{1, 0., 1.}});
    histos.add("hCounterZEMselC", "", kTH1D, {{1, 0., 1.}});
    histos.add("hCounterVCHselC", "", kTH1D, {{1, 0., 1.}});

    histos.add("hMultT0AperChannel", "", kTH2D, {axisMultChannelT0A, axisChannelsT0A});
    histos.add("hMultT0CperChannel", "", kTH2D, {axisMultChannelT0C, axisChannelsT0C});
  }

  void process(BCsRun3 const& bcs, aod::Zdcs const&, aod::FT0s const&, aod::FV0As const&)
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

      hCalibV0A = reinterpret_cast<TH1*>(callst->FindObject("hCalibZeqFV0"));
      if (hCalibV0A == nullptr) {
        LOGF(info, "hCalibZeqFV0 histogram is not available for run=%d at timestamp=%llu", runNumber, ts);
        return;
      }

      hCalibT0C = reinterpret_cast<TH1*>(callst->FindObject("hCalibZeqFT0C"));
      if (hCalibT0C == nullptr) {
        LOGF(info, "hCalibZeqFT0C histogram is not available for run=%d at timestamp=%llu", runNumber, ts);
        return;
      }

      if (runNumber >= 500000) {
        auto runInfo = o2::parameters::AggregatedRunInfo::buildAggregatedRunInfo(o2::ccdb::BasicCCDBManager::instance(), runNumber);
        auto tsSOR = runInfo.sor;
        auto tsEOR = runInfo.eor;
        minSec = floor(tsSOR / 1000.);
        maxSec = ceil(tsEOR / 1000.);
      }

      int nTimeBins = static_cast<int>((maxSec - minSec) / confTimeBinWidthInSec);
      double timeInterval = nTimeBins * confTimeBinWidthInSec;

      const AxisSpec axisBCs{nBCsPerOrbit, 0., static_cast<double>(nBCsPerOrbit), ""};
      const AxisSpec axisSeconds{nTimeBins, 0, timeInterval, "seconds"};
      histos.add("hSecondsBcsTCEselB", "", kTH2D, {axisSeconds, axisBCs});
      histos.add("hSecondsBcsZNAselB", "", kTH2D, {axisSeconds, axisBCs});
      histos.add("hSecondsBcsZNCselB", "", kTH2D, {axisSeconds, axisBCs});
      histos.add("hSecondsBcsZEMselB", "", kTH2D, {axisSeconds, axisBCs});
      histos.add("hSecondsBcsVCHselB", "", kTH2D, {axisSeconds, axisBCs});
      histos.add("hSecondsBcsTCEselA", "", kTH2D, {axisSeconds, axisBCs});
      histos.add("hSecondsBcsZNAselA", "", kTH2D, {axisSeconds, axisBCs});
      histos.add("hSecondsBcsZNCselA", "", kTH2D, {axisSeconds, axisBCs});
      histos.add("hSecondsBcsZEMselA", "", kTH2D, {axisSeconds, axisBCs});
      histos.add("hSecondsBcsVCHselA", "", kTH2D, {axisSeconds, axisBCs});
      histos.add("hSecondsBcsTCEselC", "", kTH2D, {axisSeconds, axisBCs});
      histos.add("hSecondsBcsZNAselC", "", kTH2D, {axisSeconds, axisBCs});
      histos.add("hSecondsBcsZNCselC", "", kTH2D, {axisSeconds, axisBCs});
      histos.add("hSecondsBcsZEMselC", "", kTH2D, {axisSeconds, axisBCs});
      histos.add("hSecondsBcsVCHselC", "", kTH2D, {axisSeconds, axisBCs});
    }

    for (const auto& bc : bcs) {
      int64_t ts = bc.timestamp();
      auto inputMask = bc.inputMask();
      double secFromSOR = ts / 1000. - minSec;
      int bcInOrbit = bc.globalBC() % nBCsPerOrbit;
      bool maskB = bcPatternB[bcInOrbit];
      bool maskA = bcPatternA[bcInOrbit];
      bool maskC = bcPatternC[bcInOrbit];
      bool tvxCTP = TESTBIT(inputMask, 2);
      bool tscCTP = TESTBIT(inputMask, 3);
      bool tceCTP = TESTBIT(inputMask, 4);
      bool vchCTP = TESTBIT(inputMask, 9);
      bool zemCTP = TESTBIT(inputMask, 24);
      bool zncCTP = TESTBIT(inputMask, 25);

      LOGP(debug, "CTP: tvx={} tsc={} tce={} vch={} zem={} znc={}", tvxCTP, tscCTP, tceCTP, vchCTP, zemCTP, zncCTP);

      bool tvx = bc.has_ft0() ? TESTBIT(bc.ft0().triggerMask(), o2::ft0::Triggers::bitVertex) : 0;
      bool tsc = bc.has_ft0() ? TESTBIT(bc.ft0().triggerMask(), o2::ft0::Triggers::bitSCen) : 0;
      bool tce = bc.has_ft0() ? TESTBIT(bc.ft0().triggerMask(), o2::ft0::Triggers::bitCen) : 0;
      bool vch = bc.has_fv0a() ? TESTBIT(bc.fv0a().triggerMask(), o2::fit::Triggers::bitTrgCharge) : 0;

      float meanTimeZNA = 0.;
      float meanTimeZNC = 0.;
      bool zna = bc.has_zdc() ? std::fabs(bc.zdc().timeZNA() - meanTimeZNA) < 2 : 0; // ns
      bool znc = bc.has_zdc() ? std::fabs(bc.zdc().timeZNC() - meanTimeZNC) < 2 : 0; // ns
      bool zem = zna || znc;
      bool zac = zna && znc;

      // check if TCE triggers from FT0 mask and CTP are consistent
      if (tce != tceCTP) {
        LOGP(warning, "TCEfromFT0={} TCEfromCTP={}", tce, tceCTP);
      }

      if (bc.has_zdc()) {
        float timeZNA = bc.zdc().timeZNA();
        float timeZNC = bc.zdc().timeZNC();
        float multZNA = bc.zdc().energyCommonZNA();
        float multZNC = bc.zdc().energyCommonZNC();
        histos.fill(HIST("hMultZNA"), multZNA);
        histos.fill(HIST("hMultZNC"), multZNC);
        histos.fill(HIST("hTimeZNA"), timeZNA);
        histos.fill(HIST("hTimeZNC"), timeZNC);
        if (tvx) { // TVX
          histos.fill(HIST("hTimeZNAselTVX"), timeZNA);
          histos.fill(HIST("hTimeZNCselTVX"), timeZNC);
        }
        if (maskB) { // B-mask
          histos.fill(HIST("hTimeZNAselB"), timeZNA);
          histos.fill(HIST("hTimeZNCselB"), timeZNC);
        }
        if (maskA) { // A-mask
          histos.fill(HIST("hTimeZNAselA"), timeZNA);
          histos.fill(HIST("hTimeZNCselA"), timeZNC);
        }
        if (maskC) { // C-mask
          histos.fill(HIST("hTimeZNAselC"), timeZNA);
          histos.fill(HIST("hTimeZNCselC"), timeZNC);
        }

        // B-mask
        if (zna && maskB) {
          histos.get<TH1>(HIST("hCounterZNAselB"))->Fill(srun, 1);
          histos.fill(HIST("hSecondsBcsZNAselB"), secFromSOR, bcInOrbit);
        }
        if (znc && maskB) {
          histos.get<TH1>(HIST("hCounterZNCselB"))->Fill(srun, 1);
          histos.fill(HIST("hSecondsBcsZNCselB"), secFromSOR, bcInOrbit);
        }
        if (zem && maskB) {
          histos.get<TH1>(HIST("hCounterZEMselB"))->Fill(srun, 1);
          histos.fill(HIST("hSecondsBcsZEMselB"), secFromSOR, bcInOrbit);
        }

        // A-mask
        if (zna && maskA) {
          histos.get<TH1>(HIST("hCounterZNAselA"))->Fill(srun, 1);
          histos.fill(HIST("hSecondsBcsZNAselA"), secFromSOR, bcInOrbit);
        }
        if (znc && maskA) {
          histos.get<TH1>(HIST("hCounterZNCselA"))->Fill(srun, 1);
          histos.fill(HIST("hSecondsBcsZNCselA"), secFromSOR, bcInOrbit);
        }
        if (zem && maskA) {
          histos.get<TH1>(HIST("hCounterZEMselA"))->Fill(srun, 1);
          histos.fill(HIST("hSecondsBcsZEMselA"), secFromSOR, bcInOrbit);
        }

        // C-mask
        if (zna && maskC) {
          histos.get<TH1>(HIST("hCounterZNAselC"))->Fill(srun, 1);
          histos.fill(HIST("hSecondsBcsZNAselC"), secFromSOR, bcInOrbit);
        }
        if (znc && maskC) {
          histos.get<TH1>(HIST("hCounterZNCselC"))->Fill(srun, 1);
          histos.fill(HIST("hSecondsBcsZNCselC"), secFromSOR, bcInOrbit);
        }
        if (zem && maskC) {
          histos.get<TH1>(HIST("hCounterZEMselC"))->Fill(srun, 1);
          histos.fill(HIST("hSecondsBcsZEMselC"), secFromSOR, bcInOrbit);
        }
      }

      if (bc.has_fv0a()) {
        float multV0A = 0.;
        for (auto mult : bc.fv0a().amplitude()) {
          multV0A += mult;
        }
        float centV0A = hCalibV0A->GetBinContent(hCalibV0A->FindFixBin(multV0A));
        histos.fill(HIST("hMultV0A"), multV0A);
        histos.fill(HIST("hCentV0A"), centV0A);

        if (tvx && maskB) {
          histos.fill(HIST("hMultV0AselTVXB"), multV0A);
          histos.fill(HIST("hCentV0AselTVXB"), centV0A);
        }

        if (tvx && maskB && vch) {
          histos.fill(HIST("hMultVCHselTVXB"), multV0A);
          histos.fill(HIST("hCentVCHselTVXB"), centV0A);
          histos.get<TH1>(HIST("hCounterVCHselB"))->Fill(srun, 1);
          histos.fill(HIST("hSecondsBcsVCHselB"), secFromSOR, bcInOrbit);
        }

        if (tvx && maskA && vch) {
          histos.get<TH1>(HIST("hCounterVCHselA"))->Fill(srun, 1);
          histos.fill(HIST("hSecondsBcsVCHselA"), secFromSOR, bcInOrbit);
        }

        if (tvx && maskC && vch) {
          histos.get<TH1>(HIST("hCounterVCHselC"))->Fill(srun, 1);
          histos.fill(HIST("hSecondsBcsVCHselC"), secFromSOR, bcInOrbit);
        }

        if (tvx && maskB && zac) {
          histos.fill(HIST("hMultV0AselZACTVXB"), multV0A);
          histos.fill(HIST("hCentV0AselZACTVXB"), centV0A);
        }

        if (tvx && maskB && vch && zac) {
          histos.fill(HIST("hMultVCHselZACTVXB"), multV0A);
          histos.fill(HIST("hCentVCHselZACTVXB"), centV0A);
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

      if (tvx && maskB) {
        histos.fill(HIST("hMultT0CselTVXB"), multT0C);
        histos.fill(HIST("hCentT0CselTVXB"), centT0C);
      }

      if (tsc) {
        histos.fill(HIST("hMultT0MselTSC"), multT0M);
      }

      if (tce) {
        histos.fill(HIST("hMultT0MselTCE"), multT0M);
        histos.fill(HIST("hMultT0CselTCE"), multT0C);
        histos.fill(HIST("hCentT0CselTCE"), centT0C);
      }

      if (tce && tvx) {
        histos.fill(HIST("hMultT0CselTVXTCE"), multT0C);
        histos.fill(HIST("hCentT0CselTVXTCE"), centT0C);
      }

      if (tce && tvx && maskB) {
        histos.fill(HIST("hMultT0CselTVXTCEB"), multT0C);
        histos.fill(HIST("hCentT0CselTVXTCEB"), centT0C);
        histos.get<TH1>(HIST("hCounterTCEselB"))->Fill(srun, 1);
        histos.fill(HIST("hSecondsBcsTCEselB"), secFromSOR, bcInOrbit);
      }
    } // bc loop
  } // process
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<LumiQaTask>(cfgc)};
}
