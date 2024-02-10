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
#include "CommonConstants/LHCConstants.h"
#include "DataFormatsFIT/Triggers.h"
#include "DataFormatsParameters/GRPLHCIFData.h"

using namespace o2::framework;
using namespace o2::framework::expressions;

struct UpcVetoAnalysis {
  int32_t fRun{0};

  float ampLo{0.f};
  float ampHi{20000.f};
  float damp{100.f};
  int32_t nAmpThr;
  int32_t nBCRanges{5};

  std::bitset<o2::constants::lhc::LHCMaxBunches> bcPatternB;

  HistogramRegistry hr{"hr", {}, OutputObjHandlingPolicy::AnalysisObject};

  void getBCPattern(const o2::aod::BCs::iterator& bc)
  {
    fRun = bc.runNumber();
    o2::ccdb::CcdbApi ccdb_api;
    ccdb_api.init("http://alice-ccdb.cern.ch");
    std::map<string, string> metadata, headers;
    headers = ccdb_api.retrieveHeaders(Form("RCT/Info/RunInformation/%i", fRun), metadata, -1);
    int64_t ts = std::atol(headers["SOR"].c_str());
    LOGP(info, "ts={}", ts);
    auto grplhcif = ccdb_api.retrieveFromTFileAny<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", metadata, ts);
    bcPatternB = grplhcif->getBunchFilling().getBCPattern();
  }

  auto findClosestBCIter(int64_t globalBC, std::vector<int64_t>& bcs)
  {
    auto it = std::lower_bound(bcs.begin(), bcs.end(), globalBC);
    auto bc1 = *it;
    auto it1 = it;
    if (it != bcs.begin())
      --it;
    auto it2 = it;
    auto bc2 = *it;
    auto dbc1 = bc1 >= globalBC ? bc1 - globalBC : globalBC - bc1;
    auto dbc2 = bc2 >= globalBC ? bc2 - globalBC : globalBC - bc2;
    return (dbc1 <= dbc2) ? it1 : it2;
  }

  template <typename T1, typename T2>
  T2 findMaxInRange(int32_t range,
                    int64_t cent,
                    const typename std::vector<T1>::iterator startIt,
                    std::vector<T1>& vPosi,
                    std::vector<T2>& vVals)
  {
    auto pos = cent;
    auto it = startIt;
    T2 max = -999;
    while (pos - cent < range && it != vPosi.end()) {
      auto idx = it - vPosi.begin();
      if (vVals[idx] > max)
        max = vVals[idx];
      ++it;
      pos = *it;
    }
    it = startIt;
    pos = cent;
    while (cent - pos < range && it >= vPosi.begin()) {
      --it;
      pos = *it;
      auto idx = it - vPosi.begin();
      if (vVals[idx] > max)
        max = vVals[idx];
    }
    return max;
  }

  void init(InitContext&)
  {
    const AxisSpec axisTrigCounter{1, 0., 1., ""};
    hr.add("hCounterTCE", "", kTH1F, {axisTrigCounter});

    nAmpThr = (ampHi - ampLo) / damp;
    const AxisSpec axisCounters{100, 0.5, 100.5, ""};
    const AxisSpec axisAmp{nAmpThr, ampLo, ampHi, ""};
    const AxisSpec axisAmp2{1000, ampLo, ampHi * 10.f, ""};
    const AxisSpec axisBCRanges{nBCRanges, 10., 60.};
    hr.add("hSelBCsTOR", "", kTH1F, {axisCounters});
    hr.add("hSelBCsTSC", "", kTH1F, {axisCounters});
    hr.add("hSelBCsTVX", "", kTH1F, {axisCounters});
    hr.add("hSelBCsV0A", "", kTH1F, {axisCounters});
    hr.add("hSelBCAmpT0M", "", kTH2F, {axisAmp, axisBCRanges});
    hr.add("hSelBCAmpV0A", "", kTH2F, {axisAmp, axisBCRanges});
    hr.add("hAmpT0M", "", kTH1F, {axisAmp});
    hr.add("hAmpV0A", "", kTH1F, {axisAmp});
  }

  void process(const o2::aod::BCs& bcs,
               const o2::aod::FT0s& ft0s,
               const o2::aod::FV0As& fv0as)
  {
    const o2::aod::BCs::iterator& fbc = bcs.begin();
    if (fbc.runNumber() != fRun)
      getBCPattern(fbc);

    auto nBCs = bcs.size();

    std::vector<int64_t> gbcs(nBCs, -1);
    std::vector<int64_t> gbcsWithTOR;
    std::vector<int64_t> gbcsWithT0;
    std::vector<int64_t> gbcsWithTSC;
    std::vector<int64_t> gbcsWithTVX;
    std::vector<int64_t> gbcsWithV0A;
    std::vector<float> ampsT0M;
    std::vector<float> ampsV0A;

    int64_t nBCsPat = 0;
    for (const o2::aod::BCs::iterator& bc : bcs) {
      if (bcPatternB[bc.globalBC() % 3564] == 0)
        continue;
      gbcs[bc.globalIndex()] = bc.globalBC();
      nBCsPat++;
      hr.get<TH1>(HIST("hSelBCsV0A"))->AddBinContent(1, 1);
      hr.get<TH1>(HIST("hSelBCsTOR"))->AddBinContent(1, 1);
      hr.get<TH1>(HIST("hSelBCsTSC"))->AddBinContent(1, 1);
      hr.get<TH1>(HIST("hSelBCsTVX"))->AddBinContent(1, 1);
    }

    for (int32_t i = 1; i <= nBCRanges; ++i) {
      auto bcRange = i * 10;
      hr.get<TH2>(HIST("hSelBCAmpT0M"))->Fill(0.f, bcRange, nBCsPat);
      hr.get<TH2>(HIST("hSelBCAmpV0A"))->Fill(0.f, bcRange, nBCsPat);
    }

    for (const o2::aod::FV0As::iterator& fv0a : fv0as) {
      int64_t globalBC = gbcs[fv0a.bc().globalIndex()];
      if (globalBC == -1)
        continue;
      if (std::abs(fv0a.time()) > 15.f)
        continue;
      gbcsWithV0A.push_back(globalBC);
      const auto& amps = fv0a.amplitude();
      auto amp = std::accumulate(amps.begin(), amps.end(), 0.f);
      ampsV0A.push_back(amp);
      hr.get<TH1>(HIST("hAmpV0A"))->Fill(amp);
    }

    for (const o2::aod::FT0s::iterator& ft0 : ft0s) {
      if (TESTBIT(ft0.triggerMask(), o2::fit::Triggers::bitVertex) &&
          TESTBIT(ft0.triggerMask(), o2::fit::Triggers::bitCen))
        hr.get<TH1>(HIST("hCounterTCE"))->Fill(0);
      int64_t globalBC = gbcs[ft0.bc().globalIndex()];
      if (globalBC == -1)
        continue;
      gbcsWithT0.push_back(globalBC);
      const auto& ampsA = ft0.amplitudeA();
      const auto& ampsC = ft0.amplitudeC();
      auto ampA = std::accumulate(ampsA.begin(), ampsA.end(), 0.f);
      auto ampC = std::accumulate(ampsC.begin(), ampsC.end(), 0.f);
      ampsT0M.push_back(ampA + ampC);
      hr.get<TH1>(HIST("hAmpT0M"))->Fill(ampA + ampC);
      if (!(std::abs(ft0.timeA()) > 2.f && std::abs(ft0.timeC()) > 2.f)) // TOR
        gbcsWithTOR.push_back(globalBC);
      if (TESTBIT(ft0.triggerMask(), o2::fit::Triggers::bitVertex)) // TVX
        gbcsWithTVX.push_back(globalBC);
      if (TESTBIT(ft0.triggerMask(), o2::fit::Triggers::bitVertex) &&
          (TESTBIT(ft0.triggerMask(), o2::fit::Triggers::bitCen) ||
           TESTBIT(ft0.triggerMask(), o2::fit::Triggers::bitSCen))) // TVX & (TSC | TCE)
        gbcsWithTSC.push_back(globalBC);
    }

    for (auto gbc : gbcs) {
      if (gbc == -1)
        continue;
      auto bcV0Ait = findClosestBCIter(gbc, gbcsWithV0A);
      auto bcTORit = findClosestBCIter(gbc, gbcsWithTOR);
      auto bcT0Mit = findClosestBCIter(gbc, gbcsWithT0);
      auto bcTSCit = findClosestBCIter(gbc, gbcsWithTSC);
      auto bcTVXit = findClosestBCIter(gbc, gbcsWithTVX);
      auto bcV0A = bcV0Ait != gbcsWithV0A.end() ? *bcV0Ait : gbc + 999;
      auto bcTOR = bcTORit != gbcsWithTOR.end() ? *bcTORit : gbc + 999;
      auto bcTSC = bcTSCit != gbcsWithTSC.end() ? *bcTSCit : gbc + 999;
      auto bcTVX = bcTVXit != gbcsWithTVX.end() ? *bcTVXit : gbc + 999;
      auto bcDistV0A = std::abs(gbc - bcV0A);
      auto bcDistTOR = std::abs(gbc - bcTOR);
      auto bcDistTSC = std::abs(gbc - bcTSC);
      auto bcDistTVX = std::abs(gbc - bcTVX);
      for (int32_t vetoBC = 1; vetoBC <= 99; ++vetoBC) {
        if (bcDistV0A < vetoBC)
          hr.get<TH1>(HIST("hSelBCsV0A"))->AddBinContent(vetoBC + 1, 1);
        if (bcDistTOR < vetoBC)
          hr.get<TH1>(HIST("hSelBCsTOR"))->AddBinContent(vetoBC + 1, 1);
        if (bcDistTSC < vetoBC)
          hr.get<TH1>(HIST("hSelBCsTSC"))->AddBinContent(vetoBC + 1, 1);
        if (bcDistTVX < vetoBC)
          hr.get<TH1>(HIST("hSelBCsTVX"))->AddBinContent(vetoBC + 1, 1);
      }
      for (int32_t iBCRange = 1; iBCRange <= 5; ++iBCRange) {
        auto bcRange = 5 + iBCRange * 10;
        auto maxAmpV0A = findMaxInRange(bcRange, gbc, bcV0Ait, gbcsWithV0A, ampsV0A);
        auto maxAmpT0M = findMaxInRange(bcRange, gbc, bcT0Mit, gbcsWithT0, ampsT0M);
        for (int32_t iamp = 1; iamp <= nAmpThr; ++iamp) {
          auto amp = damp * iamp;
          if (maxAmpV0A < amp)
            hr.get<TH2>(HIST("hSelBCAmpT0M"))->Fill(amp, bcRange + 1, 1);
          if (maxAmpT0M < amp)
            hr.get<TH2>(HIST("hSelBCAmpV0A"))->Fill(amp, bcRange + 1, 1);
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UpcVetoAnalysis>(cfgc)};
}
