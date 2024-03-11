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
  Configurable<int> fMaxBcRange{"maxBcRange", 50, "Max bc range for veto analysis"};

  int32_t fRun{0};

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

  void init(InitContext&)
  {
    const AxisSpec axisAmp{101, 0., 1000., ""};
    const AxisSpec axisBCs{51, 0., 50., ""};
    const AxisSpec axisT{101, -10., 10., ""};
    const AxisSpec axisCounters{10, 0., 10., ""};

    hr.add("hTimeZNA", "", kTH1F, {axisT});
    hr.add("hTimeZNC", "", kTH1F, {axisT});

    hr.add("hCountBCsZDC", "", kTH1F, {axisCounters});
    auto hCountBCsZDC = hr.get<TH1>(HIST("hCountBCsZDC"));
    hCountBCsZDC->GetXaxis()->SetBinLabel(1, "Total");
    hCountBCsZDC->GetXaxis()->SetBinLabel(2, "0n0n");
    hCountBCsZDC->GetXaxis()->SetBinLabel(3, "0nXn");
    hCountBCsZDC->GetXaxis()->SetBinLabel(4, "Xn0n");
    hCountBCsZDC->GetXaxis()->SetBinLabel(5, "XnXn");
    hCountBCsZDC->GetXaxis()->SetBinLabel(6, "ZNA");
    hCountBCsZDC->GetXaxis()->SetBinLabel(7, "ZNC");

    hr.add("hSelBCAmpV0A_total", "", kTH2F, {axisAmp, axisBCs});
    hr.add("hSelBCAmpV0A_0n0n", "", kTH2F, {axisAmp, axisBCs});
    hr.add("hSelBCAmpV0A_0nXn", "", kTH2F, {axisAmp, axisBCs});
    hr.add("hSelBCAmpV0A_Xn0n", "", kTH2F, {axisAmp, axisBCs});
    hr.add("hSelBCAmpV0A_XnXn", "", kTH2F, {axisAmp, axisBCs});

    hr.add("hSelBCAmpT0A_total", "", kTH2F, {axisAmp, axisBCs});
    hr.add("hSelBCAmpT0A_0n0n", "", kTH2F, {axisAmp, axisBCs});
    hr.add("hSelBCAmpT0A_0nXn", "", kTH2F, {axisAmp, axisBCs});
    hr.add("hSelBCAmpT0A_Xn0n", "", kTH2F, {axisAmp, axisBCs});
    hr.add("hSelBCAmpT0A_XnXn", "", kTH2F, {axisAmp, axisBCs});
  }

  void process(const o2::aod::BCs& bcs,
               const o2::aod::FT0s& ft0s,
               const o2::aod::FV0As& fv0s,
               const o2::aod::Zdcs& zdcs)
  {
    const o2::aod::BCs::iterator& fbc = bcs.begin();
    if (fbc.runNumber() != fRun)
      getBCPattern(fbc);

    auto nBCs = bcs.size();

    std::vector<int64_t> gbcs(nBCs, -1);

    int64_t nBCsPat_total = 0;
    int64_t nBCsPat_0n0n = 0;
    int64_t nBCsPat_0nXn = 0;
    int64_t nBCsPat_Xn0n = 0;
    int64_t nBCsPat_XnXn = 0;
    int64_t nBCsPat_ZNA = 0;
    int64_t nBCsPat_ZNC = 0;

    int nBcs = bcs.size();
    int nFt0s = ft0s.size();
    int nFv0s = fv0s.size();
    int nZdcs = zdcs.size();

    for (auto i = 0; i < nBcs; ++i) {
      const auto& bc = bcs.iteratorAt(i);
      if (bcPatternB[bc.globalBC() % 3564] == 0)
        continue;
      gbcs[i] = bc.globalBC();
      nBCsPat_total++;
    }

    std::map<int64_t, float> mapGlobalBcWithZNA{};
    std::map<int64_t, float> mapGlobalBcWithZNC{};
    for (auto i = 0; i < nZdcs; ++i) {
      const auto& zdc = zdcs.iteratorAt(i);
      auto bcId = zdc.bcId();
      int64_t globalBC = gbcs[bcId];
      if (globalBC == -1)
        continue;
      float timeA = zdc.timeZNA();
      float timeC = zdc.timeZNC();
      hr.get<TH1>(HIST("hTimeZNA"))->Fill(timeA);
      hr.get<TH1>(HIST("hTimeZNC"))->Fill(timeC);
      bool hasZNA = !(std::abs(timeA) > 2.f);
      bool hasZNC = !(std::abs(timeC) > 2.f);
      mapGlobalBcWithZNA[globalBC] = timeA;
      mapGlobalBcWithZNC[globalBC] = timeC;
      if (!hasZNA && !hasZNC)
        nBCsPat_0n0n++;
      if (!hasZNA && hasZNC)
        nBCsPat_0nXn++;
      if (hasZNA && !hasZNC)
        nBCsPat_Xn0n++;
      if (hasZNA && hasZNC)
        nBCsPat_XnXn++;
      if (hasZNA)
        nBCsPat_ZNA++;
      if (hasZNC)
        nBCsPat_ZNC++;
    }

    // looking for bcs without any zdc signals --> adding to 0n0n
    for (auto i = 0; i < nBcs; ++i) {
      auto globalBC = gbcs[i];
      if (globalBC == -1)
        continue;
      auto itA = mapGlobalBcWithZNA.find(globalBC);
      auto itC = mapGlobalBcWithZNC.find(globalBC);
      if (itA == mapGlobalBcWithZNA.end() && itC == mapGlobalBcWithZNC.end())
        nBCsPat_0n0n++;
    }

    hr.get<TH1>(HIST("hCountBCsZDC"))->Fill("Total", static_cast<double>(nBCsPat_total));
    hr.get<TH1>(HIST("hCountBCsZDC"))->Fill("0n0n", static_cast<double>(nBCsPat_0n0n));
    hr.get<TH1>(HIST("hCountBCsZDC"))->Fill("0nXn", static_cast<double>(nBCsPat_0nXn));
    hr.get<TH1>(HIST("hCountBCsZDC"))->Fill("Xn0n", static_cast<double>(nBCsPat_Xn0n));
    hr.get<TH1>(HIST("hCountBCsZDC"))->Fill("XnXn", static_cast<double>(nBCsPat_XnXn));
    hr.get<TH1>(HIST("hCountBCsZDC"))->Fill("ZNA", static_cast<double>(nBCsPat_ZNA));
    hr.get<TH1>(HIST("hCountBCsZDC"))->Fill("ZNC", static_cast<double>(nBCsPat_ZNC));

    std::map<int64_t, float> mapGlobalBcWithV0A{};
    for (auto i = 0; i < nFv0s; ++i) {
      const auto& fv0 = fv0s.iteratorAt(i);
      if (!TESTBIT(fv0.triggerMask(), o2::fit::Triggers::bitA))
        continue;
      auto bcId = fv0.bcId();
      int64_t globalBC = gbcs[bcId];
      if (globalBC == -1)
        continue;
      if (std::abs(fv0.time()) > 15.f)
        continue;
      //
      const auto& amps = fv0.amplitude();
      auto amp = std::accumulate(amps.begin(), amps.end(), 0.f);
      mapGlobalBcWithV0A[globalBC] = amp;
    }

    std::map<int64_t, float> mapGlobalBcWithT0A{};
    for (auto i = 0; i < nFt0s; ++i) {
      const auto& ft0 = ft0s.iteratorAt(i);
      if (!TESTBIT(ft0.triggerMask(), o2::fit::Triggers::bitVertex))
        continue;
      auto bcId = ft0.bcId();
      int64_t globalBC = gbcs[bcId];
      if (globalBC == -1)
        continue;
      if (std::abs(ft0.timeA()) > 2.f)
        continue;
      //
      const auto& amps = ft0.amplitudeA();
      auto amp = std::accumulate(amps.begin(), amps.end(), 0.f);
      mapGlobalBcWithT0A[globalBC] = amp;
    }

    hr.get<TH2>(HIST("hSelBCAmpV0A_total"))->Fill(0.f, 0.f, static_cast<double>(nBCsPat_total));
    hr.get<TH2>(HIST("hSelBCAmpV0A_0n0n"))->Fill(0.f, 0.f, static_cast<double>(nBCsPat_0n0n));
    hr.get<TH2>(HIST("hSelBCAmpV0A_0nXn"))->Fill(0.f, 0.f, static_cast<double>(nBCsPat_0nXn));
    hr.get<TH2>(HIST("hSelBCAmpV0A_Xn0n"))->Fill(0.f, 0.f, static_cast<double>(nBCsPat_Xn0n));
    hr.get<TH2>(HIST("hSelBCAmpV0A_XnXn"))->Fill(0.f, 0.f, static_cast<double>(nBCsPat_XnXn));

    hr.get<TH2>(HIST("hSelBCAmpT0A_total"))->Fill(0.f, 0.f, static_cast<double>(nBCsPat_total));
    hr.get<TH2>(HIST("hSelBCAmpT0A_0n0n"))->Fill(0.f, 0.f, static_cast<double>(nBCsPat_0n0n));
    hr.get<TH2>(HIST("hSelBCAmpT0A_0nXn"))->Fill(0.f, 0.f, static_cast<double>(nBCsPat_0nXn));
    hr.get<TH2>(HIST("hSelBCAmpT0A_Xn0n"))->Fill(0.f, 0.f, static_cast<double>(nBCsPat_Xn0n));
    hr.get<TH2>(HIST("hSelBCAmpT0A_XnXn"))->Fill(0.f, 0.f, static_cast<double>(nBCsPat_XnXn));

    for (auto gbc : gbcs) {
      if (gbc == -1)
        continue;
      auto itZNA = mapGlobalBcWithZNA.find(gbc);
      auto itZNC = mapGlobalBcWithZNC.find(gbc);
      float timeZNA = -999.f;
      float timeZNC = -999.f;
      if (itZNA != mapGlobalBcWithZNA.end())
        timeZNA = itZNA->second;
      if (itZNC != mapGlobalBcWithZNC.end())
        timeZNC = itZNC->second;
      bool hasZNA = !(std::abs(timeZNA) > 2.f);
      bool hasZNC = !(std::abs(timeZNC) > 2.f);
      for (auto r = 1; r < fMaxBcRange; ++r) {
        auto maxV0A = -999.f;
        auto s = gbc - r;
        auto e = gbc + r;
        auto it = mapGlobalBcWithV0A.lower_bound(s);
        while (it->first <= e && it != mapGlobalBcWithV0A.end()) {
          auto a = it->second;
          if (a > maxV0A)
            maxV0A = a;
          ++it;
        }
        for (auto iamp = 1; iamp <= 100; ++iamp) {
          float amp = 10.f * iamp;
          if (maxV0A < amp) {
            hr.get<TH2>(HIST("hSelBCAmpV0A_total"))->Fill(amp, r, 1);
            if (!hasZNA && !hasZNC)
              hr.get<TH2>(HIST("hSelBCAmpV0A_0n0n"))->Fill(amp, r, 1);
            if (!hasZNA && hasZNC)
              hr.get<TH2>(HIST("hSelBCAmpV0A_0nXn"))->Fill(amp, r, 1);
            if (hasZNA && !hasZNC)
              hr.get<TH2>(HIST("hSelBCAmpV0A_Xn0n"))->Fill(amp, r, 1);
            if (hasZNA && hasZNC)
              hr.get<TH2>(HIST("hSelBCAmpV0A_XnXn"))->Fill(amp, r, 1);
          }
        }
      }
      for (auto r = 1; r < fMaxBcRange; ++r) {
        auto maxT0A = -999.f;
        auto s = gbc - r;
        auto e = gbc + r;
        auto it = mapGlobalBcWithT0A.lower_bound(s);
        while (it->first <= e && it != mapGlobalBcWithT0A.end()) {
          auto a = it->second;
          if (a > maxT0A)
            maxT0A = a;
          ++it;
        }
        for (auto iamp = 1; iamp <= 100; ++iamp) {
          float amp = 10.f * iamp;
          if (maxT0A < amp)
            hr.get<TH2>(HIST("hSelBCAmpT0A_total"))->Fill(amp, r, 1);
          if (!hasZNA && !hasZNC)
            hr.get<TH2>(HIST("hSelBCAmpT0A_0n0n"))->Fill(amp, r, 1);
          if (!hasZNA && hasZNC)
            hr.get<TH2>(HIST("hSelBCAmpT0A_0nXn"))->Fill(amp, r, 1);
          if (hasZNA && !hasZNC)
            hr.get<TH2>(HIST("hSelBCAmpT0A_Xn0n"))->Fill(amp, r, 1);
          if (hasZNA && hasZNC)
            hr.get<TH2>(HIST("hSelBCAmpT0A_XnXn"))->Fill(amp, r, 1);
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
