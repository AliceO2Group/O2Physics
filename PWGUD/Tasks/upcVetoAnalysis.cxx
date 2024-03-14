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
#include "DataFormatsParameters/GRPECSObject.h"

using namespace o2::framework;
using namespace o2::framework::expressions;

struct UpcVetoAnalysis {
  // amplitude axis
  int nbAmp = 101;
  double minAmp = 0.;
  double maxAmp = 1000.;

  // bcs axis
  int nbBcs = 51;
  double minBcs = 0.;
  double maxBcs = 50.;

  // zdc time axis
  int nbT = 101;
  double minT = -10.;
  double maxT = 10.;

  HistogramRegistry hr{"hr", {}, OutputObjHandlingPolicy::AnalysisObject};

  OutputObj<TH1D> hCountBCs{TH1D("hCountBCs", "", 10, 0., 10.)};

  OutputObj<TH1D> hTimeZNA{TH1D("hTimeZNA", "", nbT, minT, maxT)};
  OutputObj<TH1D> hTimeZNC{TH1D("hTimeZNC", "", nbT, minT, maxT)};

  OutputObj<TH2D> hSelBCAmpV0A_Total{TH2D("hSelBCAmpV0A_Total", "", nbAmp, minAmp, maxAmp, nbBcs, minBcs, maxBcs)};
  OutputObj<TH2D> hSelBCAmpV0A_0n0n{TH2D("hSelBCAmpV0A_0n0n", "", nbAmp, minAmp, maxAmp, nbBcs, minBcs, maxBcs)};
  OutputObj<TH2D> hSelBCAmpV0A_0nXn{TH2D("hSelBCAmpV0A_0nXn", "", nbAmp, minAmp, maxAmp, nbBcs, minBcs, maxBcs)};
  OutputObj<TH2D> hSelBCAmpV0A_Xn0n{TH2D("hSelBCAmpV0A_Xn0n", "", nbAmp, minAmp, maxAmp, nbBcs, minBcs, maxBcs)};
  OutputObj<TH2D> hSelBCAmpV0A_XnXn{TH2D("hSelBCAmpV0A_XnXn", "", nbAmp, minAmp, maxAmp, nbBcs, minBcs, maxBcs)};

  OutputObj<TH2D> hSelBCAmpT0A_Total{TH2D("hSelBCAmpT0A_Total", "", nbAmp, minAmp, maxAmp, nbBcs, minBcs, maxBcs)};
  OutputObj<TH2D> hSelBCAmpT0A_0n0n{TH2D("hSelBCAmpT0A_0n0n", "", nbAmp, minAmp, maxAmp, nbBcs, minBcs, maxBcs)};
  OutputObj<TH2D> hSelBCAmpT0A_0nXn{TH2D("hSelBCAmpT0A_0nXn", "", nbAmp, minAmp, maxAmp, nbBcs, minBcs, maxBcs)};
  OutputObj<TH2D> hSelBCAmpT0A_Xn0n{TH2D("hSelBCAmpT0A_Xn0n", "", nbAmp, minAmp, maxAmp, nbBcs, minBcs, maxBcs)};
  OutputObj<TH2D> hSelBCAmpT0A_XnXn{TH2D("hSelBCAmpT0A_XnXn", "", nbAmp, minAmp, maxAmp, nbBcs, minBcs, maxBcs)};

  TH1* hOrbitTotal;
  TH1* hOrbit0n0n;
  TH1* hOrbit0nXn;
  TH1* hOrbitXn0n;
  TH1* hOrbitXnXn;
  TH1* hOrbitZNA;
  TH1* hOrbitZNC;
  TH1* hOrbitTCE;
  TH1* hOrbitV0A;
  TH1* hOrbitT0A;

  TH1* hSecondsTotal;
  TH1* hSeconds0n0n;
  TH1* hSeconds0nXn;
  TH1* hSecondsXn0n;
  TH1* hSecondsXnXn;
  TH1* hSecondsZNA;
  TH1* hSecondsZNC;
  TH1* hSecondsTCE;
  TH1* hSecondsV0A;
  TH1* hSecondsT0A;

  int32_t fRun{0};

  int64_t fTsOrbitReset{0};
  int64_t fTsSOR{0};
  int64_t fTsEOR{0};
  int64_t fBcSOR{0};
  int64_t fMinOrbit{0};
  int32_t fNOrbits{0};
  int64_t fNOrbitsPerTF{0};
  int64_t fNBcsPerTF{0};

  std::bitset<o2::constants::lhc::LHCMaxBunches> bcPatternB{};
  int32_t fNBcsB{0};

  void init(InitContext&)
  {
    int il = 1;
    for (const auto* lbl : {"All", "Total", "0n0n", "0nXn", "Xn0n", "XnXn", "ZNA", "ZNC", "TCE"}) {
      hCountBCs->GetXaxis()->SetBinLabel(il++, lbl);
    }
  }

  void getBcInfo(const o2::aod::BCs::iterator& bc)
  {
    fRun = bc.runNumber();
    o2::ccdb::CcdbApi ccdb_api;
    ccdb_api.init("http://alice-ccdb.cern.ch");
    std::map<std::string, std::string> metadata;
    std::map<std::string, std::string> headers = ccdb_api.retrieveHeaders(Form("RCT/Info/RunInformation/%i", fRun),
                                                                          metadata, -1);

    fTsSOR = std::atol(headers["SOR"].c_str());
    fTsEOR = std::atol(headers["EOR"].c_str());
    LOGP(info, "fTsSOR={}, fTsEOR={}", fTsSOR, fTsEOR);

    auto grplhcif = ccdb_api.retrieveFromTFileAny<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", metadata,
                                                                                fTsSOR);
    bcPatternB = grplhcif->getBunchFilling().getBCPattern();
    fNBcsB = 0;
    for (auto i = 0; i < o2::constants::lhc::LHCMaxBunches; ++i) {
      if (bcPatternB[i])
        fNBcsB++;
    }
    LOGP(info, "fNBcsB={}", fNBcsB);

    auto ctpx = ccdb_api.retrieveFromTFileAny<std::vector<Long64_t>>("CTP/Calib/OrbitReset", metadata, fTsSOR);
    fTsOrbitReset = (*ctpx)[0]; // us
    LOGP(info, "fTsOrbitReset={} us", fTsOrbitReset);

    auto grpecs = ccdb_api.retrieveFromTFileAny<o2::parameters::GRPECSObject>("GLO/Config/GRPECS", metadata, fTsSOR);
    fNOrbitsPerTF = grpecs->getNHBFPerTF(); // assuming 1 orbit = 1 HBF

    int64_t orbitSOR = (fTsSOR * 1000 - fTsOrbitReset) / o2::constants::lhc::LHCOrbitMUS;
    int64_t orbitEOR = (fTsEOR * 1000 - fTsOrbitReset) / o2::constants::lhc::LHCOrbitMUS;

    orbitSOR = orbitSOR / fNOrbitsPerTF * fNOrbitsPerTF;
    orbitEOR = orbitEOR / fNOrbitsPerTF * fNOrbitsPerTF;

    fMinOrbit = orbitSOR;
    fNOrbits = orbitEOR - orbitSOR;

    fBcSOR = orbitSOR * o2::constants::lhc::LHCMaxBunches;
    fNBcsPerTF = fNOrbitsPerTF * o2::constants::lhc::LHCMaxBunches;
  }

  void collectVeto(const std::vector<int64_t>& gbcs,
                   const std::vector<int64_t>& gbcsOrbit,
                   const std::vector<int64_t>& gbcsTs,
                   const char* det,
                   const std::map<int64_t, float>& mapGlobalBcWithFit,
                   const std::map<int64_t, float>& mapGlobalBcWithZNA,
                   const std::map<int64_t, float>& mapGlobalBcWithZNC)
  {
    TH1* hOrbit = nullptr;
    TH1* hSeconds = nullptr;
    TH1* hTotal = nullptr;
    TH1* h0n0n = nullptr;
    TH1* h0nXn = nullptr;
    TH1* hXn0n = nullptr;
    TH1* hXnXn = nullptr;
    if (std::strcmp(det, "V0A") == 0) {
      hOrbit = hOrbitV0A;
      hSeconds = hSecondsV0A;
      hTotal = hSelBCAmpV0A_Total.object.get();
      h0n0n = hSelBCAmpV0A_0n0n.object.get();
      h0nXn = hSelBCAmpV0A_0nXn.object.get();
      hXn0n = hSelBCAmpV0A_Xn0n.object.get();
      hXnXn = hSelBCAmpV0A_XnXn.object.get();
    }
    if (std::strcmp(det, "T0A") == 0) {
      hOrbit = hOrbitT0A;
      hSeconds = hSecondsT0A;
      hTotal = hSelBCAmpT0A_Total.object.get();
      h0n0n = hSelBCAmpT0A_0n0n.object.get();
      h0nXn = hSelBCAmpT0A_0nXn.object.get();
      hXn0n = hSelBCAmpT0A_Xn0n.object.get();
      hXnXn = hSelBCAmpT0A_XnXn.object.get();
    }
    for (auto i = 0; i < gbcs.size(); ++i) {
      auto gbc = gbcs[i];
      if (gbc == -1)
        continue;
      auto orbit = gbcsOrbit[i];
      auto ts = gbcsTs[i];
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
      for (auto r : {1, 5, 10}) {
        auto maxAmp = -999.f;
        auto s = gbc - r;
        auto e = gbc + r;
        auto it = mapGlobalBcWithFit.lower_bound(s);
        while (it->first <= e && it != mapGlobalBcWithFit.end()) {
          auto a = it->second;
          if (a > maxAmp)
            maxAmp = a;
          ++it;
        }
        for (auto iamp = 1; iamp <= 100; ++iamp) {
          float amp = 10.f * iamp;
          if (maxAmp < amp) {
            auto bin = hTotal->FindFixBin(amp, r);
            hTotal->AddBinContent(bin, 1);
            if (!hasZNA && !hasZNC)
              h0n0n->AddBinContent(bin, 1);
            if (!hasZNA && hasZNC)
              h0nXn->AddBinContent(bin, 1);
            if (hasZNA && !hasZNC)
              hXn0n->AddBinContent(bin, 1);
            if (hasZNA && hasZNC)
              hXnXn->AddBinContent(bin, 1);
          }
        }
        if (r == 1 && maxAmp < 100.f) {
          hOrbit->Fill(orbit - fMinOrbit);
          hSeconds->Fill(ts);
        }
      }
    }
  }

  void process(const o2::aod::BCs& bcs,
               const o2::aod::FT0s& ft0s,
               const o2::aod::FV0As& fv0s,
               const o2::aod::Zdcs& zdcs)
  {
    int nBcs = bcs.size();
    int nFt0s = ft0s.size();
    int nFv0s = fv0s.size();
    int nZdcs = zdcs.size();

    const o2::aod::BCs::iterator& fbc = bcs.begin();
    if (fbc.runNumber() != fRun) {
      getBcInfo(fbc);

      const AxisSpec axisOrbits{fNOrbits / static_cast<int>(fNOrbitsPerTF), 0., static_cast<double>(fNOrbits), ""};
      hOrbitTotal = hr.add<TH1>("hOrbitTotal", "", kTH1D, {axisOrbits}).get();
      hOrbit0n0n = hr.add<TH1>("hOrbit0n0n", "", kTH1D, {axisOrbits}).get();
      hOrbit0nXn = hr.add<TH1>("hOrbit0nXn", "", kTH1D, {axisOrbits}).get();
      hOrbitXn0n = hr.add<TH1>("hOrbitXn0n", "", kTH1D, {axisOrbits}).get();
      hOrbitXnXn = hr.add<TH1>("hOrbitXnXn", "", kTH1D, {axisOrbits}).get();
      hOrbitZNA = hr.add<TH1>("hOrbitZNA", "", kTH1D, {axisOrbits}).get();
      hOrbitZNC = hr.add<TH1>("hOrbitZNC", "", kTH1D, {axisOrbits}).get();
      hOrbitTCE = hr.add<TH1>("hOrbitTCE", "", kTH1D, {axisOrbits}).get();
      hOrbitV0A = hr.add<TH1>("hOrbitV0A", "", kTH1D, {axisOrbits}).get();
      hOrbitT0A = hr.add<TH1>("hOrbitT0A", "", kTH1D, {axisOrbits}).get();

      double minSec = std::floor(fTsSOR / 1000.);
      double maxSec = std::ceil(fTsEOR / 1000.);
      const AxisSpec axisSeconds{static_cast<int>(maxSec - minSec), minSec, maxSec, "seconds"};
      hSecondsTotal = hr.add<TH1>("hSecondsTotal", "", kTH1D, {axisSeconds}).get();
      hSeconds0n0n = hr.add<TH1>("hSeconds0n0n", "", kTH1D, {axisSeconds}).get();
      hSeconds0nXn = hr.add<TH1>("hSeconds0nXn", "", kTH1D, {axisSeconds}).get();
      hSecondsXn0n = hr.add<TH1>("hSecondsXn0n", "", kTH1D, {axisSeconds}).get();
      hSecondsXnXn = hr.add<TH1>("hSecondsXnXn", "", kTH1D, {axisSeconds}).get();
      hSecondsZNA = hr.add<TH1>("hSecondsZNA", "", kTH1D, {axisSeconds}).get();
      hSecondsZNC = hr.add<TH1>("hSecondsZNC", "", kTH1D, {axisSeconds}).get();
      hSecondsTCE = hr.add<TH1>("hSecondsTCE", "", kTH1D, {axisSeconds}).get();
      hSecondsV0A = hr.add<TH1>("hSecondsV0A", "", kTH1D, {axisSeconds}).get();
      hSecondsT0A = hr.add<TH1>("hSecondsT0A", "", kTH1D, {axisSeconds}).get();
    }

    std::unordered_set<uint64_t> orbitNumbers{}; // non-empty unique orbit numbers

    std::vector<int64_t> gbcs(nBcs, -1);
    std::vector<int64_t> gbcsTs(nBcs, -1);
    std::vector<int64_t> gbcsOrbit(nBcs, -1);

    int64_t nBCsPat_all = 0;
    int64_t nBCsPat_total = 0;
    int64_t nBCsPat_0n0n = 0;
    int64_t nBCsPat_0nXn = 0;
    int64_t nBCsPat_Xn0n = 0;
    int64_t nBCsPat_XnXn = 0;
    int64_t nBCsPat_ZNA = 0;
    int64_t nBCsPat_ZNC = 0;
    int64_t nBCsPat_TCE = 0;

    for (auto i = 0; i < nBcs; ++i) {
      const auto& bc = bcs.iteratorAt(i);
      if (bcPatternB[bc.globalBC() % o2::constants::lhc::LHCMaxBunches] == 0)
        continue;
      gbcs[i] = bc.globalBC();
      uint64_t orbit = bc.globalBC() / o2::constants::lhc::LHCMaxBunches;
      uint64_t ts = (bc.globalBC() * o2::constants::lhc::LHCBunchSpacingMUS + fTsOrbitReset) / 1000000; // seconds
      gbcsOrbit[i] = orbit;
      gbcsTs[i] = ts;
      orbitNumbers.insert(orbit);
      hOrbitTotal->Fill(orbit - fMinOrbit);
      hSecondsTotal->Fill(ts);
      nBCsPat_total++;
    }

    nBCsPat_all = orbitNumbers.size() * o2::constants::lhc::LHCMaxBunches;

    std::map<int64_t, float> mapGlobalBcWithZNA{};
    std::map<int64_t, float> mapGlobalBcWithZNC{};
    for (auto i = 0; i < nZdcs; ++i) {
      const auto& zdc = zdcs.iteratorAt(i);
      auto bcId = zdc.bcId();
      int64_t globalBC = gbcs[bcId];
      if (globalBC == -1)
        continue;
      auto orbit = gbcsOrbit[bcId];
      auto ts = gbcsTs[bcId];
      float timeA = zdc.timeZNA();
      float timeC = zdc.timeZNC();
      hTimeZNA->Fill(timeA);
      hTimeZNC->Fill(timeC);
      bool hasZNA = !(std::abs(timeA) > 2.f);
      bool hasZNC = !(std::abs(timeC) > 2.f);
      mapGlobalBcWithZNA[globalBC] = timeA;
      mapGlobalBcWithZNC[globalBC] = timeC;
      if (!hasZNA && !hasZNC) {
        nBCsPat_0n0n++;
        hOrbit0n0n->Fill(orbit - fMinOrbit);
        hSeconds0n0n->Fill(ts);
      }
      if (!hasZNA && hasZNC) {
        nBCsPat_0nXn++;
        hOrbit0nXn->Fill(orbit - fMinOrbit);
        hSeconds0nXn->Fill(ts);
      }
      if (hasZNA && !hasZNC) {
        nBCsPat_Xn0n++;
        hOrbitXn0n->Fill(orbit - fMinOrbit);
        hSecondsXn0n->Fill(ts);
      }
      if (hasZNA && hasZNC) {
        nBCsPat_XnXn++;
        hOrbitXnXn->Fill(orbit - fMinOrbit);
        hSecondsXnXn->Fill(ts);
      }
      if (hasZNA) {
        nBCsPat_ZNA++;
        hOrbitZNA->Fill(orbit - fMinOrbit);
        hSecondsZNA->Fill(ts);
      }
      if (hasZNC) {
        nBCsPat_ZNC++;
        hOrbitZNC->Fill(orbit - fMinOrbit);
        hSecondsZNC->Fill(ts);
      }
    }

    // looking for bcs without any zdc signals --> adding to 0n0n
    for (auto i = 0; i < nBcs; ++i) {
      auto globalBC = gbcs[i];
      if (globalBC == -1)
        continue;
      auto itA = mapGlobalBcWithZNA.find(globalBC);
      auto itC = mapGlobalBcWithZNC.find(globalBC);
      if (itA == mapGlobalBcWithZNA.end() && itC == mapGlobalBcWithZNC.end()) {
        auto orbit = gbcsOrbit[i];
        auto ts = gbcsTs[i];
        hOrbit0n0n->Fill(orbit - fMinOrbit);
        hSeconds0n0n->Fill(ts);
        nBCsPat_0n0n++;
      }
    }

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
      const auto& amps = fv0.amplitude();
      auto amp = std::accumulate(amps.begin(), amps.end(), 0.f);
      mapGlobalBcWithV0A[globalBC] = amp;
    }

    std::map<int64_t, float> mapGlobalBcWithT0A{};
    for (auto i = 0; i < nFt0s; ++i) {
      const auto& ft0 = ft0s.iteratorAt(i);
      if (!TESTBIT(ft0.triggerMask(), o2::fit::Triggers::bitVertex))
        continue;
      if (TESTBIT(ft0.triggerMask(), o2::fit::Triggers::bitCen)) {
        auto orbit = gbcsOrbit[i];
        auto ts = gbcsTs[i];
        hOrbitTCE->Fill(orbit - fMinOrbit);
        hSecondsTCE->Fill(ts);
        nBCsPat_TCE++;
      }
      auto bcId = ft0.bcId();
      int64_t globalBC = gbcs[bcId];
      if (globalBC == -1)
        continue;
      if (std::abs(ft0.timeA()) > 2.f)
        continue;
      const auto& amps = ft0.amplitudeA();
      auto amp = std::accumulate(amps.begin(), amps.end(), 0.f);
      mapGlobalBcWithT0A[globalBC] = amp;
    }

    hCountBCs->Fill("All", static_cast<double>(nBCsPat_all));
    hCountBCs->Fill("Total", static_cast<double>(nBCsPat_total));
    hCountBCs->Fill("0n0n", static_cast<double>(nBCsPat_0n0n));
    hCountBCs->Fill("0nXn", static_cast<double>(nBCsPat_0nXn));
    hCountBCs->Fill("Xn0n", static_cast<double>(nBCsPat_Xn0n));
    hCountBCs->Fill("XnXn", static_cast<double>(nBCsPat_XnXn));
    hCountBCs->Fill("ZNA", static_cast<double>(nBCsPat_ZNA));
    hCountBCs->Fill("ZNC", static_cast<double>(nBCsPat_ZNC));
    hCountBCs->Fill("TCE", static_cast<double>(nBCsPat_TCE));

    hSelBCAmpV0A_Total.object->Fill(0.f, 0.f, static_cast<double>(nBCsPat_total));
    hSelBCAmpV0A_0n0n.object->Fill(0.f, 0.f, static_cast<double>(nBCsPat_0n0n));
    hSelBCAmpV0A_0nXn.object->Fill(0.f, 0.f, static_cast<double>(nBCsPat_0nXn));
    hSelBCAmpV0A_Xn0n.object->Fill(0.f, 0.f, static_cast<double>(nBCsPat_Xn0n));
    hSelBCAmpV0A_XnXn.object->Fill(0.f, 0.f, static_cast<double>(nBCsPat_XnXn));

    hSelBCAmpT0A_Total.object->Fill(0.f, 0.f, static_cast<double>(nBCsPat_total));
    hSelBCAmpT0A_0n0n.object->Fill(0.f, 0.f, static_cast<double>(nBCsPat_0n0n));
    hSelBCAmpT0A_0nXn.object->Fill(0.f, 0.f, static_cast<double>(nBCsPat_0nXn));
    hSelBCAmpT0A_Xn0n.object->Fill(0.f, 0.f, static_cast<double>(nBCsPat_Xn0n));
    hSelBCAmpT0A_XnXn.object->Fill(0.f, 0.f, static_cast<double>(nBCsPat_XnXn));

    collectVeto(gbcs, gbcsOrbit, gbcsTs, "V0A", mapGlobalBcWithV0A, mapGlobalBcWithZNA, mapGlobalBcWithZNC);
    collectVeto(gbcs, gbcsOrbit, gbcsTs, "T0A", mapGlobalBcWithT0A, mapGlobalBcWithZNA, mapGlobalBcWithZNC);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UpcVetoAnalysis>(cfgc)};
}
