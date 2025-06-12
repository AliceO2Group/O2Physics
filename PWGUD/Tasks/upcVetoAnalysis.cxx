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
#include <map>
#include <string>
#include <unordered_set>
#include <vector>
#include <unordered_map>
#include <memory>
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
  bool fIsCreated = false;

  // amplitude axis
  int fNbAmp = 101;
  double fMinAmp = 0.;
  double fMaxAmp = 1000;

  // bcs axis
  int fNbBcs = 52;
  double fMinBcs = -1.;
  double fMaxBcs = 50.;

  // zdc time axis
  int fNbT = 101;
  double fMinT = -10.;
  double fMaxT = 10.;

  HistogramRegistry hr{"hr", {}, OutputObjHandlingPolicy::AnalysisObject};

  OutputObj<TH1D> hCountBCs{TH1D("hCountBCs", "", 10, 0., 10.)};

  OutputObj<TH1D> hTimeZNA{TH1D("hTimeZNA", "", fNbT, fMinT, fMaxT)};
  OutputObj<TH1D> hTimeZNC{TH1D("hTimeZNC", "", fNbT, fMinT, fMaxT)};

  OutputObj<TH2D> hSelBCAmpV0A_Total{TH2D("hSelBCAmpV0A_Total", "", fNbAmp, fMinAmp, fMaxAmp, fNbBcs, fMinBcs, fMaxBcs)};
  OutputObj<TH2D> hSelBCAmpV0A_0n0n{TH2D("hSelBCAmpV0A_0n0n", "", fNbAmp, fMinAmp, fMaxAmp, fNbBcs, fMinBcs, fMaxBcs)};
  OutputObj<TH2D> hSelBCAmpV0A_0nXn{TH2D("hSelBCAmpV0A_0nXn", "", fNbAmp, fMinAmp, fMaxAmp, fNbBcs, fMinBcs, fMaxBcs)};
  OutputObj<TH2D> hSelBCAmpV0A_Xn0n{TH2D("hSelBCAmpV0A_Xn0n", "", fNbAmp, fMinAmp, fMaxAmp, fNbBcs, fMinBcs, fMaxBcs)};
  OutputObj<TH2D> hSelBCAmpV0A_XnXn{TH2D("hSelBCAmpV0A_XnXn", "", fNbAmp, fMinAmp, fMaxAmp, fNbBcs, fMinBcs, fMaxBcs)};

  OutputObj<TH2D> hSelBCAmpT0A_Total{TH2D("hSelBCAmpT0A_Total", "", fNbAmp, fMinAmp, fMaxAmp, fNbBcs, fMinBcs, fMaxBcs)};
  OutputObj<TH2D> hSelBCAmpT0A_0n0n{TH2D("hSelBCAmpT0A_0n0n", "", fNbAmp, fMinAmp, fMaxAmp, fNbBcs, fMinBcs, fMaxBcs)};
  OutputObj<TH2D> hSelBCAmpT0A_0nXn{TH2D("hSelBCAmpT0A_0nXn", "", fNbAmp, fMinAmp, fMaxAmp, fNbBcs, fMinBcs, fMaxBcs)};
  OutputObj<TH2D> hSelBCAmpT0A_Xn0n{TH2D("hSelBCAmpT0A_Xn0n", "", fNbAmp, fMinAmp, fMaxAmp, fNbBcs, fMinBcs, fMaxBcs)};
  OutputObj<TH2D> hSelBCAmpT0A_XnXn{TH2D("hSelBCAmpT0A_XnXn", "", fNbAmp, fMinAmp, fMaxAmp, fNbBcs, fMinBcs, fMaxBcs)};

  TH1* hOrbitTotal;
  TH1* hOrbit0n0n;
  TH1* hOrbit0nXn;
  TH1* hOrbitXn0n;
  TH1* hOrbitXnXn;
  TH1* hOrbitZNA;
  TH1* hOrbitZNC;
  TH1* hOrbitTCE;

  TH1* hOrbitV0A_Total;
  TH1* hOrbitV0A_0n0n;
  TH1* hOrbitV0A_0nXn;
  TH1* hOrbitV0A_Xn0n;
  TH1* hOrbitV0A_XnXn;
  TH1* hOrbitT0A_Total;
  TH1* hOrbitT0A_0n0n;
  TH1* hOrbitT0A_0nXn;
  TH1* hOrbitT0A_Xn0n;
  TH1* hOrbitT0A_XnXn;

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
  std::vector<int64_t> bbcs{}; // b-mask bcs inside orbit

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
    auto soreor = o2::ccdb::BasicCCDBManager::getRunDuration(ccdb_api, fRun);
    fTsSOR = soreor.first;
    fTsEOR = soreor.second;
    LOGP(info, "fTsSOR={}, fTsEOR={}", fTsSOR, fTsEOR);

    std::map<std::string, std::string> metadata;
    auto grplhcif = ccdb_api.retrieveFromTFileAny<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", metadata, fTsSOR);
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
    if (fbc.runNumber() != fRun && !fIsCreated) {
      fIsCreated = true;
      getBcInfo(fbc);

      for (auto i = 0; i < o2::constants::lhc::LHCMaxBunches; ++i) {
        if (bcPatternB[i] == 0)
          continue;
        bbcs.push_back(i);
      }

      const AxisSpec axisOrbits{fNOrbits / static_cast<int>(fNOrbitsPerTF), 0., static_cast<double>(fNOrbits), ""};
      hOrbitTotal = hr.add<TH1>("hOrbitTotal", "", kTH1D, {axisOrbits}).get();
      hOrbit0n0n = hr.add<TH1>("hOrbit0n0n", "", kTH1D, {axisOrbits}).get();
      hOrbit0nXn = hr.add<TH1>("hOrbit0nXn", "", kTH1D, {axisOrbits}).get();
      hOrbitXn0n = hr.add<TH1>("hOrbitXn0n", "", kTH1D, {axisOrbits}).get();
      hOrbitXnXn = hr.add<TH1>("hOrbitXnXn", "", kTH1D, {axisOrbits}).get();
      hOrbitZNA = hr.add<TH1>("hOrbitZNA", "", kTH1D, {axisOrbits}).get();
      hOrbitZNC = hr.add<TH1>("hOrbitZNC", "", kTH1D, {axisOrbits}).get();
      hOrbitTCE = hr.add<TH1>("hOrbitTCE", "", kTH1D, {axisOrbits}).get();
      hOrbitV0A_Total = hr.add<TH1>("hOrbitV0A_Total", "", kTH1D, {axisOrbits}).get();
      hOrbitV0A_0n0n = hr.add<TH1>("hOrbitV0A_0n0n", "", kTH1D, {axisOrbits}).get();
      hOrbitV0A_0nXn = hr.add<TH1>("hOrbitV0A_0nXn", "", kTH1D, {axisOrbits}).get();
      hOrbitV0A_Xn0n = hr.add<TH1>("hOrbitV0A_Xn0n", "", kTH1D, {axisOrbits}).get();
      hOrbitV0A_XnXn = hr.add<TH1>("hOrbitV0A_XnXn", "", kTH1D, {axisOrbits}).get();
      hOrbitT0A_Total = hr.add<TH1>("hOrbitT0A_Total", "", kTH1D, {axisOrbits}).get();
      hOrbitT0A_0n0n = hr.add<TH1>("hOrbitT0A_0n0n", "", kTH1D, {axisOrbits}).get();
      hOrbitT0A_0nXn = hr.add<TH1>("hOrbitT0A_0nXn", "", kTH1D, {axisOrbits}).get();
      hOrbitT0A_Xn0n = hr.add<TH1>("hOrbitT0A_Xn0n", "", kTH1D, {axisOrbits}).get();
      hOrbitT0A_XnXn = hr.add<TH1>("hOrbitT0A_XnXn", "", kTH1D, {axisOrbits}).get();
    }

    std::unordered_set<uint64_t> orbitNumbers{}; // non-empty unique orbit numbers

    std::vector<int64_t> gbcs(nBcs, -1);
    std::vector<int64_t> gbcsOrbit(nBcs, -1);

    std::vector<float> vBcIdsWithV0A(nBcs, -999.f);
    std::vector<float> vBcIdsWithT0A(nBcs, -999.f);

    std::vector<bool> vBcIdsWithZNA(nBcs, false);
    std::vector<bool> vBcIdsWithZNC(nBcs, false);

    std::unordered_map<int64_t, int64_t> bgbcZdc{}; // fired b-mask bcIds with zdc

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
      gbcsOrbit[i] = orbit;
      orbitNumbers.insert(orbit);
      hOrbitTotal->Fill(orbit - fMinOrbit);
      nBCsPat_total++;
    }

    nBCsPat_all = orbitNumbers.size() * fNBcsB;

    for (auto i = 0; i < nZdcs; ++i) {
      const auto& zdc = zdcs.iteratorAt(i);
      auto bcId = zdc.bcId();
      int64_t globalBC = gbcs[bcId];
      if (globalBC == -1)
        continue;
      auto orbit = gbcsOrbit[bcId] - fMinOrbit;
      float timeA = zdc.timeZNA();
      float timeC = zdc.timeZNC();
      hTimeZNA->Fill(timeA);
      hTimeZNC->Fill(timeC);
      bool hasZNA = !(std::abs(timeA) > 2.f);
      bool hasZNC = !(std::abs(timeC) > 2.f);
      vBcIdsWithZNA[bcId] = hasZNA;
      vBcIdsWithZNC[bcId] = hasZNC;
      bgbcZdc[globalBC] = bcId;
      if (!hasZNA && !hasZNC) {
        nBCsPat_0n0n++;
        hOrbit0n0n->Fill(orbit);
      }
      if (!hasZNA && hasZNC) {
        nBCsPat_0nXn++;
        hOrbit0nXn->Fill(orbit);
      }
      if (hasZNA && !hasZNC) {
        nBCsPat_Xn0n++;
        hOrbitXn0n->Fill(orbit);
      }
      if (hasZNA && hasZNC) {
        nBCsPat_XnXn++;
        hOrbitXnXn->Fill(orbit);
      }
      if (hasZNA) {
        nBCsPat_ZNA++;
        hOrbitZNA->Fill(orbit);
      }
      if (hasZNC) {
        nBCsPat_ZNC++;
        hOrbitZNC->Fill(orbit);
      }
    }

    // looking for bcs without any zdc signals --> adding to 0n0n
    for (int32_t i = 0; i < nBcs; ++i) {
      auto globalBC = gbcs[i];
      if (globalBC == -1)
        continue;
      if (!vBcIdsWithZNA[i] && !vBcIdsWithZNC[i]) {
        auto orbit = gbcsOrbit[i] - fMinOrbit;
        hOrbit0n0n->Fill(orbit);
        nBCsPat_0n0n++;
      }
    }

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
      vBcIdsWithV0A[bcId] = amp;
    }

    for (auto i = 0; i < nFt0s; ++i) {
      const auto& ft0 = ft0s.iteratorAt(i);
      if (!TESTBIT(ft0.triggerMask(), o2::fit::Triggers::bitVertex))
        continue;
      if (TESTBIT(ft0.triggerMask(), o2::fit::Triggers::bitCen)) {
        auto orbit = gbcsOrbit[i] - fMinOrbit;
        hOrbitTCE->Fill(orbit);
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
      vBcIdsWithT0A[bcId] = amp;
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

    hSelBCAmpV0A_Total.object->Fill(0.f, -1.f, static_cast<double>(nBCsPat_all));
    hSelBCAmpV0A_0n0n.object->Fill(0.f, -1.f, static_cast<double>(nBCsPat_0n0n));
    hSelBCAmpV0A_0nXn.object->Fill(0.f, -1.f, static_cast<double>(nBCsPat_0nXn));
    hSelBCAmpV0A_Xn0n.object->Fill(0.f, -1.f, static_cast<double>(nBCsPat_Xn0n));
    hSelBCAmpV0A_XnXn.object->Fill(0.f, -1.f, static_cast<double>(nBCsPat_XnXn));

    hSelBCAmpT0A_Total.object->Fill(0.f, -1.f, static_cast<double>(nBCsPat_all));
    hSelBCAmpT0A_0n0n.object->Fill(0.f, -1.f, static_cast<double>(nBCsPat_0n0n));
    hSelBCAmpT0A_0nXn.object->Fill(0.f, -1.f, static_cast<double>(nBCsPat_0nXn));
    hSelBCAmpT0A_Xn0n.object->Fill(0.f, -1.f, static_cast<double>(nBCsPat_Xn0n));
    hSelBCAmpT0A_XnXn.object->Fill(0.f, -1.f, static_cast<double>(nBCsPat_XnXn));

    std::shared_ptr<TH2> hSelBCAmpV0A = nullptr;
    TH1* hOrbitV0A = nullptr;
    std::shared_ptr<TH2> hSelBCAmpT0A = nullptr;
    TH1* hOrbitT0A = nullptr;
    for (auto orb : orbitNumbers) {
      for (auto bbc : bbcs) {
        auto orbit = orb - fMinOrbit;
        auto gbc = orb * o2::constants::lhc::LHCMaxBunches + bbc;
        auto itBcId = bgbcZdc.find(gbc);
        bool hasZNA = false;
        bool hasZNC = false;
        if (itBcId != bgbcZdc.end()) {
          auto zdc_idx = itBcId->second;
          hasZNA = vBcIdsWithZNA[zdc_idx];
          hasZNC = vBcIdsWithZNC[zdc_idx];
        }
        if (!hasZNA && !hasZNC) {
          hSelBCAmpV0A = hSelBCAmpV0A_0n0n.object;
          hOrbitV0A = hOrbitV0A_0n0n;
          hSelBCAmpT0A = hSelBCAmpT0A_0n0n.object;
          hOrbitT0A = hOrbitT0A_0n0n;
        }
        if (!hasZNA && hasZNC) {
          hSelBCAmpV0A = hSelBCAmpV0A_0nXn.object;
          hOrbitV0A = hOrbitV0A_0nXn;
          hSelBCAmpT0A = hSelBCAmpT0A_0nXn.object;
          hOrbitT0A = hOrbitT0A_0nXn;
        }
        if (hasZNA && !hasZNC) {
          hSelBCAmpV0A = hSelBCAmpV0A_Xn0n.object;
          hOrbitV0A = hOrbitV0A_Xn0n;
          hSelBCAmpT0A = hSelBCAmpT0A_Xn0n.object;
          hOrbitT0A = hOrbitT0A_Xn0n;
        }
        if (hasZNA && hasZNC) {
          hSelBCAmpV0A = hSelBCAmpV0A_XnXn.object;
          hOrbitV0A = hOrbitV0A_XnXn;
          hSelBCAmpT0A = hSelBCAmpT0A_XnXn.object;
          hOrbitT0A = hOrbitT0A_XnXn;
        }
        for (auto r = 0; r <= 10; ++r) {
          auto maxAmpV0A = -999.f;
          auto maxAmpT0A = -999.f;
          int64_t s = gbc - r;
          int64_t e = gbc + r;
          auto lower = std::lower_bound(gbcs.begin(), gbcs.end(), s);
          if (lower != gbcs.end()) {
            auto idx = std::distance(gbcs.begin(), lower);
            while (gbcs[idx] >= s && gbcs[idx] <= e && idx < std::ssize(gbcs)) {
              auto aV0A = vBcIdsWithV0A[idx];
              auto aT0A = vBcIdsWithT0A[idx];
              if (aV0A > maxAmpV0A)
                maxAmpV0A = aV0A;
              if (aT0A > maxAmpT0A)
                maxAmpT0A = aT0A;
              ++idx;
            }
          }
          int biny = r + 2;
          for (auto iamp = 1; iamp < 101; ++iamp) {
            float amp = 10.f * iamp;
            int binx = iamp + 1;
            int bin = binx + (101 + 2) * biny;
            if (maxAmpV0A < amp) {
              hSelBCAmpV0A_Total->AddBinContent(bin);
              hSelBCAmpV0A->AddBinContent(bin);
            }
            if (maxAmpT0A < amp) {
              hSelBCAmpT0A_Total->AddBinContent(bin);
              hSelBCAmpT0A->AddBinContent(bin);
            }
          }
          if (r == 0 && maxAmpV0A < 100.f) {
            hOrbitV0A_Total->Fill(orbit);
            hOrbitV0A->Fill(orbit);
          }
          if (r == 0 && maxAmpT0A < 100.f) {
            hOrbitT0A_Total->Fill(orbit);
            hOrbitT0A->Fill(orbit);
          }
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
