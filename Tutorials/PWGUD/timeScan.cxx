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

/// \file timeScan.cxx
/// \brief Example task to scan data vs time TF-by-TF in single AO2D.root
/// Run as: o2-analysis-ud-time-structure --aod-file AO2D.root
/// Filling 1D histos (fill1D=1) might require large shm size: --shm-segment-size 14000000000
/// \author Evgeny Kryshen

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/LHCConstants.h>
#include <DataFormatsParameters/AggregatedRunInfo.h>
#include <DataFormatsParameters/GRPLHCIFData.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>
#include <TString.h>

#include <bitset>
#include <cstdint>
#include <map>
#include <memory>
#include <numeric>

using namespace o2;
using namespace o2::framework;

struct TimeScan {
  // Need several histogram registries since each HistogramRegistry has a limit of 512 histos
  HistogramRegistry histos{"histos", {}};
  HistogramRegistry histosCollidingBCs{"histosCollidingBCs", {}};
  HistogramRegistry histosBCs{"histosBCs", {}};
  HistogramRegistry histosCols{"histosCols", {}};
  HistogramRegistry histosFT0A{"histosFT0A", {}};
  HistogramRegistry histosFT0C{"histosFT0C", {}};
  HistogramRegistry histosFV0A{"histosFV0A", {}};
  HistogramRegistry histosFDDA{"histosFDDA", {}};
  HistogramRegistry histosFDDC{"histosFDDC", {}};
  HistogramRegistry histosZNA{"histosZNA", {}};
  HistogramRegistry histosZNC{"histosZNC", {}};
  Service<o2::ccdb::BasicCCDBManager> ccdb{};

  int32_t cachedRunNumber = -1;
  int64_t bcSOR = 0;
  int32_t nBCsPerTF = 32 * o2::constants::lhc::LHCMaxBunches;
  std::bitset<o2::constants::lhc::LHCMaxBunches> collidingBCs;

  std::map<int, std::shared_ptr<TH1>> mapHistoCollidingBCs;
  std::map<int, std::shared_ptr<TH1>> mapHistoBCs;
  std::map<int, std::shared_ptr<TH1>> mapHistoCols;
  std::map<int, std::shared_ptr<TH1>> mapHistoFT0A;
  std::map<int, std::shared_ptr<TH1>> mapHistoFT0C;
  std::map<int, std::shared_ptr<TH1>> mapHistoFV0A;
  std::map<int, std::shared_ptr<TH1>> mapHistoFDDA;
  std::map<int, std::shared_ptr<TH1>> mapHistoFDDC;
  std::map<int, std::shared_ptr<TH1>> mapHistoZNA;
  std::map<int, std::shared_ptr<TH1>> mapHistoZNC;

  Configurable<int> fill1D{"fill1D", 0, "Fill 1D histos: 0=off, 1=on"};

  void init(InitContext&)
  {
    ccdb->setLocalObjectValidityChecking();
    const AxisSpec axisBCinTF{nBCsPerTF, 0., static_cast<double>(nBCsPerTF), "BC in TF"};
    const AxisSpec axisTF{1, 0., 1., "TF"};
    const AxisSpec axisNumContrib{1000, 0., 8000., "Number of contributors"};
    const AxisSpec axisSumAmpFV0A{2000, 0., 200000., "V0A amplitude"};
    const AxisSpec axisSumAmpFT0A{2000, 0., 200000., "T0A amplitude"};
    const AxisSpec axisSumAmpFT0C{2000, 0., 70000., "T0C amplitude"};
    const AxisSpec axisSumAmpFDDA{2000, 0., 40000., "FDA amplitude"};
    const AxisSpec axisSumAmpFDDC{2000, 0., 40000., "FDC amplitude"};
    const AxisSpec axisSumAmpZNA{3000, 0., 300., "ZNA amplitude"};
    const AxisSpec axisSumAmpZNC{3000, 0., 300., "ZNC amplitude"};

    histos.add("hCollidingBCs", "Colliding BCs", HistType::kTH2D, {axisTF, axisBCinTF});
    histos.add("hBCs", "Non-empty BCs", HistType::kTH2D, {axisTF, axisBCinTF});
    histos.add("hCols", "Number of contributors", HistType::kTH2D, {axisTF, axisBCinTF});
    histos.add("hFT0A", "FT0A amplitude", HistType::kTH2D, {axisTF, axisBCinTF});
    histos.add("hFT0C", "FT0C amplitude", HistType::kTH2D, {axisTF, axisBCinTF});
    histos.add("hFV0A", "FV0A amplitude", HistType::kTH2D, {axisTF, axisBCinTF});
    histos.add("hFDDA", "FDDA amplitude", HistType::kTH2D, {axisTF, axisBCinTF});
    histos.add("hFDDC", "FDDC amplitude", HistType::kTH2D, {axisTF, axisBCinTF});
    histos.add("hZNA", "ZNA amplitude", HistType::kTH2D, {axisTF, axisBCinTF});
    histos.add("hZNC", "ZNC amplitude", HistType::kTH2D, {axisTF, axisBCinTF});

    histos.add("hNumContrib", "Number of contributors", HistType::kTH1D, {axisNumContrib});
    histos.add("hSumAmpFV0A", "FV0A amplitude", HistType::kTH1D, {axisSumAmpFV0A});
    histos.add("hSumAmpFT0A", "FT0A amplitude", HistType::kTH1D, {axisSumAmpFT0A});
    histos.add("hSumAmpFT0C", "FT0C amplitude", HistType::kTH1D, {axisSumAmpFT0C});
    histos.add("hSumAmpFDDA", "FDDA amplitude", HistType::kTH1D, {axisSumAmpFDDA});
    histos.add("hSumAmpFDDC", "FDDC amplitude", HistType::kTH1D, {axisSumAmpFDDC});
    histos.add("hSumAmpZNA", "ZNA amplitude", HistType::kTH1D, {axisSumAmpZNA});
    histos.add("hSumAmpZNC", "ZNC amplitude", HistType::kTH1D, {axisSumAmpZNC});
  }

  void process(aod::BCs const& bcs, aod::Collisions const& collisions, aod::FT0s const& ft0s, aod::FV0As const& fv0s, aod::FDDs const& fdds, aod::Zdcs const& zdcs)
  {
    LOGP(info, "Processing new DF");

    int32_t runNumber = bcs.begin().runNumber();
    if (runNumber != cachedRunNumber) {
      LOGP(debug, "Updating run info (called once per run)");
      cachedRunNumber = runNumber;
      const auto runInfo = o2::parameters::AggregatedRunInfo::buildAggregatedRunInfo(ccdb->instance(), runNumber);
      bcSOR = runInfo.orbitSOR * o2::constants::lhc::LHCMaxBunches;
      nBCsPerTF = runInfo.orbitsPerTF * o2::constants::lhc::LHCMaxBunches;
      collidingBCs = runInfo.grpLHC->getBunchFilling().getBCPattern();
    }

    for (const auto& bc : bcs) {
      const int64_t globalBC = bc.globalBC();
      const int32_t tfID = (globalBC - bcSOR) / nBCsPerTF;
      if (!mapHistoCollidingBCs.contains(tfID)) {
        LOGP(debug, "Fill maps of 1D histograms for new TF");
        const AxisSpec axisBCinTF{nBCsPerTF, 0., static_cast<double>(nBCsPerTF), "BC in TF"};
        mapHistoCollidingBCs[tfID] = histosCollidingBCs.add<TH1>(Form("%d/hCollidingBCs", tfID), "Colliding BCs", HistType::kTH1D, {axisBCinTF});
        if (fill1D > 0) {
          mapHistoBCs[tfID] = histosBCs.add<TH1>(Form("%d/hBCs", tfID), "Non-empty BCs", HistType::kTH1D, {axisBCinTF});
          mapHistoCols[tfID] = histosCols.add<TH1>(Form("%d/hCols", tfID), "Number of contributors", HistType::kTH1D, {axisBCinTF});
          mapHistoFT0A[tfID] = histosFT0A.add<TH1>(Form("%d/hFT0A", tfID), "FT0A amplitude", HistType::kTH1D, {axisBCinTF});
          mapHistoFT0C[tfID] = histosFT0C.add<TH1>(Form("%d/hFT0C", tfID), "FT0C amplitude", HistType::kTH1D, {axisBCinTF});
          mapHistoFV0A[tfID] = histosFV0A.add<TH1>(Form("%d/hFV0A", tfID), "FV0A amplitude", HistType::kTH1D, {axisBCinTF});
          mapHistoFDDA[tfID] = histosFDDA.add<TH1>(Form("%d/hFDDA", tfID), "FDDA amplitude", HistType::kTH1D, {axisBCinTF});
          mapHistoFDDC[tfID] = histosFDDC.add<TH1>(Form("%d/hFDDC", tfID), "FDDC amplitude", HistType::kTH1D, {axisBCinTF});
          mapHistoZNA[tfID] = histosZNA.add<TH1>(Form("%d/hZNA", tfID), "ZNA amplitude", HistType::kTH1D, {axisBCinTF});
          mapHistoZNC[tfID] = histosZNC.add<TH1>(Form("%d/hZNC", tfID), "ZNC amplitude", HistType::kTH1D, {axisBCinTF});
        }

        LOGP(debug, "Filling colliding bcs");
        for (int32_t bcInTF = 0; bcInTF < nBCsPerTF; bcInTF++) {
          if (collidingBCs[bcInTF % o2::constants::lhc::LHCMaxBunches]) {
            mapHistoCollidingBCs[tfID]->Fill(bcInTF);
            histos.get<TH2>(HIST("hCollidingBCs"))->Fill(Form("%d", tfID), bcInTF, 1.);
          }
        }
      }
    }

    LOGP(debug, "Filling bc activity");
    for (const auto& bc : bcs) {
      const int64_t globalBC = bc.globalBC();
      const int32_t tfID = (globalBC - bcSOR) / nBCsPerTF;
      const int32_t bcInTF = (globalBC - bcSOR) % nBCsPerTF;
      LOGP(debug, "globalBC={} tfID={} bcInTF={}", globalBC, tfID, bcInTF);
      histos.get<TH2>(HIST("hBCs"))->Fill(Form("%d", tfID), bcInTF, 1.);
      if (fill1D > 0) {
        mapHistoBCs[tfID]->Fill(bcInTF);
      }
    }

    LOGP(debug, "Filling collision activity");
    for (const auto& collision : collisions) {
      const int64_t globalBC = collision.bc().globalBC();
      const int32_t tfID = (globalBC - bcSOR) / nBCsPerTF;
      const int32_t bcInTF = (globalBC - bcSOR) % nBCsPerTF;
      const int32_t nContributors = collision.numContrib();
      LOGP(debug, "globalBC={} nContributors={}", globalBC, nContributors);
      histos.get<TH2>(HIST("hCols"))->Fill(Form("%d", tfID), bcInTF, nContributors);
      histos.get<TH1>(HIST("hNumContrib"))->Fill(nContributors);
      if (fill1D > 0) {
        mapHistoCols[tfID]->Fill(bcInTF, nContributors);
      }
    }

    LOGP(debug, "Filling FT0 activity");
    for (const auto& ft0 : ft0s) {
      const int64_t globalBC = ft0.bc().globalBC();
      const int32_t tfID = (globalBC - bcSOR) / nBCsPerTF;
      const int32_t bcInTF = (globalBC - bcSOR) % nBCsPerTF;
      const float sumAmpFT0A = ft0.sumAmpA();
      const float sumAmpFT0C = ft0.sumAmpC();
      LOGP(debug, "globalBC={} sumAmpFT0A={} sumAmpFT0C", globalBC, sumAmpFT0A, sumAmpFT0C);
      histos.get<TH2>(HIST("hFT0A"))->Fill(Form("%d", tfID), bcInTF, sumAmpFT0A);
      histos.get<TH2>(HIST("hFT0C"))->Fill(Form("%d", tfID), bcInTF, sumAmpFT0C);
      histos.get<TH1>(HIST("hSumAmpFT0A"))->Fill(sumAmpFT0A);
      histos.get<TH1>(HIST("hSumAmpFT0C"))->Fill(sumAmpFT0C);
      if (fill1D > 0) {
        mapHistoFT0A[tfID]->Fill(bcInTF, sumAmpFT0A);
        mapHistoFT0C[tfID]->Fill(bcInTF, sumAmpFT0C);
      }
    }

    LOGP(debug, "Filling FV0 activity");
    for (const auto& fv0 : fv0s) {
      const int64_t globalBC = fv0.bc().globalBC();
      const int32_t tfID = (globalBC - bcSOR) / nBCsPerTF;
      const int32_t bcInTF = (globalBC - bcSOR) % nBCsPerTF;
      const float sumAmpFV0A = std::accumulate(fv0.amplitude().begin(), fv0.amplitude().end(), 0.f);
      LOGP(debug, "globalBC={} sumAmpFV0A={}", globalBC, sumAmpFV0A);
      histos.get<TH2>(HIST("hFV0A"))->Fill(Form("%d", tfID), bcInTF, sumAmpFV0A);
      histos.get<TH1>(HIST("hSumAmpFV0A"))->Fill(sumAmpFV0A);
      if (fill1D > 0) {
        mapHistoFV0A[tfID]->Fill(bcInTF, sumAmpFV0A);
      }
    }

    LOGP(debug, "Filling FDD activity");
    for (const auto& fdd : fdds) {
      const int64_t globalBC = fdd.bc().globalBC();
      const int32_t tfID = (globalBC - bcSOR) / nBCsPerTF;
      const int32_t bcInTF = (globalBC - bcSOR) % nBCsPerTF;
      const float sumAmpFDDA = std::accumulate(std::begin(fdd.chargeA()), std::end(fdd.chargeA()), 0.f);
      const float sumAmpFDDC = std::accumulate(std::begin(fdd.chargeC()), std::end(fdd.chargeC()), 0.f);
      LOGP(debug, "globalBC={} sumAmpFDDA={} sumAmpFDDC={}", globalBC, sumAmpFDDA, sumAmpFDDC);
      histos.get<TH2>(HIST("hFDDA"))->Fill(Form("%d", tfID), bcInTF, sumAmpFDDA);
      histos.get<TH2>(HIST("hFDDC"))->Fill(Form("%d", tfID), bcInTF, sumAmpFDDC);
      histos.get<TH1>(HIST("hSumAmpFDDA"))->Fill(sumAmpFDDA);
      histos.get<TH1>(HIST("hSumAmpFDDC"))->Fill(sumAmpFDDC);
      if (fill1D > 0) {
        mapHistoFDDA[tfID]->Fill(bcInTF, sumAmpFDDA);
        mapHistoFDDC[tfID]->Fill(bcInTF, sumAmpFDDC);
      }
    }

    LOGP(debug, "Filling ZDC activity");
    for (const auto& zdc : zdcs) {
      const int64_t globalBC = zdc.bc().globalBC();
      const int32_t tfID = (globalBC - bcSOR) / nBCsPerTF;
      const int32_t bcInTF = (globalBC - bcSOR) % nBCsPerTF;
      const float sumAmpZNA = zdc.energyCommonZNA() > 0 ? zdc.energyCommonZNA() : 0;
      const float sumAmpZNC = zdc.energyCommonZNC() > 0 ? zdc.energyCommonZNC() : 0;
      LOGP(debug, "globalBC={} sumAmpZNA={} sumAmpZNC={}", globalBC, sumAmpZNA, sumAmpZNC);
      histos.get<TH2>(HIST("hZNA"))->Fill(Form("%d", tfID), bcInTF, sumAmpZNA);
      histos.get<TH2>(HIST("hZNC"))->Fill(Form("%d", tfID), bcInTF, sumAmpZNC);
      histos.get<TH1>(HIST("hSumAmpZNA"))->Fill(sumAmpZNA);
      histos.get<TH1>(HIST("hSumAmpZNC"))->Fill(sumAmpZNC);
      if (fill1D > 0) {
        mapHistoZNA[tfID]->Fill(bcInTF, sumAmpZNA);
        mapHistoZNC[tfID]->Fill(bcInTF, sumAmpZNC);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& context) { return WorkflowSpec{adaptAnalysisTask<TimeScan>(context)}; }
