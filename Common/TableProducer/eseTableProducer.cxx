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

/// \file eseTableProducer.cxx
/// \brief Producer for the ESE table
///
/// \author Joachim C. K. B. Hansen

#include <CCDB/BasicCCDBManager.h>

#include <chrono>
#include <string>
#include <algorithm>
#include <numeric>
#include <vector>

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/runDataProcessing.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"

#include "Common/DataModel/EseTable.h"
#include "Common/DataModel/Qvectors.h"
#include "FFitWeights.h"

using namespace o2;
using namespace o2::framework;

using CollWithMults = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::QvectorFT0CVecs, aod::QvectorFT0AVecs, aod::QvectorFV0AVecs>;

struct EseTableProducer {
  Produces<o2::aod::QPercentileFT0Cs> qPercsFT0C;
  Produces<o2::aod::QPercentileFT0As> qPercsFT0A;
  Produces<o2::aod::QPercentileFV0As> qPercsFV0A;
  Produces<o2::aod::QPercentileTPCs> qPercsTPC;
  Produces<o2::aod::FEseCols> fEseCol;

  OutputObj<FFitWeights> FFitObj{FFitWeights("weights")};
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  Configurable<bool> cfgESE{"cfgESE", 1, "ese activation step: false = no ese, true = evaluate splines and fill table"};
  Configurable<std::string> cfgEsePath{"cfgEsePath", "Users/j/joachiha/ESE/local/ffitsplines", "CCDB path for ese splines"};
  Configurable<bool> cfgFT0C{"cfgFT0C", 1, "FT0C flag"};
  Configurable<bool> cfgFT0A{"cfgFT0A", 0, "FT0A flag"};
  Configurable<bool> cfgFV0A{"cfgFV0A", 0, "FV0A flag"};
  Configurable<bool> cfgTPC{"cfgTPC", 0, "TPC flag"};
  Configurable<std::vector<int>> cfgLoopHarmonics{"cfgLoopHarmonics", {2, 3}, "Harmonics to loop over when filling and evaluating splines"};
  Configurable<std::string> cfgCentEst{"cfgCentEst", "FT0C", "centrality estimator"};

  Configurable<std::vector<float>> cfgaxisqn{"cfgaxisqn", {500, 0, 25}, "q_n amplitude range"};
  Configurable<int> cfgnResolution{"cfgnResolution", 3000, "resolution of splines"};

  int runNumber{-1};

  FFitWeights* splines{nullptr};

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  void init(o2::framework::InitContext&)
  {

    LOGF(info, "ESETable::init()");

    registry.add("hEventCounter", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
    registry.add("hESEstat", "ese status;ese status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    std::vector<std::string> vecStr{};
    if (cfgFT0C)
      vecStr.push_back("FT0C");
    if (cfgFT0A)
      vecStr.push_back("FT0A");
    if (cfgFV0A)
      vecStr.push_back("FV0A");
    if (cfgTPC)
      vecStr.push_back("TPC");

    std::vector<std::pair<int, std::string>> veccfg;
    for (std::size_t i{0}; i < cfgLoopHarmonics->size(); i++) {
      for (const auto& j : vecStr) {
        veccfg.push_back({cfgLoopHarmonics->at(i), j});
      }
    }
    FFitObj->SetBinAxis(cfgaxisqn->at(0), cfgaxisqn->at(1), cfgaxisqn->at(2));
    FFitObj->SetResolution(cfgnResolution);
    FFitObj->SetQnType(veccfg);
    FFitObj->Init();
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {

    auto timestamp = bc.timestamp();

    if (cfgESE) {
      splines = ccdb->getForTimeStamp<FFitWeights>(cfgEsePath, timestamp);
      if (!splines)
        LOGF(fatal, "failed loading splines with ese flag");
      LOGF(info, "successfully loaded splines");
    }
  }

  float Calcqn(const float& Qx, const float& Qy, const float& Mult)
  {
    float dqn{0.0f};
    float qn{0.0f};

    dqn = Qx * Qx + Qy * Qy;
    qn = TMath::Sqrt(dqn) / TMath::Sqrt(Mult);
    return qn;
  }

  bool validQvec(const float& qVec)
  {
    if (qVec == 999. || qVec == -999) {
      return false;
    } else {
      return true;
    }
  };

  void doSpline(float& splineVal, int& eseval, const float& centr, const float& nHarm, const char* pf, const auto& QX, const auto& QY, const auto& Ampl)
  {
    if (validQvec(QX[nHarm - 2]) && validQvec(QY[nHarm - 2]) && Ampl > 1e-8) {
      float qnval = Calcqn(QX[nHarm - 2] * Ampl, QY[nHarm - 2] * Ampl, Ampl);
      FFitObj->Fill(centr, qnval, nHarm, pf);
      if (cfgESE) {
        splineVal = splines->EvalSplines(centr, qnval, nHarm, pf);
        eseval = cfgFT0C ? 1 : 0;
      }
    }
  }

  template <typename T>
  void calculateESE(T const& collision,
                    std::vector<float>& qnpFT0C,
                    std::vector<float>& qnpFT0A,
                    std::vector<float>& qnpFV0A,
                    std::vector<float>& qnpTPC,
                    std::vector<int>& fIsEseAvailable)
  {

    const float centrality = collision.centFT0C();
    registry.fill(HIST("hESEstat"), 0.5);
    for (std::size_t i{0}; i < cfgLoopHarmonics->size(); i++) {
      float splineValFT0C{-1.0};
      float splineValFT0A{-1.0};
      float splineValFV0A{-1.0};
      qnpTPC.push_back(-1.0); /* not implemented yet */
      int eseAvailable{0};

      int nHarm = cfgLoopHarmonics->at(i);
      if (cfgFT0C) {
        const auto QxFT0C_Qvec = collision.qvecFT0CReVec();
        const auto QyFT0C_Qvec = collision.qvecFT0CImVec();
        const auto SumAmplFT0C = collision.sumAmplFT0C();
        doSpline(splineValFT0C, eseAvailable, centrality, nHarm, "FT0C", QxFT0C_Qvec, QyFT0C_Qvec, SumAmplFT0C);
        if (i == 0)
          registry.fill(HIST("hESEstat"), 1.5);
      }
      qnpFT0C.push_back(splineValFT0C);
      fIsEseAvailable.push_back(eseAvailable);

      if (cfgFT0A) {
        const auto QxFT0A_Qvec = collision.qvecFT0AReVec();
        const auto QyFT0A_Qvec = collision.qvecFT0AImVec();
        const auto SumAmplFT0A = collision.sumAmplFT0A();
        doSpline(splineValFT0A, eseAvailable, centrality, nHarm, "FT0A", QxFT0A_Qvec, QyFT0A_Qvec, SumAmplFT0A);
        if (i == 0)
          registry.fill(HIST("hESEstat"), 2.5);
      }
      qnpFT0A.push_back(splineValFT0A);

      if (cfgFV0A) {
        const auto QxFV0A_Qvec = collision.qvecFV0AReVec();
        const auto QyFV0A_Qvec = collision.qvecFV0AImVec();
        const auto SumAmplFV0A = collision.sumAmplFV0A();
        doSpline(splineValFV0A, eseAvailable, centrality, nHarm, "FV0A", QxFV0A_Qvec, QyFV0A_Qvec, SumAmplFV0A);
        if (i == 0)
          registry.fill(HIST("hESEstat"), 3.5);
      }
      qnpFV0A.push_back(splineValFV0A);
    }
  };

  void processESE(CollWithMults::iterator const& collision, aod::BCsWithTimestamps const&, aod::FV0As const&, aod::FT0s const&)
  {
    float counter{0.5};
    registry.fill(HIST("hEventCounter"), counter++);

    std::vector<float> qnpFT0C{};
    std::vector<float> qnpFT0A{};
    std::vector<float> qnpFV0A{};
    std::vector<float> qnpTPC{};

    std::vector<int> fIsEseAvailable{};

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int currentRun = bc.runNumber();
    if (runNumber != currentRun) {
      runNumber = currentRun;
      initCCDB(bc);
    }
    registry.fill(HIST("hEventCounter"), counter++);
    calculateESE(collision, qnpFT0C, qnpFT0A, qnpFV0A, qnpTPC, fIsEseAvailable);

    qPercsFT0C(qnpFT0C);
    qPercsFT0A(qnpFT0A);
    qPercsFV0A(qnpFV0A);
    qPercsTPC(qnpTPC);
    fEseCol(fIsEseAvailable);
    registry.fill(HIST("hEventCounter"), counter++);
  }
  PROCESS_SWITCH(EseTableProducer, processESE, "proccess q vectors to calculate reduced q-vector", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<EseTableProducer>(cfgc)}; }
