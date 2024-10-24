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

using CollWithMults = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::QvectorFT0CVecs>;

struct EseTableProducer {
  Produces<o2::aod::QPercentileFT0Cs> qPercsFT0C;
  Produces<o2::aod::FEseCols> fEseCol;

  OutputObj<FFitWeights> FFitObj{FFitWeights("weights")};
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  Configurable<bool> cfgESE{"cfgESE", 1, "ese activation step: false = no ese, true = evaluate splines and fill table"};
  Configurable<std::string> cfgEsePath{"cfgEsePath", "Users/j/joachiha/ESE/local/splines", "CCDB path for ese splines"};
  Configurable<bool> cfgFT0C{"cfgFT0C", 1, "FT0C flag"};
  Configurable<std::vector<int>> cfgLoopHarmonics{"cfgLoopHarmonics", {2, 3}, "Harmonics to loop over when filling and evaluating splines"};
  Configurable<std::vector<std::string>> vecStr{"vecStr", {"FT0C"}, "Qn type"};

  Configurable<std::vector<float>> cfgaxisqnFT0C{"cfgaxisqnFT0C", {500, 0, 25}, "q_n amplitude range"};
  Configurable<int> cfgnResolution{"cfgnResolution", 3000, "resolution of splines"};

  int runNumber{-1};

  FFitWeights* splines{nullptr};

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  void init(o2::framework::InitContext&)
  {

    LOGF(info, "ESETable::init()");

    registry.add("hEventCounter", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
    registry.add("hESEstat", "ese status;ese status;entries", {HistType::kTH1F, {{2, 0.0, 2.0}}});

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    std::vector<std::pair<int, std::string>> veccfg;
    for (std::size_t i{0}; i < cfgLoopHarmonics->size(); i++) 
    {
      for (std::size_t j{0}; j < vecStr->size(); j++) 
      {
        veccfg.push_back({cfgLoopHarmonics->at(i), vecStr->at(j)});
      }
    }
    FFitObj->SetBinAxis(cfgaxisqnFT0C->at(0), cfgaxisqnFT0C->at(1), cfgaxisqnFT0C->at(2));
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

  template <typename T>
  void calculateESE(T const& collision, std::vector<float>& qnpFT0C, std::vector<int>& fIsEseAvailable)
  {

    if (cfgFT0C) {
      const auto QxFT0C_Qvec = collision.qvecFT0CReVec();
      const auto QyFT0C_Qvec = collision.qvecFT0CImVec();
      const auto SumAmplFT0C = collision.sumAmplFT0C();
      for (std::size_t i{0}; i < cfgLoopHarmonics->size(); i++) {
        int nHarm = cfgLoopHarmonics->at(i);
        std::cout<< "nharm " << nHarm << std::endl;
        if (validQvec(QxFT0C_Qvec[nHarm - 2]) && validQvec(QyFT0C_Qvec[nHarm - 2]) && SumAmplFT0C > 1e-8) {
          float qnval = Calcqn(QxFT0C_Qvec[nHarm - 2]*SumAmplFT0C, QyFT0C_Qvec[nHarm - 2]*SumAmplFT0C, SumAmplFT0C);
          FFitObj->Fill(collision.centFT0C(), qnval, nHarm, "FT0C");
          if (cfgESE) {
            float qnCent{splines->EvalSplines(collision.centFT0C(), qnval, nHarm, "FT0C")};
            qnpFT0C.push_back(qnCent);
            fIsEseAvailable.push_back(1);
            registry.fill(HIST("hESEstat"), 1.5);
          } else {
            qnpFT0C.push_back(-1);
            fIsEseAvailable.push_back(0);
            registry.fill(HIST("hESEstat"), .5);
          }
        } else {
          qnpFT0C.push_back(-1);
          fIsEseAvailable.push_back(0);
          registry.fill(HIST("hESEstat"), .5);
        }
      }
    } else {
      for (std::size_t i{0}; i < cfgLoopHarmonics->size(); i++) {
        qnpFT0C.push_back(-1);
        fIsEseAvailable.push_back(0);
        registry.fill(HIST("hESEstat"), .5);
      }
    }
  };

  void processESE(CollWithMults::iterator const& collision, aod::BCsWithTimestamps const&, aod::FV0As const&, aod::FT0s const&)
  {
    float counter{0.5};
    registry.fill(HIST("hEventCounter"), counter++);

    std::vector<float> qnpFT0C{};
    std::vector<int> fIsEseAvailable{};

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int currentRun = bc.runNumber();
    if (runNumber != currentRun) {
      runNumber = currentRun;
      initCCDB(bc);
    }
    registry.fill(HIST("hEventCounter"), counter++);
    calculateESE(collision, qnpFT0C, fIsEseAvailable);

    qPercsFT0C(qnpFT0C);
    fEseCol(fIsEseAvailable);
    registry.fill(HIST("hEventCounter"), counter++);
  }
  PROCESS_SWITCH(EseTableProducer, processESE, "proccess q vectors to calculate reduced q-vector", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<EseTableProducer>(cfgc)}; }
