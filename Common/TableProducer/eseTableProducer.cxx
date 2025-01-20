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
#include <tuple>
#include <unordered_map>

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

using CollWithMults = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs, aod::Qvectors>;

struct EseTableProducer {
  Produces<o2::aod::QPercentileFT0Cs> qPercsFT0C;
  Produces<o2::aod::QPercentileFT0As> qPercsFT0A;
  Produces<o2::aod::QPercentileFV0As> qPercsFV0A;
  Produces<o2::aod::QPercentileTPCalls> qPercsTPCall;
  Produces<o2::aod::QPercentileTPCnegs> qPercsTPCneg;
  Produces<o2::aod::QPercentileTPCposs> qPercsTPCpos;

  OutputObj<FFitWeights> FFitObj{FFitWeights("weights")};
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  Configurable<bool> cfgESE{"cfgESE", 1, "ese activation step: false = no ese, true = evaluate qSelection and fill table"};
  Configurable<std::string> cfgEsePath{"cfgEsePath", "Users/j/joachiha/ESE/local/ffitsplines", "CCDB path for ese splines"};
  Configurable<std::vector<std::string>> cfgDetectors{"cfgDetectors", {"FT0C"}, "detectors to loop over: ['FT0C', 'FT0A', 'FV0A', 'TPCall', 'TPCneg', 'TPCpos']"};
  Configurable<std::vector<int>> cfgLoopHarmonics{"cfgLoopHarmonics", {2, 3}, "Harmonics to loop over when filling and evaluating q-Selection"};

  Configurable<std::vector<float>> cfgaxisqn{"cfgaxisqn", {500, 0, 25}, "q_n amplitude range"};
  Configurable<int> cfgnResolution{"cfgnResolution", 3000, "resolution of q-Selection"};

  Configurable<int> cfgnTotalSystem{"cfgnTotalSystem", 7, "total qvector number // look in Qvector table for this number"};
  Configurable<int> cfgnCorrLevel{"cfgnCorrLevel", 3, "QVector step: 0 = no corr, 1 = rect, 2 = twist, 3 = full"};

  int runNumber{-1};

  enum class DetID { FT0C,
                     FT0A,
                     FT0M,
                     FV0A,
                     TPCpos,
                     TPCneg,
                     TPCall };

  std::unordered_map<std::string, DetID> detMap = {
    {"FT0C", DetID::FT0C},
    {"FT0A", DetID::FT0A},
    {"FT0M", DetID::FT0M},
    {"FV0A", DetID::FV0A},
    {"TPCpos", DetID::TPCpos},
    {"TPCneg", DetID::TPCneg},
    {"TPCall", DetID::TPCall}};

  FFitWeights* qSelection{nullptr};

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

    std::vector<std::pair<int, std::string>> veccfg;
    for (std::size_t i{0}; i < cfgLoopHarmonics->size(); i++) {
      for (std::size_t j{0}; j < cfgDetectors->size(); j++) {
        veccfg.push_back({cfgLoopHarmonics->at(i), cfgDetectors->at(j)});
      }
    }
    FFitObj->setBinAxis(cfgaxisqn->at(0), cfgaxisqn->at(1), cfgaxisqn->at(2));
    FFitObj->setResolution(cfgnResolution);
    FFitObj->setQnType(veccfg);
    FFitObj->init();
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {

    auto timestamp = bc.timestamp();

    if (cfgESE) {
      qSelection = ccdb->getForTimeStamp<FFitWeights>(cfgEsePath, timestamp);
      if (!qSelection)
        LOGF(fatal, "failed loading qSelection with ese flag");
      LOGF(info, "successfully loaded qSelection");
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

  constexpr int detIDN(const DetID id)
  {
    switch (id) {
      case DetID::FT0C:
        return 0;
      case DetID::FT0A:
        return 1;
      case DetID::FT0M:
        return 2;
      case DetID::FV0A:
        return 3;
      case DetID::TPCpos:
        return 4;
      case DetID::TPCneg:
        return 5;
      case DetID::TPCall:
        return 6;
    }
    return -1;
  }

  void doSpline(float& splineVal, const float& centr, const float& nHarm, const char* pf, const auto& QX, const auto& QY, const auto& sumAmpl)
  {
    if (sumAmpl > 1e-8) {
      float qnval = Calcqn(QX * sumAmpl, QY * sumAmpl, sumAmpl);
      FFitObj->fillWeights(centr, qnval, nHarm, pf);
      if (cfgESE) {
        splineVal = qSelection->eval(centr, qnval, nHarm, pf);
      }
    }
  }

  template <typename C>
  std::tuple<const float, const float, const float> getVectors(const C& col, const int& nHarm, const DetID& id)
  {
    const int detId = detIDN(id);
    const int detInd{detId * 4 + cfgnTotalSystem * 4 * (nHarm - 2)};
    const auto Qx{col.qvecRe()[detInd + cfgnCorrLevel]};
    const auto Qy{col.qvecIm()[detInd + cfgnCorrLevel]};
    const auto sumAmpl{col.qvecAmp()[detId]};
    return {Qx, Qy, sumAmpl};
  }

  template <typename T>
  void calculateESE(T const& collision,
                    std::vector<float>& qnpFT0C,
                    std::vector<float>& qnpFT0A,
                    std::vector<float>& qnpFV0A,
                    std::vector<float>& qnpTPCall,
                    std::vector<float>& qnpTPCneg,
                    std::vector<float>& qnpTPCpos)
  {
    const float centrality = collision.centFT0C();
    float counter{0.5};
    registry.fill(HIST("hESEstat"), counter++);

    std::unordered_map<std::string, std::vector<float>*> vMap{
      {"FT0C", &qnpFT0C},
      {"FT0A", &qnpFT0A},
      {"FV0A", &qnpFV0A},
      {"TPCall", &qnpTPCall},
      {"TPCneg", &qnpTPCneg},
      {"TPCpos", &qnpTPCpos}};

    for (std::size_t j{0}; j < cfgDetectors->size(); j++) {
      const auto det{cfgDetectors->at(j)};
      const auto iter{detMap.find(det)};
      float splineVal{-1.0};

      if (iter != detMap.end()) {
        for (std::size_t i{0}; i < cfgLoopHarmonics->size(); i++) {
          const int nHarm{cfgLoopHarmonics->at(i)};
          const auto [qxt, qyt, st] = getVectors(collision, nHarm, iter->second);
          doSpline(splineVal, centrality, nHarm, det.c_str(), qxt, qyt, st);
          if (i == 0)
            registry.fill(HIST("hESEstat"), counter++);

          if (vMap.find(det) != vMap.end()) {
            vMap[det]->push_back(splineVal);
          }
        }
      }
    }
  };

  void processESE(CollWithMults::iterator const& collision, aod::BCsWithTimestamps const&, aod::FV0As const&, aod::FT0s const&)
  {
    float counter{0.5};
    registry.fill(HIST("hEventCounter"), counter++);

    std::vector<float> qnpFT0C{};
    std::vector<float> qnpFT0A{};
    std::vector<float> qnpFV0A{};
    std::vector<float> qnpTPCall{};
    std::vector<float> qnpTPCneg{};
    std::vector<float> qnpTPCpos{};

    auto bc{collision.bc_as<aod::BCsWithTimestamps>()};
    int currentRun{bc.runNumber()};
    if (runNumber != currentRun) {
      runNumber = currentRun;
      initCCDB(bc);
    }
    registry.fill(HIST("hEventCounter"), counter++);
    calculateESE(collision, qnpFT0C, qnpFT0A, qnpFV0A, qnpTPCall, qnpTPCneg, qnpTPCpos);

    qPercsFT0C(qnpFT0C);
    qPercsFT0A(qnpFT0A);
    qPercsFV0A(qnpFV0A);
    qPercsTPCall(qnpTPCall);
    qPercsTPCneg(qnpTPCneg);
    qPercsTPCpos(qnpTPCpos);
    registry.fill(HIST("hEventCounter"), counter++);
  }
  PROCESS_SWITCH(EseTableProducer, processESE, "proccess q vectors to calculate reduced q-vector", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<EseTableProducer>(cfgc)}; }
