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

#include "FFitWeights.h"

#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EseTable.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TMath.h>

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

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

  Produces<o2::aod::MeanPt> meanPts;
  Produces<o2::aod::MeanPtShapes> meanPtShapes;

  OutputObj<FFitWeights> weightsFFit{FFitWeights("weights")};
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  Configurable<bool> cfgESE{"cfgESE", 1, "ese activation step: false = no ese, true = evaluate qSelection and fill table"};
  Configurable<int> cfgMeanPt{"cfgMeanPt", 0, "lvl, 0: First profile, 1: Second profile, 2: fill table"};
  Configurable<std::string> cfgEsePath{"cfgEsePath", "Users/j/joachiha/ESE/local/ffitsplines", "CCDB path for ese splines"};
  Configurable<std::vector<std::string>> cfgDetectors{"cfgDetectors", {"FT0C"}, "detectors to loop over: ['FT0C', 'FT0A', 'FV0A', 'TPCall', 'TPCneg', 'TPCpos']"};
  Configurable<std::vector<int>> cfgLoopHarmonics{"cfgLoopHarmonics", {2, 3}, "Harmonics to loop over when filling and evaluating q-Selection"};

  Configurable<std::vector<float>> cfgaxisqn{"cfgaxisqn", {500, 0, 25}, "q_n amplitude range"};
  Configurable<int> cfgnResolution{"cfgnResolution", 3000, "resolution of q-Selection"};

  Configurable<int> cfgnTotalSystem{"cfgnTotalSystem", 7, "total qvector number // look in Qvector table for this number"};
  Configurable<int> cfgnCorrLevel{"cfgnCorrLevel", 3, "QVector step: 0 = no corr, 1 = rect, 2 = twist, 3 = full"};

  int runNumber{-1};

  static constexpr float ThresholdAmplitude{1e-8f};
  static constexpr int Step1{1};
  static constexpr int Step2{2};

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

  FFitWeights* eventShape{nullptr};

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  struct Config {
    TH1D* mEfficiency = nullptr;
    bool correctionsLoaded = false;
  } cfg;

  Configurable<float> cfgVtxZ{"cfgVtxZ", 10.0f, "max z vertex position"};
  Configurable<float> cfgEta{"cfgEta", 0.8f, "max eta"};
  Configurable<float> cfgPtmin{"cfgPtmin", 0.2f, "min pt"};
  Configurable<float> cfgPtmax{"cfgPtmax", 5.0f, "max pt"};
  Configurable<float> cfgChi2PrITSCls{"cfgChi2PrITSCls", 4.0f, "max chi2 per ITS cluster"};
  Configurable<float> cfgChi2PrTPCCls{"cfgChi2PrTPCCls", 2.5f, "max chi2 per TPC cluster"};
  Configurable<float> cfgDCAz{"cfgDCAz", 2.0f, "max DCAz cut"};
  Configurable<std::string> cfgEfficiency{"cfgEfficiency", "", "CCDB path to efficiency object"};

  // o2::framework::expressions::Filter collisionFilter = nabs(aod::collision::posZ) < cfgVtxZ;
  o2::framework::expressions::Filter trackFilter = nabs(aod::track::eta) < cfgEta && aod::track::pt > cfgPtmin&& aod::track::pt < cfgPtmax && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == static_cast<uint8_t>(true))) && (aod::track::itsChi2NCl < cfgChi2PrITSCls) && (aod::track::tpcChi2NCl < cfgChi2PrTPCCls) && nabs(aod::track::dcaZ) < cfgDCAz;

  Preslice<aod::Tracks> perCollision = aod::track::collisionId;
  using GFWTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>>;

  void init(o2::framework::InitContext&)
  {

    LOGF(info, "ESETable::init()");

    registry.add("hEventCounter", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
    registry.add("hESEstat", "ese status;ese status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
    registry.add("hMeanPtStat", "", {HistType::kTH1F, {{4, 0.0, 4.0}}});

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
    weightsFFit->setBinAxis(cfgaxisqn->at(0), cfgaxisqn->at(1), cfgaxisqn->at(2));
    weightsFFit->setResolution(cfgnResolution);
    weightsFFit->setQnType(veccfg);
    weightsFFit->init();
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {

    auto timestamp = bc.timestamp();

    if (cfgESE) {
      eventShape = ccdb->getForTimeStamp<FFitWeights>(cfgEsePath, timestamp);
      if (!eventShape)
        LOGF(fatal, "failed loading qSelection with ese flag");
      LOGF(info, "successfully loaded qSelection");
    }
    if (!cfgEfficiency.value.empty()) {
      cfg.mEfficiency = ccdb->getForTimeStamp<TH1D>(cfgEfficiency, timestamp);
      if (cfg.mEfficiency == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram from %s", cfgEfficiency.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram from %s (%p)", cfgEfficiency.value.c_str(), (void*)cfg.mEfficiency);
      cfg.correctionsLoaded = true;
    }
  }

  float calcRedqn(const float& Qx, const float& Qy, const float& Mult)
  {
    float dqn{0.0f};
    float qn{0.0f};

    dqn = Qx * Qx + Qy * Qy;
    qn = std::sqrt(dqn) / std::sqrt(Mult);
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
    if (sumAmpl > ThresholdAmplitude) {
      float qnval = calcRedqn(QX * sumAmpl, QY * sumAmpl, sumAmpl);
      weightsFFit->fillWeights(centr, qnval, nHarm, pf);
      if (cfgESE) {
        splineVal = eventShape->eval(centr, qnval, nHarm, pf);
      }
    }
  }

  template <typename C>
  std::tuple<const float, const float, const float> getVectors(const C& col, const int& nHarm, const DetID& id)
  {
    const int detId = detIDN(id);
    const int detInd{detId * 4 + cfgnTotalSystem * 4 * (nHarm - 2)};
    const auto fQx{col.qvecRe()[detInd + cfgnCorrLevel]};
    const auto fQy{col.qvecIm()[detInd + cfgnCorrLevel]};
    const auto sumAmpl{col.qvecAmp()[detId]};
    return {fQx, fQy, sumAmpl};
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

  template <typename TTracks, typename Cent>
  std::pair<double, double> calculateMeanPt(TTracks const& tracks, Cent const& centrality)
  {
    double meanPtEvent{0.0};
    double effEvent{0.0};
    for (const auto& track : tracks) {
      double weff = getEfficiency(track);
      effEvent += weff;
      meanPtEvent += track.pt() * weff;
    }
    if (meanPtEvent == 0.0)
      return std::make_pair(0.0, 0.0);
    double mean = meanPtEvent / effEvent;
    weightsFFit->fillPt(centrality, mean, effEvent, true);
    return std::make_pair(mean, effEvent);
  }

  template <typename TTrack>
  double getEfficiency(TTrack track)
  {
    double eff = 1.;
    if (cfg.mEfficiency)
      eff = cfg.mEfficiency->GetBinContent(cfg.mEfficiency->FindBin(track.pt()));
    if (eff == 0)
      return -1.;
    else
      return 1. / eff;
  }

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
  PROCESS_SWITCH(EseTableProducer, processESE, "process q vectors to calculate reduced q-vector", true);

  void processMeanPt(soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::CentFT0Ms, aod::CentNTPVs, aod::CentNGlobals, aod::CentMFTs>::iterator const& collision, aod::BCsWithTimestamps const&, GFWTracks const& tracks)
  {

    std::vector<float> meanPt{-1};
    std::vector<float> meanPtShape{-1};
    if (collision.posZ() < -cfgVtxZ || collision.posZ() > cfgVtxZ) {
      meanPts(meanPt);
      meanPtShapes(meanPtShape);
      return;
    }

    registry.fill(HIST("hMeanPtStat"), 0.5);
    const auto centrality = collision.centFT0C();
    const auto mean = calculateMeanPt(tracks, centrality);

    if (cfgMeanPt == 0) {
      registry.fill(HIST("hMeanPtStat"), 1.5);
    } else {
      const auto avgpt = eventShape->getPtMult(centrality);
      if (mean.first == 0.0) {
        registry.fill(HIST("hMeanPtStat"), cfgMeanPt == Step1 ? 2.5 : 3.5);
      } else {
        const auto binval = (mean.first - avgpt) / avgpt;
        weightsFFit->fillPt(centrality, binval, mean.second, false);
        meanPt[0] = binval;
        if (cfgMeanPt == Step1) {
          registry.fill(HIST("hMeanPtStat"), 2.5);
        } else if (cfgMeanPt == Step2) {
          registry.fill(HIST("hMeanPtStat"), 3.5);
          const auto value = eventShape->evalPt(centrality, binval);
          meanPtShape[0] = value;
        }
      }
    }

    meanPts(meanPt);
    meanPtShapes(meanPtShape);
  }
  PROCESS_SWITCH(EseTableProducer, processMeanPt, "process mean pt selection", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<EseTableProducer>(cfgc)}; }
