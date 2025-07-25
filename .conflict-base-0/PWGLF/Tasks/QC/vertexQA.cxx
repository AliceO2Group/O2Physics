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

#include <algorithm>
#include <cmath>
#include <deque>
#include <string>
#include <vector>
#include <utility>

#include "CCDB/BasicCCDBManager.h"
#include "Common/CCDB/ctpRateFetcher.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using BCcoll = std::pair<aod::BC, aod::Collision>;

namespace
{
constexpr double LHCRFFreq = 400.789e6;
constexpr double LHCBunchSpacingNS = 10 * 1.e9 / LHCRFFreq;
double deltaTimeColl(BCcoll const bccoll1, BCcoll const bccoll2)
{
  auto coll1 = std::get<aod::Collision>(bccoll1);
  auto coll2 = std::get<aod::Collision>(bccoll2);
  auto bc1 = std::get<aod::BC>(bccoll1);
  auto bc2 = std::get<aod::BC>(bccoll2);
  int64_t tmpDT = static_cast<int64_t>(bc1.globalBC()) - static_cast<int64_t>(bc2.globalBC());
  double deltaT = tmpDT * LHCBunchSpacingNS + coll1.collisionTime() - coll2.collisionTime();
  return deltaT;
}
} // namespace

namespace o2::aod
{
DECLARE_SOA_TABLE(VtxQAtable, "AOD", "VTXQATABLE",
                  bc::GlobalBC,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  collision::CovXX,
                  collision::CovXY,
                  collision::CovYY,
                  collision::CovZZ,
                  collision::CollisionTime,
                  collision::CollisionTimeRes,
                  collision::NumContrib)
} // namespace o2::aod

struct vertexQA {
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  ctpRateFetcher mRateFetcher;

  Produces<o2::aod::VtxQAtable> vtxQAtable;

  Configurable<int> storeTree{"storeTree", 1000, "Store in tree collisions from BC's with more than 'storeTree' vertices, for in-depth analysis"};

  Configurable<uint64_t> nCollMax{"nCollMax", 20, "Maximum size of collision buffer"};
  Configurable<double> nSigmaZ{"nSigmaZ", 1000., "Number of sigmas for z of vertices"};
  Configurable<double> nSigmaR{"nSigmaR", 1000., "Number of sigmas for transverse displacement of vertices"};
  Configurable<double> nSigmaT{"nSigmaT", 1000., "Number of sigmas for time of vertices"};
  Configurable<double> maxTime{"maxTime", 100000., "Maximum time difference between split vertices in ns"};
  Configurable<double> minTimeDiff{"minTimeDiff", 14000., "Minimum time difference (ITS Rof peak)"};
  Configurable<double> maxTimeDiff{"maxTimeDiff", 15500., "Maximum time difference (ITS Rof peak)"};
  Configurable<double> minNcontrib{"minNcontrib", 500., "Minimum number of contributors (duplicate vertices)"};

  ConfigurableAxis xVtxAxis{"xVtxBins", {100, -0.1f, 0.1f}, "Binning for the vertex x (y) in cm"};
  ConfigurableAxis zVtxAxis{"zVtxBins", {100, -20.f, 20.f}, "Binning for the vertex z in cm"};
  ConfigurableAxis tVtxAxis{"tVtxBins", {100, -20.f, 20.f}, "Binning for the vertex t in ns"};

  ConfigurableAxis xDiffVtxAxis{"xDiffVtxBins", {200, -0.05f, 0.05f}, "Binning for the vertex x (y) distance in cm"};
  ConfigurableAxis xyDiffVtxAxis{"xyDiffVtxBins", {200, 0.f, 0.05f}, "Binning for the vertex xy distance in cm"};
  ConfigurableAxis zDiffVtxAxis{"zDiffVtxBins", {200, -10.f, 10.f}, "Binning for the vertex z distance in cm"};
  ConfigurableAxis zDiffVtxAxis2{"zDiffVtxBins2", {1000, -0.5f, 0.5f}, "Binning for the vertex z distance in cm (2)"};
  ConfigurableAxis tDiffVtxAxis{"tDiffVtxBins", {300, 0.f, 30.f}, "Binning for the vertex t distance in ns"};
  ConfigurableAxis tDiffVtxAxisExtend{"tDiffVtxBinsExtend", {1000, 0.f, 100000.f}, "Binning for the vertex t distance in ns, extended range"};
  ConfigurableAxis tDiffVtxAxisExtendSigned{"tDiffVtxBinsExtendSigned", {1000, -100000.f, 100000.f}, "Binning for the vertex t distance in ns, extended range, signed"};

  ConfigurableAxis nVtxAxis{"nVtxBins", {11, -0.5, 10.5}, "Binning for the number of reconstructed vertices per BC"};
  ConfigurableAxis nVtxTwoAxis{"nVtxTwoBins", {2, 0.5, 2.5}, "Binning for the number of reconstructed vertices vs. time"};
  ConfigurableAxis nContribAxis{"nContribBins", {1000, 0, 5000}, "Binning for number of contributors to PV"};
  ConfigurableAxis nContribDiffAxis{"nContribDiffBins", {1000, -5000, 5000}, "Binning for the difference in number of contributors to PV"};

  ConfigurableAxis irBinning{"IRbinning", {500, 0, 100}, "Binning for the interaction rate (kHz)"};
  Configurable<std::string> irSource{"irSource", "ZNC hadronic", "Source of the interaction rate"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  std::deque<BCcoll> colls;

  int64_t mFirstBCid = -1;

  void init(InitContext const&)
  {
    histos.add<TH1>("nVtxHistogram", ";#it{N}_{vtx}^{rec};Entries", HistType::kTH1F, {nVtxAxis});

    histos.add<TH1>("xVtxHistogram", ";#it{x}_{vtx} (cm);Entries", HistType::kTH1F, {xVtxAxis});
    histos.add<TH1>("yVtxHistogram", ";#it{y}_{vtx} (cm);Entries", HistType::kTH1F, {xVtxAxis});
    histos.add<TH1>("zVtxHistogram", ";#it{z}_{vtx} (cm);Entries", HistType::kTH1F, {zVtxAxis});
    histos.add<TH1>("tVtxHistogram", ";#it{t}_{vtx} (ns);Entries", HistType::kTH1F, {tVtxAxis});

    histos.add<TH1>("xDistVtxHistogram", ";#Delta#it{x}_{vtx} (cm);Entries", HistType::kTH1F, {xDiffVtxAxis});
    histos.add<TH1>("yDistVtxHistogram", ";#Delta#it{y}_{vtx} (cm);Entries", HistType::kTH1F, {xDiffVtxAxis});
    histos.add<TH1>("xyDistVtxHistogram", ";#Delta#it{xy}_{vtx} (cm);Entries", HistType::kTH1F, {xyDiffVtxAxis});
    histos.add<TH1>("zDistVtxHistogram", ";#Delta#it{z}_{vtx} (cm);Entries", HistType::kTH1F, {zDiffVtxAxis});
    histos.add<TH1>("tDistVtxHistogram", ";#Delta#it{t}_{vtx} (ns);Entries", HistType::kTH1F, {tDiffVtxAxis});

    histos.add<TH2>("xTwoVtxHistogram", ";#it{x}_{vtx}^{1} (cm);#it{x}_{vtx}^{2} (cm)", HistType::kTH2F, {xVtxAxis, xVtxAxis});
    histos.add<TH2>("yTwoVtxHistogram", ";#it{y}_{vtx}^{1} (cm);#it{y}_{vtx}^{2} (cm)", HistType::kTH2F, {xVtxAxis, xVtxAxis});
    histos.add<TH2>("zTwoVtxHistogram", ";#it{z}_{vtx}^{1} (cm);#it{z}_{vtx}^{2} (cm)", HistType::kTH2F, {zVtxAxis, zVtxAxis});
    histos.add<TH2>("tTwoVtxHistogram", ";#it{t}_{vtx}^{1} (ns);#it{t}_{vtx}^{2} (ns)", HistType::kTH2F, {tVtxAxis, tVtxAxis});

    histos.add<TH1>("nVtxTimeSeriesHistogram", ";#it{N}_{vtx}^{rec};Entries", HistType::kTH1F, {nVtxTwoAxis});
    histos.add<TH1>("zDistVtxTimeSeriesHistogram", ";#Delta#it{z}_{vtx} (cm);Entries", HistType::kTH1F, {zDiffVtxAxis2});
    histos.add<TH1>("tDistVtxTimeSeriesHistogram", ";#Delta#it{t}_{vtx} (ns);Entries", HistType::kTH1F, {tDiffVtxAxisExtend});
    histos.add<TH2>("tCollTwoVtxTimeSeriesHistogram", ";#Delta#it{t}_{coll}^{1} (ns);#Delta#it{t}_{coll}^{2} (ns)", HistType::kTH2F, {tVtxAxis, tVtxAxis});
    histos.add<TH2>("zDistVsTDistVtxTimeSeriesHistogram", ";#Delta#it{t}_{vtx} (ns);#Delta#it{z}_{vtx} (cm)", HistType::kTH2F, {tDiffVtxAxisExtend, zDiffVtxAxis2});
    histos.add<TH2>("nContribTwoVtxTimeSeriesHistogram", ";#it{N}_{contrib}^{1};#it{N}_{contrib}^{2}", HistType::kTH2F, {nContribAxis, nContribAxis});
    histos.add<TH2>("nContribVsTDistTimeSeriesHistogram", ";#Delta#it{t}_{vtx} (ns);#Delta#it{N}_{contrib}", HistType::kTH2F, {tDiffVtxAxisExtendSigned, nContribDiffAxis});

    histos.add<TH2>("nContribITSRofTimeSeriesHistogram", ";#it{N}_{contrib}^{1};#it{N}_{contrib}^{2}", HistType::kTH2F, {nContribAxis, nContribAxis});
    histos.add<TH1>("tDiffDuplicateTimeSeriesHistogram", ";#Delta#it{t}_{vtx} (ns);Entries", HistType::kTH1F, {tDiffVtxAxisExtend});
    histos.add<TH2>("tIRvsCollisionRateHistogram", Form(";IR from %s (kHz);IR from reconstructed vertices (kHz)", irSource.value.data()), HistType::kTH2D, {irBinning, irBinning});

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setFatalWhenNull(false);
  }

  void process(aod::BC const& bc, aod::Collisions const& collisions)
  {
    auto collSize{collisions.size()};
    std::vector<double> collPosX;
    std::vector<double> collPosY;
    std::vector<double> collPosZ;
    std::vector<double> collPosT;
    std::vector<double> collCovXX;
    std::vector<double> collCovXY;
    std::vector<double> collCovYY;
    std::vector<double> collCovZZ;
    std::vector<double> collResT;
    std::vector<uint16_t> collContribs;

    for (auto& collision : collisions) {
      auto posX{collision.posX()};
      auto posY{collision.posY()};
      auto posZ{collision.posZ()};
      auto posT{collision.collisionTime()};

      collPosX.push_back(posX);
      collPosY.push_back(posY);
      collPosZ.push_back(posZ);
      collPosT.push_back(posT);
      collCovXX.push_back(collision.covXX());
      collCovXY.push_back(collision.covXY());
      collCovYY.push_back(collision.covYY());
      collCovZZ.push_back(collision.covZZ());
      collResT.push_back(collision.collisionTimeRes());
      collContribs.push_back(collision.numContrib());

      histos.fill(HIST("xVtxHistogram"), posX);
      histos.fill(HIST("yVtxHistogram"), posY);
      histos.fill(HIST("zVtxHistogram"), posZ);
      histos.fill(HIST("tVtxHistogram"), posT);

      auto compareCollisions = [&](BCcoll other) {
        auto coll2 = std::get<aod::Collision>(other);
        double testZ = std::abs(collision.posZ() - coll2.posZ()) / std::sqrt(collision.covZZ() + coll2.covZZ());
        double distR = std::hypot(collision.posX() - coll2.posX(), collision.posY() - coll2.posY());
        double drdxCol = (collision.posX() - coll2.posX()) / distR;
        double drdxOther = -(collision.posX() - coll2.posX()) / distR;
        double drdyCol = (collision.posY() - coll2.posY()) / distR;
        double drdyOther = -(collision.posY() - coll2.posY()) / distR;
        double covR = std::pow(drdxCol, 2) * collision.covXX() + std::pow(drdxOther, 2) * coll2.covXX() + std::pow(drdyCol, 2) * collision.covYY() + std::pow(drdxOther, 2) * coll2.covYY() + 2 * drdxCol * drdyCol * collision.covXY() + 2 * drdxOther * drdyOther * coll2.covXY();
        double testR = distR / std::sqrt(covR);
        double deltaT = deltaTimeColl(BCcoll{bc, collision}, other);
        double testT = std::abs(deltaT) / std::sqrt(std::pow(collision.collisionTimeRes(), 2) + std::pow(coll2.collisionTimeRes(), 2));
        return (testT < nSigmaT && testZ < nSigmaZ && testR < nSigmaR && std::abs(deltaT) < maxTime);
      };

      if (colls.size() > 0) {
        auto id = std::find_if(colls.begin(), colls.end(), compareCollisions);
        if (id != colls.end()) {
          auto coll2 = std::get<aod::Collision>(*id);
          double deltaT = deltaTimeColl(BCcoll{bc, collision}, *id);
          // Keep into account that we have "past" collisions in the queue
          histos.fill(HIST("zDistVtxTimeSeriesHistogram"), collision.posZ() - coll2.posZ());
          histos.fill(HIST("tDistVtxTimeSeriesHistogram"), std::abs(deltaT));
          histos.fill(HIST("tCollTwoVtxTimeSeriesHistogram"), coll2.collisionTime(), collision.collisionTime());
          histos.fill(HIST("zDistVsTDistVtxTimeSeriesHistogram"), std::abs(deltaT), collision.posZ() - coll2.posZ());
          histos.fill(HIST("nContribTwoVtxTimeSeriesHistogram"), coll2.numContrib(), collision.numContrib());
          histos.fill(HIST("nContribVsTDistTimeSeriesHistogram"), deltaT, collision.numContrib() - coll2.numContrib());
          histos.fill(HIST("nVtxTimeSeriesHistogram"), 2);
          if (std::abs(deltaT) > minTimeDiff && std::abs(deltaT) < maxTimeDiff) {
            histos.fill(HIST("nContribITSRofTimeSeriesHistogram"), coll2.numContrib(), collision.numContrib());
          }
          if (std::abs(collision.numContrib() - coll2.numContrib()) < 3 * std::sqrt(collision.numContrib() + coll2.numContrib()) && collision.numContrib() > minNcontrib && coll2.numContrib() > minNcontrib) {
            histos.fill(HIST("tDiffDuplicateTimeSeriesHistogram"), std::abs(deltaT));
          }
          colls.erase(id);
        } else {
          colls.emplace_back(bc, collision);
        }
      } else {
        colls.emplace_back(bc, collision);
      }
    }

    histos.fill(HIST("nVtxHistogram"), collSize);

    if (collSize >= 2) {
      histos.fill(HIST("xDistVtxHistogram"), collPosX[0] - collPosX[1]);
      histos.fill(HIST("yDistVtxHistogram"), collPosY[0] - collPosY[1]);
      histos.fill(HIST("xyDistVtxHistogram"), std::hypot(collPosX[0] - collPosX[1], collPosY[0] - collPosY[1]));
      histos.fill(HIST("zDistVtxHistogram"), collPosZ[0] - collPosZ[1]);
      histos.fill(HIST("tDistVtxHistogram"), collPosT[0] - collPosT[1]);

      histos.fill(HIST("xTwoVtxHistogram"), collPosX[0], collPosX[1]);
      histos.fill(HIST("yTwoVtxHistogram"), collPosY[0], collPosY[1]);
      histos.fill(HIST("zTwoVtxHistogram"), collPosZ[0], collPosZ[1]);
      histos.fill(HIST("tTwoVtxHistogram"), collPosT[0], collPosT[1]);
    }

    if (colls.size() > nCollMax) {
      auto delta = colls.size() - nCollMax;
      for (uint64_t iC{0}; iC < delta; ++iC) {
        histos.fill(HIST("nVtxTimeSeriesHistogram"), 1);
        colls.pop_front();
      }
    }

    if (collSize > storeTree) {
      for (int i{0}; i < collSize; ++i) {
        vtxQAtable(bc.globalBC(), collPosX[i], collPosY[i], collPosZ[i], collCovXX[i], collCovXY[i], collCovYY[i], collCovZZ[i], collPosT[i], collResT[i], collContribs[i]);
      }
    }
  }
  PROCESS_SWITCH(vertexQA, process, "Standard vertex QA", true);

  void processIR(aod::BCsWithTimestamps const& bcs, aod::Collisions const& collisions)
  {
    if (collisions.size() <= 2) {
      return;
    }

    std::vector<int64_t> jumps{0ll};
    int64_t lastBC = bcs.rawIteratorAt(0).globalBC();
    for (auto bc : bcs) {
      if (bc.globalBC() - lastBC > 3564 * 32) { // 32 orbits
        jumps.push_back(bc.globalIndex());
        lastBC = bc.globalBC();
      }
    }
    uint64_t jumpsSentinel{1};
    std::vector<int64_t> collisionsIndices{0ll};
    for (auto col : collisions) {
      if (jumpsSentinel == jumps.size()) {
        break;
      }
      if (col.bcId() > jumps[jumpsSentinel]) {
        collisionsIndices.push_back(col.globalIndex());
        jumpsSentinel++;
      }
    }
    jumps.push_back(bcs.size());
    collisionsIndices.push_back(collisions.size());

    for (size_t i{0}; i < jumps.size() - 1; ++i) {
      auto startBC = bcs.rawIteratorAt(jumps[i]);
      auto endBC = bcs.rawIteratorAt(jumps[i + 1] - 1);
      double startIR = mRateFetcher.fetch(ccdb.service, startBC.timestamp(), startBC.runNumber(), irSource.value);
      double endIR = mRateFetcher.fetch(ccdb.service, endBC.timestamp(), endBC.runNumber(), irSource.value);
      double deltaT = (endBC.globalBC() - startBC.globalBC()) * LHCBunchSpacingNS * 1.e-9;
      double collisionRate = (collisionsIndices[i + 1] - collisionsIndices[i]) / deltaT; /// -1 to remove the bias of the collisions at extremities?
      double ir = (startIR + endIR) * 0.5;
      histos.fill(HIST("tIRvsCollisionRateHistogram"), ir * 1.e-3, collisionRate * 1.e-3);
    }
  }
  PROCESS_SWITCH(vertexQA, processIR, "Checks on interaction rate", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<vertexQA>(cfgc)};
}
