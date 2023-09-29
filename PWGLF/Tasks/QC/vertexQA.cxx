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

#include <cmath>
#include <vector>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
DECLARE_SOA_TABLE(VtxQAtable, "AOD", "VTXQATABLE",
                  bc::GlobalBC,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  collision::CollisionTime,
                  collision::NumContrib)
}
struct vertexQA {
  Produces<o2::aod::VtxQAtable> vtxQAtable;

  Configurable<bool> storeTree{"storeTree", false, "Store tree for in-depth analysis"};

  ConfigurableAxis xVtxAxis{"xVtxBins", {100, -0.1f, 0.1f}, "Binning for the vertex x (y) in cm"};
  ConfigurableAxis zVtxAxis{"zVtxBins", {100, -20.f, 20.f}, "Binning for the vertex z in cm"};
  ConfigurableAxis tVtxAxis{"tVtxBins", {100, -50.f, 50.f}, "Binning for the vertex t in ns"};

  ConfigurableAxis xDiffVtxAxis{"xDiffVtxBins", {100, -0.1f, 0.1f}, "Binning for the vertex x (y) distance in cm"};
  ConfigurableAxis xyDiffVtxAxis{"xyDiffVtxBins", {100, 0.f, 0.1f}, "Binning for the vertex xy distance in cm"};
  ConfigurableAxis zDiffVtxAxis{"zDiffVtxBins", {100, -20.f, 20.f}, "Binning for the vertex z distance in cm"};
  ConfigurableAxis tDiffVtxAxis{"tDiffVtxBins", {100, -50.f, 50.f}, "Binning for the vertex t distance in ns"};

  ConfigurableAxis nVtxAxis{"nVtxBins", {11, -0.5, 10.5}, "Binning for the number of reconstructed vertices per BC"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

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
  }

  void process(aod::BC const& bc, aod::Collisions const& collisions)
  {
    auto collSize{collisions.size()};
    std::vector<double> collPosX;
    std::vector<double> collPosY;
    std::vector<double> collPosZ;
    std::vector<double> collPosT;
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
      collContribs.push_back(collision.numContrib());

      histos.fill(HIST("xVtxHistogram"), posX);
      histos.fill(HIST("yVtxHistogram"), posY);
      histos.fill(HIST("zVtxHistogram"), posZ);
      histos.fill(HIST("tVtxHistogram"), posT);
    }

    if (collSize == 2) {
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

    histos.fill(HIST("nVtxHistogram"), collSize);
    if (storeTree && collSize > 1) {
      for (int i{0}; i < collSize; ++i) {
        vtxQAtable(bc.globalBC(), collPosX[i], collPosY[i], collPosZ[i], collPosT[i], collContribs[i]);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<vertexQA>(cfgc)};
}
