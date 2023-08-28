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
/// \author Junlee Kim (jikim1290@gmail.com)
/// \since November 2021
// usage: o2-analysis-timestamp -b --aod-file AO2D.root --configuration
// json://./config.json | o2-analysis-trackextension -b |
// o2-analysis-trackselection -b --isRun3 0 | o2-analysis-mm-lumi -b
// --configuration json://./config.json

#include <array>
#include <cmath>
#include <iostream>
#include <vector>

#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonUtils/NameConf.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"

#include "DetectorsVertexing/PVertexer.h"

#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/PrimaryVertex.h"
#include "ReconstructionDataFormats/Vertex.h"

#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"

#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"

#include "CommonConstants/GeomConstants.h"

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"

#include "DataFormatsCalibration/MeanVertexObject.h"

namespace o2::aod
{
namespace full
{
DECLARE_SOA_COLUMN(TimeStamp, timeStamp, uint64_t);
DECLARE_SOA_COLUMN(VertexX, vertexX, double);
DECLARE_SOA_COLUMN(VertexY, vertexY, double);
DECLARE_SOA_COLUMN(VertexZ, vertexZ, double);

DECLARE_SOA_COLUMN(VertexXX, vertexXX, double);
DECLARE_SOA_COLUMN(VertexYY, vertexYY, double);
DECLARE_SOA_COLUMN(VertexXY, vertexXY, double);

DECLARE_SOA_COLUMN(VertexChi2, vertexChi2, double);
DECLARE_SOA_COLUMN(NContrib, nContrib, int);
} // namespace full
DECLARE_SOA_TABLE(EventInfo, "AOD", "EventInfo", full::TimeStamp, full::VertexX,
                  full::VertexY, full::VertexZ,

                  full::VertexXX, full::VertexYY, full::VertexXY,

                  full::VertexChi2, full::NContrib);
} // namespace o2::aod

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct lumiTask {
  Produces<o2::aod::EventInfo> rowEventInfo;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  const char* ccdbpath_grp = "GLO/Config/GRPMagField";
  const char* ccdburl = "http://alice-ccdb.cern.ch";
  int mRunNumber;

  Configurable<uint64_t> ftts{"ftts", 1530319778000,
                              "First time of time stamp"};
  Configurable<int> nContribMax{"nContribMax", 2500,
                                "Maximum number of contributors"};
  Configurable<int> nContribMin{"nContribMin", 10,
                                "Minimum number of contributors"};

  HistogramRegistry histos{
    "histos",
    {
      {"vertexx", "", {HistType::kTH1F, {{1000, -1, 1, "x"}}}},     //
      {"vertexy", "", {HistType::kTH1F, {{1000, -1, 1, "y"}}}},     //
      {"timestamp", "", {HistType::kTH1F, {{20000, 0, 2e7, "t"}}}}, //
      {"vertexx_timestamp",
       "",
       {HistType::kTH2F, {{1000, 0, 4e6, "t"}, {2000, -1, 1, "x"}}}}, //
      {"vertexy_timestamp",
       "",
       {HistType::kTH2F, {{1000, 0, 4e6, "t"}, {2000, -1, 1, "y"}}}},     //
      {"chisquare", "", {HistType::kTH1F, {{1000, 0, 100, "#chi^{2}"}}}}, //

      {"vertexx_Refitted", "", {HistType::kTH1F, {{1000, -1, 1, "x"}}}}, //
      {"vertexy_Refitted", "", {HistType::kTH1F, {{1000, -1, 1, "y"}}}}, //
      {"vertexx_Refitted_timestamp",
       "",
       {HistType::kTH2F, {{1000, 0, 4e6, "t"}, {2000, -1, 1, "x"}}}}, //
      {"vertexy_Refitted_timestamp",
       "",
       {HistType::kTH2F, {{1000, 0, 4e6, "t"}, {2000, -1, 1, "y"}}}}, //
      {"chisquare_Refitted",
       "",
       {HistType::kTH1F, {{1000, 0, 100, "#chi^{2}"}}}}, //

      {"vertexx_Refitted_vertexx",
       "",
       {HistType::kTH2F, {{1000, -1, 1, "x"}, {1000, -1, 1, "rx"}}}}, //
      {"vertexy_Refitted_vertexy",
       "",
       {HistType::kTH2F, {{1000, -1, 1, "y"}, {1000, -1, 1, "ry"}}}} //
    }};
  bool doPVrefit = true;

  void init(InitContext&)
  {
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    uint64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);
    mRunNumber = 0;
  }

  void process(aod::Collision const& collision, aod::BCsWithTimestamps const&,
               o2::soa::Join<o2::aod::Tracks, o2::aod::TrackSelection,
                             o2::aod::TracksCov, o2::aod::TracksExtra,
                             o2::aod::TracksDCA> const& tracks,
               o2::soa::Join<o2::aod::Tracks, o2::aod::TracksCov,
                             o2::aod::TracksExtra> const& unfiltered_tracks)
  {

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    uint64_t relTS = bc.timestamp() - ftts;

    std::vector<int64_t> vec_globID_contr = {};
    std::vector<o2::track::TrackParCov> vec_TrkContributos = {};

    int nContrib = 0;
    int nNonContrib = 0;
    for (const auto& unfiltered_track : unfiltered_tracks) {
      if (!unfiltered_track.hasITS()) {
        nNonContrib++;
        continue;
      }
      if (unfiltered_track.pt() < 0.8 || unfiltered_track.itsNCls() < 5) {
        nNonContrib++;
        continue;
      }
      vec_globID_contr.push_back(unfiltered_track.globalIndex());
      vec_TrkContributos.push_back(getTrackParCov(unfiltered_track));
      nContrib++;
    }

    std::vector<bool> vec_useTrk_PVrefit(vec_globID_contr.size(), true);

    if (mRunNumber != bc.runNumber()) {
      o2::parameters::GRPMagField* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(ccdbpath_grp, bc.timestamp());
      if (grpo != nullptr) {
        o2::base::Propagator::initFieldFromGRP(grpo);
      } else {
        LOGF(fatal,
             "GRP object is not available in CCDB for run=%d at timestamp=%llu",
             bc.runNumber(), bc.timestamp());
      }
      mRunNumber = bc.runNumber();
    }

    o2::dataformats::VertexBase Pvtx;
    Pvtx.setX(collision.posX());
    Pvtx.setY(collision.posY());
    Pvtx.setZ(collision.posZ());
    Pvtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(),
                collision.covXZ(), collision.covYZ(), collision.covZZ());
    // configure PVertexer

    o2::dataformats::VertexBase PVbase_recalculated;

    o2::vertexing::PVertexer vertexer;
    o2::conf::ConfigurableParam::updateFromString(
      "pvertexer.useMeanVertexConstraint=false"); // we want to refit w/o
                                                  // MeanVertex constraint
    vertexer.init();
    bool PVrefit_doable = vertexer.prepareVertexRefit(vec_TrkContributos, Pvtx);
    double chi2 = -1.;
    double refitX = -9999.;
    double refitY = -9999.;
    double refitZ = -9999.;
    double refitXX = -9999.;
    double refitYY = -9999.;
    double refitXY = -9999.;

    if (doPVrefit && PVrefit_doable) {
      auto Pvtx_refitted = vertexer.refitVertex(vec_useTrk_PVrefit, Pvtx);
      chi2 = Pvtx_refitted.getChi2();
      refitX = Pvtx_refitted.getX();
      refitY = Pvtx_refitted.getY();
      refitZ = Pvtx_refitted.getZ();
      refitXX = Pvtx_refitted.getSigmaX2();
      refitYY = Pvtx_refitted.getSigmaY2();
      refitXY = Pvtx_refitted.getSigmaXY();
    }

    rowEventInfo(relTS, refitX, refitY, refitZ, refitXX, refitYY, refitXY, chi2,
                 nContrib);

    //    LOGP(info,"chi2: {}, Ncont: {}, nonctr:
    //    {}",chi2,nContrib,nNonContrib);

    histos.fill(HIST("chisquare_Refitted"), chi2);
    if (nContrib > nContribMin && nContrib < nContribMax &&
        (chi2 / nContrib) < 4.0 && chi2 > 0) {
      histos.fill(HIST("vertexx_Refitted"), refitX);
      histos.fill(HIST("vertexy_Refitted"), refitY);

      histos.fill(HIST("vertexx_Refitted_timestamp"), relTS, refitX);
      histos.fill(HIST("vertexy_Refitted_timestamp"), relTS, refitY);
    }
    histos.fill(HIST("chisquare"), collision.chi2());
    if (collision.chi2() / collision.numContrib() > 4)
      return;
    if (collision.numContrib() > nContribMax ||
        collision.numContrib() < nContribMin)
      return;

    histos.fill(HIST("vertexx"), collision.posX());
    histos.fill(HIST("vertexy"), collision.posY());
    histos.fill(HIST("timestamp"), relTS);

    histos.fill(HIST("vertexx_timestamp"), relTS, collision.posX());
    histos.fill(HIST("vertexy_timestamp"), relTS, collision.posY());

    if (nContrib > nContribMin && nContrib < nContribMax &&
        (chi2 / nContrib) < 4.0 && chi2 > 0) {
      histos.fill(HIST("vertexx_Refitted_vertexx"), collision.posX(), refitX);
      histos.fill(HIST("vertexy_Refitted_vertexy"), collision.posY(), refitY);
    }

  } // need selections
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec w{adaptAnalysisTask<lumiTask>(cfgc, TaskName{"lumi"})};
  return w;
}
