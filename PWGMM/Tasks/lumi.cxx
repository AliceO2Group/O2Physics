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
// usage: o2-analysis-timestamp -b --aod-file AO2D.root --configuration json://./config.json | o2-analysis-trackextension -b | o2-analysis-trackselection -b --isRun3 <0, 1> | o2-analysis-mm-lumi -b --configuration json://./config.json

#include <array>
#include <cmath>
#include <iostream>
#include <vector>

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"

#include "CommonUtils/NameConf.h"

#include "DetectorsBase/Propagator.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoAHelpers.h"

#include "DetectorsVertexing/PVertexer.h"

#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/PrimaryVertex.h"
#include "ReconstructionDataFormats/Vertex.h"

#include "DetectorsBase/GeometryManager.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPObject.h"

using namespace o2::framework;
using namespace o2::framework::expressions;

struct lumiTask {
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut;
  const char* ccdbpath_lut = "GLO/Param/MatLUT";
  const char* ccdbpath_geo = "GLO/Config/Geometry";
  const char* ccdbpath_grp = "GLO/GRP/GRP";
  const char* ccdburl = "http://alice-ccdb.cern.ch";
  int mRunNumber;

  Configurable<uint64_t> ftts{"ftts", 1530319778000, "First time of time stamp"};
  Configurable<int> nContribMax{"nContribMax", 2500, "Maximum number of contributors"};
  Configurable<int> nContribMin{"nContribMin", 10, "Minimum number of contributors"};

  HistogramRegistry histos{"histos", {
                                       {"vertexx", "", {HistType::kTH1F, {{1000, -1, 1, "x"}}}},                                 //
                                       {"vertexy", "", {HistType::kTH1F, {{1000, -1, 1, "y"}}}},                                 //
                                       {"timestamp", "", {HistType::kTH1F, {{20000, 0, 2e7, "t"}}}},                             //
                                       {"vertexx_timestamp", "", {HistType::kTH2F, {{1000, 0, 4e6, "t"}, {20000, -1, 1, "x"}}}}, //
                                       {"vertexy_timestamp", "", {HistType::kTH2F, {{1000, 0, 4e6, "t"}, {20000, -1, 1, "y"}}}}, //
                                       {"chisquare", "", {HistType::kTH1F, {{1000, 0, 100, "#chi^{2}"}}}},                       //

                                       {"vertexx_Refitted", "", {HistType::kTH1F, {{1000, -1, 1, "x"}}}},                                 //
                                       {"vertexy_Refitted", "", {HistType::kTH1F, {{1000, -1, 1, "y"}}}},                                 //
                                       {"vertexx_Refitted_timestamp", "", {HistType::kTH2F, {{1000, 0, 4e6, "t"}, {20000, -1, 1, "x"}}}}, //
                                       {"vertexy_Refitted_timestamp", "", {HistType::kTH2F, {{1000, 0, 4e6, "t"}, {20000, -1, 1, "y"}}}}, //
                                       {"chisquare_Refitted", "", {HistType::kTH1F, {{1000, 0, 100, "#chi^{2}"}}}},                       //

                                       {"vertexx_Refitted_vertexx", "", {HistType::kTH2F, {{1000, -1, 1, "x"}, {1000, -1, 1, "rx"}}}}, //
                                       {"vertexy_Refitted_vertexy", "", {HistType::kTH2F, {{1000, -1, 1, "y"}, {1000, -1, 1, "ry"}}}}  //
                                     }};
  bool doPVrefit = true;

  void init(InitContext&)
  {
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbpath_lut));
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>(ccdbpath_geo);
    }
    mRunNumber = 0;
  }

  void process(aod::Collision const& collision,
               aod::BCsWithTimestamps const&,
               o2::soa::Join<o2::aod::Tracks, o2::aod::TrackSelection, o2::aod::TracksCov, o2::aod::TracksExtra, o2::aod::TracksExtended> const& tracks,
               o2::soa::Join<o2::aod::Tracks, o2::aod::TracksCov, o2::aod::TracksExtra> const& unfiltered_tracks)
  {

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();

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

    o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
    if (mRunNumber != bc.runNumber()) {
      auto grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(ccdbpath_grp, bc.timestamp());
      if (grpo != nullptr) {
        o2::base::Propagator::initFieldFromGRP(grpo);
        o2::base::Propagator::Instance()->setMatLUT(lut);
        LOGF(info, "Setting magnetic field to %d kG for run %d from its GRP CCDB object", grpo->getNominalL3Field(), bc.runNumber());
      } else {
        LOGF(fatal, "GRP object is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
      }
      mRunNumber = bc.runNumber();
    }

    o2::dataformats::VertexBase Pvtx;
    Pvtx.setX(collision.posX());
    Pvtx.setY(collision.posY());
    Pvtx.setZ(collision.posZ());
    Pvtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
    // configure PVertexer

    o2::dataformats::VertexBase PVbase_recalculated;

    o2::vertexing::PVertexer vertexer;
    o2::conf::ConfigurableParam::updateFromString("pvertexer.useMeanVertexConstraint=false"); // we want to refit w/o MeanVertex constraint
    vertexer.init();
    bool PVrefit_doable = vertexer.prepareVertexRefit(vec_TrkContributos, Pvtx);
    double chi2;
    double refitX;
    double refitY;

    if (doPVrefit && PVrefit_doable) {
      auto Pvtx_refitted = vertexer.refitVertex(vec_useTrk_PVrefit, Pvtx);
      chi2 = Pvtx_refitted.getChi2();
      refitX = Pvtx_refitted.getX();
      refitY = Pvtx_refitted.getY();
    }

    //    LOGP(info,"chi2: {}, Ncont: {}, nonctr: {}",chi2,nContrib,nNonContrib);

    histos.fill(HIST("chisquare_Refitted"), chi2);
    if (nContrib > nContribMin && nContrib < nContribMax && (chi2 / nContrib) < 4.0 && chi2 > 0) {
      histos.fill(HIST("vertexx_Refitted"), refitX);
      histos.fill(HIST("vertexy_Refitted"), refitY);

      histos.fill(HIST("vertexx_Refitted_timestamp"), bc.timestamp() - ftts, refitX);
      histos.fill(HIST("vertexy_Refitted_timestamp"), bc.timestamp() - ftts, refitY);
    }
    histos.fill(HIST("chisquare"), collision.chi2());
    if (collision.chi2() / collision.numContrib() > 4)
      return;
    if (collision.numContrib() > nContribMax || collision.numContrib() < nContribMin)
      return;

    histos.fill(HIST("vertexx"), collision.posX());
    histos.fill(HIST("vertexy"), collision.posY());
    histos.fill(HIST("timestamp"), bc.timestamp() - ftts);

    histos.fill(HIST("vertexx_timestamp"), bc.timestamp() - ftts, collision.posX());
    histos.fill(HIST("vertexy_timestamp"), bc.timestamp() - ftts, collision.posY());

    if (nContrib > nContribMin && nContrib < nContribMax && (chi2 / nContrib) < 4.0 && chi2 > 0) {
      histos.fill(HIST("vertexx_Refitted_vertexx"), collision.posX(), refitX);
      histos.fill(HIST("vertexy_Refitted_vertexy"), collision.posY(), refitY);
    }

  } // need selections
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec w{
    adaptAnalysisTask<lumiTask>(cfgc, TaskName{"lumi"})};
  return w;
}
