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

//
// Task to add a table of track parameters propagated to the primary vertex
//

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/TrackPropagation.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsCommonDataFormats/NameConf.h"
#include "CCDB/CcdbApi.h"
#include "DataFormatsParameters/GRPObject.h"
#include "CCDB/BasicCCDBManager.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "DataFormatsCalibration/MeanVertexObject.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2
{
namespace analysis
{
namespace trackpropagation
{
constexpr long run3lut_timestamp = (1665695116725 + 1634159124442) / 2;
constexpr long run3grp_timestamp = (1619781650000 + 1619781529000) / 2;

} // namespace trackpropagation
} // namespace analysis
} // namespace o2

struct TrackPropagation {

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  Produces<aod::tracksPropagated> tracksPropagated;
  Produces<aod::tracksParPropagated> tracksParPropagated;
  Produces<aod::TracksExtended> tracksExtended;

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  bool fillTracksPropagated = false;
  bool fillTracksParPropagated = false;
  bool fillTracksExtended = false;

  o2::base::Propagator::MatCorrType matCorr;

  const o2::dataformats::MeanVertexObject* mVtx;

  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/Geometry", "Path of the geometry file"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> mVtxPath{"mVtxPath", "GLO/Calib/MeanVertex", "Path of the mean vertex file"};

  Configurable<bool> fillQAHists{"fillQAHists", false, "option to fill the QA histograms"};
  Configurable<bool> isForceFillTracksPropagated{"isForceFillTracksPropagated", false, "option to fill the tracksPropagated table without workflow requirement"};
  Configurable<bool> isForceFillTracksParPropagated{"isForceFillTracksParPropagated", false, "option to fill the tracksParPropagated table without workflow requirement"};
  Configurable<bool> isForceFillTracksExtended{"isForceFillTracksExtended", false, "option to fill the tracksExtended tables without workflow requirement"};

  void init(o2::framework::InitContext& initContext)
  {

    // Checking if the tables are requested in the workflow and enabling them
    auto& workflows = initContext.services().get<RunningWorkflowInfo const>();
    for (DeviceSpec device : workflows.devices) {
      for (auto input : device.inputs) {
        const std::string tableTracksPropagated = "tracksPropagated";
        if (input.matcher.binding == tableTracksPropagated) {
          fillTracksPropagated = true;
        }
        const std::string tableTracksParPropagated = "tracksParPropagated";
        if (input.matcher.binding == tableTracksParPropagated) {
          fillTracksParPropagated = true;
        }
        const std::string tableTracksExtended = "TracksExtended";
        if (input.matcher.binding == tableTracksExtended) {
          fillTracksExtended = true;
        }
      }
    }

    if (isForceFillTracksPropagated) {
      fillTracksPropagated = true;
    }
    if (isForceFillTracksParPropagated) {
      fillTracksParPropagated = true;
    }
    if (isForceFillTracksExtended) {
      fillTracksExtended = true;
    }

    if (fillQAHists) {
      registry.add("hpt", "track #it{p}_{T} (GeV/#it{c});not propagated; propagated to PV", {HistType::kTH2F, {{200, 0., 20.}, {200, 0., 20.}}});
      registry.add("hphi", "track #phi; not propagated; propagated to PV", {HistType::kTH2F, {{140, -1., 8.}, {140, -1., 8.}}});
      registry.add("heta", "track #eta; not propagated; propagated to PV", {HistType::kTH2F, {{200, -2., 2.}, {200, -2., 2.}}});
    }

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    //auto lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->getForTimeStamp<o2::base::MatLayerCylSet>(lutPath, analysis::trackpropagation::run3lut_timestamp));
    auto lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(lutPath));

    matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>(geoPath);
      o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, analysis::trackpropagation::run3grp_timestamp);
      o2::base::Propagator::initFieldFromGRP(grpo);
      o2::base::Propagator::Instance()->setMatLUT(lut);
    }

    mVtx = ccdb->get<o2::dataformats::MeanVertexObject>(mVtxPath);
  }

  void process(aod::Tracks const& tracks, aod::Collisions const&)
  {

    gpu::gpustd::array<float, 2> dcaInfo;
    for (auto& track : tracks) {

      auto trackPar = getTrackPar(track);
      if (track.has_collision()) {
        auto const& collision = track.collision();

        o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackPar, 2.f, matCorr, &dcaInfo);
        if (fillTracksPropagated) {
          tracksPropagated(track.collisionId(), track.trackType(), trackPar.getSign(), trackPar.getPt(), trackPar.getPhi(), trackPar.getEta());
        }
        if (fillTracksParPropagated) {
          tracksParPropagated(trackPar.getX(), trackPar.getAlpha(), trackPar.getY(), trackPar.getZ(), trackPar.getSnp(), trackPar.getTgl(), trackPar.getQ2Pt());
        }
        if (fillTracksExtended) {
          tracksExtended(dcaInfo[0], dcaInfo[1]);
        }

        if (fillQAHists) {
          registry.fill(HIST("hpt"), track.pt(), trackPar.getPt());
          registry.fill(HIST("hphi"), track.phi(), trackPar.getPhi());
          registry.fill(HIST("heta"), track.eta(), trackPar.getEta());
        }
      } else {

        o2::base::Propagator::Instance()->propagateToDCABxByBz({mVtx->getX(), mVtx->getY(), mVtx->getZ()}, trackPar, 2.f, matCorr, &dcaInfo); //put mean collision information here

        if (fillTracksPropagated) {
          tracksPropagated(-1, track.trackType(), trackPar.getSign(), trackPar.getPt(), trackPar.getPhi(), trackPar.getEta());
        }
        if (fillTracksParPropagated) {
          tracksParPropagated(trackPar.getX(), trackPar.getAlpha(), trackPar.getY(), trackPar.getZ(), trackPar.getSnp(), trackPar.getTgl(), trackPar.getQ2Pt());
        }
        if (fillTracksExtended) {
          tracksExtended(dcaInfo[0], dcaInfo[1]);
        }
      }
    }
  }
};

//****************************************************************************************
/**
 * Workflow definition.
 */
//****************************************************************************************
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<TrackPropagation>(cfgc)};
  return workflow;
}
