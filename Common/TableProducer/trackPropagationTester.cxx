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
// FIXME: THIS IS AN EXPERIMENTAL TASK, MEANT ONLY FOR EXPLORATORY PURPOSES.
// FIXME: PLEASE ONLY USE IT WITH EXTREME CARE. IF IN DOUBT, STICK WITH THE DEFAULT
// FIXME: TRACKPROPAGATION

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "CommonUtils/NameConf.h"
#include "CCDB/CcdbApi.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "CommonConstants/GeomConstants.h"
#include "trackSelectionRequest.h"

// The Run 3 AO2D stores the tracks at the point of innermost update. For a track with ITS this is the innermost (or second innermost)
// ITS layer. For a track without ITS, this is the TPC inner wall or for loopers in the TPC even a radius beyond that.
// In order to use the track parameters, the tracks have to be propagated to the collision vertex which is done by this task.
// The task consumes the TracksIU and TracksCovIU tables and produces Tracks and TracksCov to which then the user analysis can subscribe.
//
// This task is not needed for Run 2 converted data.
// There are two versions of the task (see process flags), one producing also the covariance matrix and the other only the tracks table.

using namespace o2;
using namespace o2::framework;
// using namespace o2::framework::expressions;

struct TrackPropagationTester {
  Produces<aod::StoredTracks> tracksParPropagated;
  Produces<aod::TracksExtension> tracksParExtensionPropagated;

  Produces<aod::StoredTracksCov> tracksParCovPropagated;
  Produces<aod::TracksCovExtension> tracksParCovExtensionPropagated;

  Produces<aod::TracksDCA> tracksDCA;

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  bool fillTracksDCA = false;
  int runNumber = -1;

  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;

  const o2::dataformats::MeanVertexObject* mVtx = nullptr;
  o2::parameters::GRPMagField* grpmag = nullptr;
  o2::base::MatLayerCylSet* lut = nullptr;

  // Track selection object in this scope: not necessarily a configurable
  trackSelectionRequest trackSels;
  // Configurable based on a struct
  // Configurable<trackSelectionRequest> trackSels{"trackSels", {}, "track selections"};

  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> mVtxPath{"mVtxPath", "GLO/Calib/MeanVertex", "Path of the mean vertex file"};
  Configurable<float> minPropagationRadius{"minPropagationDistance", o2::constants::geom::XTPCInnerRef + 0.1, "Only tracks which are at a smaller radius will be propagated, defaults to TPC inner wall"};

  // Configurables regarding what to propagate
  //  FIXME: This is dangerous and error prone for general purpose use. It is meant ONLY for testing.
  Configurable<bool> propagateUnassociated{"propagateUnassociated", false, "propagate tracks with no collision assoc"};
  Configurable<bool> propagateTPConly{"propagateTPConly", false, "propagate tracks with only TPC (no ITS, TRD, TOF)"};
  Configurable<int> minTPCClusters{"minTPCClusters", 70, "min number of TPC clusters to propagate"};
  Configurable<float> maxPropagStep{"maxPropagStep", 2.0, "max propag step"}; // to be checked systematically
                                                                              // use auto-detect configuration
  Configurable<bool> d_UseAutodetectMode{"d_UseAutodetectMode", false, "Autodetect requested track criteria"};

  bool hasEnding(std::string const& fullString, std::string const& ending)
  {
    if (fullString.length() >= ending.length()) {
      return (0 == fullString.compare(fullString.length() - ending.length(), ending.length(), ending));
    } else {
      return false;
    }
  }

  void init(o2::framework::InitContext& initContext)
  {
    const AxisSpec axisX{(int)4, 0.0f, +4.0f, "Track counter"};
    histos.add("hTrackCounter", "hTrackCounter", kTH1F, {axisX});

    if (doprocessCovariance == true && doprocessStandard == true) {
      LOGF(fatal, "Cannot enable processStandard and processCovariance at the same time. Please choose one.");
    }

    // Checking if the tables are requested in the workflow and enabling them
    auto& workflows = initContext.services().get<RunningWorkflowInfo const>();
    for (DeviceSpec const& device : workflows.devices) {
      for (auto const& input : device.inputs) {
        if (input.matcher.binding == "TracksDCA") {
          fillTracksDCA = true;
        }
      }
    }

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(lutPath));

    if (d_UseAutodetectMode) {
      LOGF(info, "*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*");
      LOGF(info, "    Track propagator self-configuration");
      LOGF(info, "*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*");
      trackSels.SetTightSelections(); // Only loosen from this point forward
      for (DeviceSpec const& device : workflows.devices) {
        // Loop over options to find track selection
        for (auto const& option : device.options) {
          if (hasEnding(option.name, ".requireTPC")) {
            bool lVal = option.defaultValue.get<bool>();
            LOGF(info, "Device %s, request TPC: %i", device.name, lVal);
            if (trackSels.getRequireTPC() == false)
              trackSels.setRequireTPC(lVal);
          }
          if (hasEnding(option.name, ".minTPCclusters")) {
            int lVal = option.defaultValue.get<int>();
            LOGF(info, "Device %s, min TPC clusters: %i", device.name, lVal);
            if (trackSels.getMinTPCClusters() > lVal)
              trackSels.setMinTPCClusters(lVal);
          }
          if (hasEnding(option.name, ".minTPCcrossedrows")) {
            int lVal = option.defaultValue.get<int>();
            LOGF(info, "Device %s, min TPC crossed rows: %i", device.name, lVal);
            if (trackSels.getMinTPCCrossedRows() > lVal)
              trackSels.setMinTPCCrossedRows(lVal);
          }
          if (hasEnding(option.name, ".minTPCcrossedrowsoverfindable")) {
            float lVal = option.defaultValue.get<float>();
            LOGF(info, "Device %s, min TPC crossed rows over findable: %.3f", device.name, lVal);
            if (trackSels.getMinTPCCrossedRowsOverFindable() > lVal)
              trackSels.setMinTPCCrossedRowsOverFindable(lVal);
          }
          if (hasEnding(option.name, ".requireITS")) {
            bool lVal = option.defaultValue.get<bool>();
            LOGF(info, "Device %s, request ITS: %i", device.name, lVal);
            if (trackSels.getRequireITS() == false)
              trackSels.setRequireITS(lVal);
          }
          if (hasEnding(option.name, ".minITSclusters")) {
            int lVal = option.defaultValue.get<int>();
            LOGF(info, "Device %s, minimum ITS clusters: %i", device.name, lVal);
            if (trackSels.getMinITSClusters() > lVal)
              trackSels.setMinITSClusters(lVal);
          }
          if (hasEnding(option.name, ".maxITSChi2percluster")) {
            float lVal = option.defaultValue.get<float>();
            LOGF(info, "Device %s, max ITS chi2/clu: %.3f", device.name, lVal);
            if (trackSels.getMaxITSChi2PerCluster() < lVal)
              trackSels.setMaxITSChi2PerCluster(lVal);
          }
        }
      }
      LOGF(info, "-+*> Automatic self-config ended. Final settings:");
      trackSels.PrintSelections();
    }
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (runNumber == bc.runNumber()) {
      return;
    }
    grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, bc.timestamp());
    LOG(info) << "Setting magnetic field to current " << grpmag->getL3Current() << " A for run " << bc.runNumber() << " from its GRPMagField CCDB object";
    o2::base::Propagator::initFieldFromGRP(grpmag);
    o2::base::Propagator::Instance()->setMatLUT(lut);
    if (propagateUnassociated)
      mVtx = ccdb->getForTimeStamp<o2::dataformats::MeanVertexObject>(mVtxPath, bc.timestamp());
    runNumber = bc.runNumber();
  }

  template <typename TTrack, typename TTrackPar>
  void FillTracksPar(TTrack& track, aod::track::TrackTypeEnum trackType, TTrackPar& trackPar)
  {
    tracksParPropagated(track.collisionId(), trackType, trackPar.getX(), trackPar.getAlpha(), trackPar.getY(), trackPar.getZ(), trackPar.getSnp(), trackPar.getTgl(), trackPar.getQ2Pt());
    tracksParExtensionPropagated(trackPar.getPt(), trackPar.getP(), trackPar.getEta(), trackPar.getPhi());
  }

  void processStandard(soa::Join<aod::TracksIU, aod::TracksExtra> const& tracks, aod::Collisions const&, aod::BCsWithTimestamps const& bcs)
  {
    if (bcs.size() == 0) {
      return;
    }
    initCCDB(bcs.begin());

    std::array<float, 2> dcaInfo;

    int lNAll = 0;
    int lNaccTPC = 0;
    int lNaccNotTPCOnly = 0;
    int lNPropagated = 0;
    bool passTPCclu = kFALSE;
    bool passNotTPCOnly = kFALSE;

    for (auto& track : tracks) {
      // Selection criteria
      passTPCclu = kFALSE;
      passNotTPCOnly = kFALSE;
      lNAll++;
      if (track.tpcNClsFound() >= minTPCClusters) {
        passTPCclu = kTRUE;
        lNaccTPC++;
      }
      if ((track.hasTPC() && !track.hasITS() && !track.hasTRD() && !track.hasTOF()) || propagateTPConly) {
        passNotTPCOnly = kTRUE;
        lNaccNotTPCOnly++;
      }

      dcaInfo[0] = 999;
      dcaInfo[1] = 999;
      aod::track::TrackTypeEnum trackType = (aod::track::TrackTypeEnum)track.trackType();
      auto trackPar = getTrackPar(track);
      // Only propagate tracks which have passed the innermost wall of the TPC (e.g. skipping loopers etc). Others fill unpropagated.
      if (track.trackType() == aod::track::TrackIU && track.x() < minPropagationRadius && passTPCclu && passNotTPCOnly) {
        if (track.has_collision()) {
          auto const& collision = track.collision();
          o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackPar, maxPropagStep, matCorr, &dcaInfo);
          trackType = aod::track::Track;
          lNPropagated++;
        } else {
          if (propagateUnassociated) {
            o2::base::Propagator::Instance()->propagateToDCABxByBz({mVtx->getX(), mVtx->getY(), mVtx->getZ()}, trackPar, maxPropagStep, matCorr, &dcaInfo);
            trackType = aod::track::Track;
            lNPropagated++;
          }
        }
      }
      FillTracksPar(track, trackType, trackPar);
      if (fillTracksDCA) {
        tracksDCA(dcaInfo[0], dcaInfo[1]);
      }
    }
    // Fill only per table (not per track). ROOT FindBin is slow
    histos.fill(HIST("hTrackCounter"), 0.5, lNAll);
    histos.fill(HIST("hTrackCounter"), 1.5, lNaccTPC);
    histos.fill(HIST("hTrackCounter"), 2.5, lNaccNotTPCOnly);
    histos.fill(HIST("hTrackCounter"), 3.5, lNPropagated);
  }
  PROCESS_SWITCH(TrackPropagationTester, processStandard, "Process without covariance", true);

  void processCovariance(soa::Join<aod::TracksIU, aod::TracksCovIU, aod::TracksExtra> const& tracks, aod::Collisions const&, aod::BCsWithTimestamps const& bcs)
  {
    if (bcs.size() == 0) {
      return;
    }
    initCCDB(bcs.begin());

    o2::dataformats::DCA dcaInfoCov;
    o2::dataformats::VertexBase vtx;

    int lNAll = 0;
    int lNaccTPC = 0;
    int lNaccNotTPCOnly = 0;
    int lNPropagated = 0;
    bool passTPCclu = kFALSE;
    bool passNotTPCOnly = kFALSE;

    for (auto& track : tracks) {
      // Selection criteria
      passTPCclu = kFALSE;
      passNotTPCOnly = kFALSE;
      lNAll++;
      if (track.tpcNClsFound() >= minTPCClusters) {
        passTPCclu = kTRUE;
        lNaccTPC++;
      }
      if ((track.hasTPC() && !track.hasITS() && !track.hasTRD() && !track.hasTOF()) || propagateTPConly) {
        passNotTPCOnly = kTRUE;
        lNaccNotTPCOnly++;
      }

      dcaInfoCov.set(999, 999, 999, 999, 999);
      auto trackParCov = getTrackParCov(track);
      aod::track::TrackTypeEnum trackType = (aod::track::TrackTypeEnum)track.trackType();
      // Only propagate tracks which have passed the innermost wall of the TPC (e.g. skipping loopers etc). Others fill unpropagated.
      if (track.trackType() == aod::track::TrackIU && track.x() < minPropagationRadius && passTPCclu && passNotTPCOnly) {
        if (track.has_collision()) {
          auto const& collision = track.collision();
          vtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
          vtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
          o2::base::Propagator::Instance()->propagateToDCABxByBz(vtx, trackParCov, maxPropagStep, matCorr, &dcaInfoCov);
          trackType = aod::track::Track;
          lNPropagated++;
        } else {
          if (propagateUnassociated) {
            vtx.setPos({mVtx->getX(), mVtx->getY(), mVtx->getZ()});
            vtx.setCov(mVtx->getSigmaX() * mVtx->getSigmaX(), 0.0f, mVtx->getSigmaY() * mVtx->getSigmaY(), 0.0f, 0.0f, mVtx->getSigmaZ() * mVtx->getSigmaZ());
            o2::base::Propagator::Instance()->propagateToDCABxByBz(vtx, trackParCov, maxPropagStep, matCorr, &dcaInfoCov);
            trackType = aod::track::Track;
            lNPropagated++;
          }
        }
      }
      FillTracksPar(track, trackType, trackParCov);
      if (fillTracksDCA) {
        tracksDCA(dcaInfoCov.getY(), dcaInfoCov.getZ());
      }
      // TODO do we keep the rho as 0? Also the sigma's are duplicated information
      tracksParCovPropagated(std::sqrt(trackParCov.getSigmaY2()), std::sqrt(trackParCov.getSigmaZ2()), std::sqrt(trackParCov.getSigmaSnp2()),
                             std::sqrt(trackParCov.getSigmaTgl2()), std::sqrt(trackParCov.getSigma1Pt2()), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
      tracksParCovExtensionPropagated(trackParCov.getSigmaY2(), trackParCov.getSigmaZY(), trackParCov.getSigmaZ2(), trackParCov.getSigmaSnpY(),
                                      trackParCov.getSigmaSnpZ(), trackParCov.getSigmaSnp2(), trackParCov.getSigmaTglY(), trackParCov.getSigmaTglZ(), trackParCov.getSigmaTglSnp(),
                                      trackParCov.getSigmaTgl2(), trackParCov.getSigma1PtY(), trackParCov.getSigma1PtZ(), trackParCov.getSigma1PtSnp(), trackParCov.getSigma1PtTgl(),
                                      trackParCov.getSigma1Pt2());
    }
    // Fill only per table (not per track). ROOT FindBin is slow
    histos.fill(HIST("hTrackCounter"), 0.5, lNAll);
    histos.fill(HIST("hTrackCounter"), 1.5, lNaccTPC);
    histos.fill(HIST("hTrackCounter"), 2.5, lNaccNotTPCOnly);
    histos.fill(HIST("hTrackCounter"), 3.5, lNPropagated);
  }
  PROCESS_SWITCH(TrackPropagationTester, processCovariance, "Process with covariance", false);
};

//****************************************************************************************
/**
 * Workflow definition.
 */
//****************************************************************************************
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<TrackPropagationTester>(cfgc)};
  return workflow;
}
