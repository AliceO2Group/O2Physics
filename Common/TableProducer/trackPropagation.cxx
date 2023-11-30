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
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "CommonConstants/GeomConstants.h"
#include "TableHelper.h"

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

struct TrackPropagation {
  Produces<aod::StoredTracks> tracksParPropagated;
  Produces<aod::TracksExtension> tracksParExtensionPropagated;

  Produces<aod::StoredTracksCov> tracksParCovPropagated;
  Produces<aod::TracksCovExtension> tracksParCovExtensionPropagated;

  Produces<aod::TracksDCA> tracksDCA;
  Produces<aod::TracksDCACov> tracksDCACov;

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  bool fillTracksDCA = false;
  bool fillTracksDCACov = false;
  int runNumber = -1;

  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;

  const o2::dataformats::MeanVertexObject* mMeanVtx = nullptr;
  o2::parameters::GRPMagField* grpmag = nullptr;
  o2::base::MatLayerCylSet* lut = nullptr;

  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> mVtxPath{"mVtxPath", "GLO/Calib/MeanVertex", "Path of the mean vertex file"};
  Configurable<float> minPropagationRadius{"minPropagationDistance", o2::constants::geom::XTPCInnerRef + 0.1, "Only tracks which are at a smaller radius will be propagated, defaults to TPC inner wall"};

  void init(o2::framework::InitContext& initContext)
  {
    int nEnabledProcesses = 0;
    if (doprocessStandard == true) {
      LOG(info) << "Enabling processStandard";
      nEnabledProcesses++;
    }
    if (doprocessCovariance == true) {
      LOG(info) << "Enabling processCovariance";
      nEnabledProcesses++;
    }
    if (doprocessStandardWithPID == true) {
      LOG(info) << "Enabling processStandardWithPID";
      nEnabledProcesses++;
    }
    if (doprocessCovarianceWithPID == true) {
      LOG(info) << "Enabling processCovarianceWithPID";
      nEnabledProcesses++;
    }
    if (nEnabledProcesses != 1) {
      LOG(fatal) << "Exactly one process flag must be set to true. Please choose one.";
    }
    // Checking if the tables are requested in the workflow and enabling them
    fillTracksDCA = isTableRequiredInWorkflow(initContext, "TracksDCA");
    fillTracksDCACov = isTableRequiredInWorkflow(initContext, "TracksDCACov");

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(lutPath));
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
    mMeanVtx = ccdb->getForTimeStamp<o2::dataformats::MeanVertexObject>(mVtxPath, bc.timestamp());
    runNumber = bc.runNumber();
  }

  // Running variables
  gpu::gpustd::array<float, 2> mDcaInfo;
  o2::dataformats::DCA mDcaInfoCov;
  o2::dataformats::VertexBase mVtx;
  o2::track::TrackParametrization<float> mTrackPar;
  o2::track::TrackParametrizationWithError<float> mTrackParCov;

  template <typename TTrack, bool fillCovMat = false, bool useTrkPid = false>
  void fillTrackTables(TTrack const& tracks,
                       aod::Collisions const&,
                       aod::BCsWithTimestamps const& bcs)
  {
    if (bcs.size() == 0) {
      return;
    }
    initCCDB(bcs.begin());

    if constexpr (fillCovMat) {
      tracksParCovPropagated.reserve(tracks.size());
      tracksParCovExtensionPropagated.reserve(tracks.size());
      if (fillTracksDCACov) {
        tracksDCACov.reserve(tracks.size());
      }
    } else {
      tracksParPropagated.reserve(tracks.size());
      tracksParExtensionPropagated.reserve(tracks.size());
      if (fillTracksDCA) {
        tracksDCA.reserve(tracks.size());
      }
    }

    for (auto& track : tracks) {
      if constexpr (fillCovMat) {
        if (fillTracksDCACov) {
          mDcaInfoCov.set(999, 999, 999, 999, 999);
        }
        setTrackParCov(track, mTrackParCov);
        if constexpr (useTrkPid) {
          mTrackParCov.setPID(track.pidForTracking());
        }
      } else {
        if (fillTracksDCA) {
          mDcaInfo[0] = 999;
          mDcaInfo[1] = 999;
        }
        setTrackPar(track, mTrackPar);
        if constexpr (useTrkPid) {
          mTrackPar.setPID(track.pidForTracking());
        }
      }
      aod::track::TrackTypeEnum trackType = (aod::track::TrackTypeEnum)track.trackType();
      // Only propagate tracks which have passed the innermost wall of the TPC (e.g. skipping loopers etc). Others fill unpropagated.
      if (track.trackType() == aod::track::TrackIU && track.x() < minPropagationRadius) {
        if (track.has_collision()) {
          auto const& collision = track.collision();
          if constexpr (fillCovMat) {
            mVtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
            mVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
            o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, mTrackParCov, 2.f, matCorr, &mDcaInfoCov);
          } else {
            o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, mTrackPar, 2.f, matCorr, &mDcaInfo);
          }
        } else {
          if constexpr (fillCovMat) {
            mVtx.setPos({mMeanVtx->getX(), mMeanVtx->getY(), mMeanVtx->getZ()});
            mVtx.setCov(mMeanVtx->getSigmaX() * mMeanVtx->getSigmaX(), 0.0f, mMeanVtx->getSigmaY() * mMeanVtx->getSigmaY(), 0.0f, 0.0f, mMeanVtx->getSigmaZ() * mMeanVtx->getSigmaZ());
            o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, mTrackParCov, 2.f, matCorr, &mDcaInfoCov);
          } else {
            o2::base::Propagator::Instance()->propagateToDCABxByBz({mMeanVtx->getX(), mMeanVtx->getY(), mMeanVtx->getZ()}, mTrackPar, 2.f, matCorr, &mDcaInfo);
          }
        }
        trackType = aod::track::Track;
      }
      if constexpr (fillCovMat) {
        tracksParPropagated(track.collisionId(), trackType, mTrackParCov.getX(), mTrackParCov.getAlpha(), mTrackParCov.getY(), mTrackParCov.getZ(), mTrackParCov.getSnp(), mTrackParCov.getTgl(), mTrackParCov.getQ2Pt());
        tracksParExtensionPropagated(mTrackParCov.getPt(), mTrackParCov.getP(), mTrackParCov.getEta(), mTrackParCov.getPhi());
        // TODO do we keep the rho as 0? Also the sigma's are duplicated information
        tracksParCovPropagated(std::sqrt(mTrackParCov.getSigmaY2()), std::sqrt(mTrackParCov.getSigmaZ2()), std::sqrt(mTrackParCov.getSigmaSnp2()),
                               std::sqrt(mTrackParCov.getSigmaTgl2()), std::sqrt(mTrackParCov.getSigma1Pt2()), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        tracksParCovExtensionPropagated(mTrackParCov.getSigmaY2(), mTrackParCov.getSigmaZY(), mTrackParCov.getSigmaZ2(), mTrackParCov.getSigmaSnpY(),
                                        mTrackParCov.getSigmaSnpZ(), mTrackParCov.getSigmaSnp2(), mTrackParCov.getSigmaTglY(), mTrackParCov.getSigmaTglZ(), mTrackParCov.getSigmaTglSnp(),
                                        mTrackParCov.getSigmaTgl2(), mTrackParCov.getSigma1PtY(), mTrackParCov.getSigma1PtZ(), mTrackParCov.getSigma1PtSnp(), mTrackParCov.getSigma1PtTgl(),
                                        mTrackParCov.getSigma1Pt2());
        if (fillTracksDCA) {
          tracksDCA(mDcaInfoCov.getY(), mDcaInfoCov.getZ());
        }
        if (fillTracksDCACov) {
          tracksDCACov(mDcaInfoCov.getSigmaY2(), mDcaInfoCov.getSigmaZ2());
        }
      } else {
        tracksParPropagated(track.collisionId(), trackType, mTrackPar.getX(), mTrackPar.getAlpha(), mTrackPar.getY(), mTrackPar.getZ(), mTrackPar.getSnp(), mTrackPar.getTgl(), mTrackPar.getQ2Pt());
        tracksParExtensionPropagated(mTrackPar.getPt(), mTrackPar.getP(), mTrackPar.getEta(), mTrackPar.getPhi());
        if (fillTracksDCA) {
          tracksDCA(mDcaInfo[0], mDcaInfo[1]);
        }
      }
    }
  }

  void processStandard(aod::StoredTracksIU const& tracks, aod::Collisions const& collisions, aod::BCsWithTimestamps const& bcs)
  {
    fillTrackTables<aod::StoredTracksIU, /*fillCovMat =*/false, /*useTrkPid =*/false>(tracks, collisions, bcs);
  }
  PROCESS_SWITCH(TrackPropagation, processStandard, "Process without covariance", true);

  void processStandardWithPID(soa::Join<aod::StoredTracksIU, aod::TracksExtra> const& tracks, aod::Collisions const& collisions, aod::BCsWithTimestamps const& bcs)
  {
    fillTrackTables<soa::Join<aod::StoredTracksIU, aod::TracksExtra>, /*fillCovMat =*/false, /*useTrkPid =*/true>(tracks, collisions, bcs);
  }
  PROCESS_SWITCH(TrackPropagation, processStandardWithPID, "Process without covariance and with PID in tracking", false);

  void processCovariance(soa::Join<aod::StoredTracksIU, aod::TracksCovIU> const& tracks, aod::Collisions const& collisions, aod::BCsWithTimestamps const& bcs)
  {
    fillTrackTables<soa::Join<aod::StoredTracksIU, aod::TracksCovIU>, /*fillCovMat =*/true, /*useTrkPid =*/false>(tracks, collisions, bcs);
  }
  PROCESS_SWITCH(TrackPropagation, processCovariance, "Process with covariance", false);

  void processCovarianceWithPID(soa::Join<aod::StoredTracksIU, aod::TracksCovIU, aod::TracksExtra> const& tracks, aod::Collisions const& collisions, aod::BCsWithTimestamps const& bcs)
  {
    fillTrackTables<soa::Join<aod::StoredTracksIU, aod::TracksCovIU, aod::TracksExtra>, /*fillCovMat =*/true, /*useTrkPid =*/false>(tracks, collisions, bcs);
  }
  PROCESS_SWITCH(TrackPropagation, processCovarianceWithPID, "Process with covariance and with PID in tracking", false);
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
