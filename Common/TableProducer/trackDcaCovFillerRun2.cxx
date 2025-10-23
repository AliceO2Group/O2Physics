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

/// \file trackDcaCovFillerRun2.cxx
/// \author Aimeric Landou <aimeric.landou@cern.ch>, CERN
/// \brief Fills DCA and DCA Cov tables for Run 2 tracks
// Run 2 AO2Ds cannot have their dcacov filled by the current track-propagation workflow as the workflow isn't designed for them, given Run 2 tracks are already propagated to the PV.
// This task fills the DCA Cov (and DCA) tables for Run 2 tracks by "propagating" the tracks (though given they are already at the PV it doesn't actually do the propagation) and retrieving the DCA and DCA cov given by the propagateToDCABxByBz function

#include "Common/Core/TableHelper.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <DataFormatsCalibration/MeanVertexObject.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DataFormatsParameters/GRPObject.h>
#include <DetectorsBase/MatLayerCylSet.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/ConfigContext.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/WorkflowSpec.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/DCA.h>
#include <ReconstructionDataFormats/TrackParametrization.h>
#include <ReconstructionDataFormats/TrackParametrizationWithError.h>
#include <ReconstructionDataFormats/Vertex.h>

#include <array>
#include <string>

using namespace o2;
using namespace o2::framework;
// using namespace o2::framework::expressions;

struct TrackDcaCovFillerRun2 {
  Produces<aod::TracksDCA> tracksDCA;
  Produces<aod::TracksDCACov> tracksDCACov;

  // Produces<aod::TrackTunerTable> tunertable;

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  bool fillTracksDCA = false;
  bool fillTracksDCACov = false;
  int runNumber = -1;

  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  const o2::dataformats::MeanVertexObject* mMeanVtx = nullptr;
  o2::parameters::GRPMagField* grpmag = nullptr;
  o2::base::MatLayerCylSet* lut = nullptr;

  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "CCDB path of the grp file (run2)"};
  Configurable<std::string> mVtxPath{"mVtxPath", "GLO/Calib/MeanVertex", "Path of the mean vertex file"};

  HistogramRegistry registry{"registry"};

  void init(o2::framework::InitContext& initContext)
  {
    // Checking if the tables are requested in the workflow and enabling them
    fillTracksDCA = isTableRequiredInWorkflow(initContext, "TracksDCA");
    fillTracksDCACov = isTableRequiredInWorkflow(initContext, "TracksDCACov");

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (runNumber == bc.runNumber()) {
      return;
    }

    // Run 2 GRP object
    o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(ccdbPathGrp, bc.timestamp());
    if (grpo == nullptr) {
      LOGF(fatal, "Run 2 GRP object (type o2::parameters::GRPObject) is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
    }
    o2::base::Propagator::initFieldFromGRP(grpo);
    LOGF(info, "Setting magnetic field to %d kG for run %d from its GRP CCDB object (type o2::parameters::GRPObject)", grpo->getNominalL3Field(), bc.runNumber());

    mMeanVtx = ccdb->getForTimeStamp<o2::dataformats::MeanVertexObject>(mVtxPath, bc.timestamp());
    runNumber = bc.runNumber();
  }

  // Running variables
  std::array<float, 2> mDcaInfo;
  o2::dataformats::DCA mDcaInfoCov;
  o2::dataformats::VertexBase mVtx;
  o2::track::TrackParametrization<float> mTrackPar;
  o2::track::TrackParametrizationWithError<float> mTrackParCov;

  template <typename TTrack, typename TParticle, bool isMc, bool fillCovMat = false>
  void fillTrackTables(TTrack const& tracks,
                       TParticle const&,
                       aod::Collisions const&,
                       aod::BCsWithTimestamps const& bcs)
  {
    if (bcs.size() == 0) {
      return;
    }
    initCCDB(bcs.begin());

    if constexpr (fillCovMat) {
      if (fillTracksDCACov) {
        tracksDCACov.reserve(tracks.size());
      }
    } else {
      if (fillTracksDCA) {
        tracksDCA.reserve(tracks.size());
      }
    }

    for (auto const& track : tracks) {
      if constexpr (fillCovMat) {
        if (fillTracksDCA || fillTracksDCACov) {
          mDcaInfoCov.set(999, 999, 999, 999, 999);
        }
        setTrackParCov(track, mTrackParCov);
      } else {
        if (fillTracksDCA) {
          mDcaInfo[0] = 999;
          mDcaInfo[1] = 999;
        }
        setTrackPar(track, mTrackPar);
      }

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

      if constexpr (fillCovMat) {
        if (fillTracksDCA) {
          tracksDCA(mDcaInfoCov.getY(), mDcaInfoCov.getZ());
        }
        if (fillTracksDCACov) {
          tracksDCACov(mDcaInfoCov.getSigmaY2(), mDcaInfoCov.getSigmaZ2());
        }
      } else {
        if (fillTracksDCA) {
          tracksDCA(mDcaInfo[0], mDcaInfo[1]);
        }
      }
    }
  }

  void processCovariance(soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov> const& tracks, aod::Collisions const& collisions, aod::BCsWithTimestamps const& bcs)
  {
    fillTrackTables</*TTrack*/ soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov>, /*Particle*/ soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov>, /*isMc = */ false, /*fillCovMat =*/true>(tracks, tracks, collisions, bcs);
  }
  PROCESS_SWITCH(TrackDcaCovFillerRun2, processCovariance, "Process with covariance", false);

  void processStandard(soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov> const& tracks, aod::Collisions const& collisions, aod::BCsWithTimestamps const& bcs)
  {
    fillTrackTables</*TTrack*/ soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov>, /*Particle*/ soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov>, /*isMc = */ false, /*fillCovMat =*/false>(tracks, tracks, collisions, bcs);
  }
  PROCESS_SWITCH(TrackDcaCovFillerRun2, processStandard, "Process without covariance", true);
};

//****************************************************************************************
/**
 * Workflow definition.
 */
//****************************************************************************************
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<TrackDcaCovFillerRun2>(cfgc)};
  return workflow;
}
