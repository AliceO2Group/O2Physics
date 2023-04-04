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

#include "CCDB/BasicCCDBManager.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "CommonConstants/GeomConstants.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/TrackFwd.h"
#include "Math/MatrixFunctions.h"
#include "Math/SMatrix.h"

#include "Field/MagneticField.h"
#include "TGeoGlobalMagField.h"

#include "bestCollisionTable.h"

using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<Double_t, 5>;

// The Run 3 AO2D stores the tracks at the point of innermost update. For a
// track with ITS this is the innermost (or second innermost) ITS layer. For a
// track without ITS, this is the TPC inner wall or for loopers in the TPC even
// a radius beyond that. In order to use the track parameters, the tracks have
// to be propagated to the collision vertex which is done by this task. The task
// consumes the TracksIU and TracksCovIU tables and produces Tracks and
// TracksCov to which then the user analysis can subscribe.
//
// This task is not needed for Run 2 converted data.
// There are two versions of the task (see process flags), one producing also
// the covariance matrix and the other only the tracks table.

// This is an alternative version of the propagation task with special treatment
// of ambiguous tracks

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::track;

struct AmbiguousTrackPropagation {
  //  Produces<aod::BestCollisions> tracksBestCollisions;
  Produces<aod::BestCollisionsFwd> fwdtracksBestCollisions;
  Produces<aod::ReassignedTracksCore> tracksReassignedCore;
  Produces<aod::ReassignedTracksExtra> tracksReassignedExtra;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  int runNumber = -1;
  float Bz = 0;                                         // Magnetic field for MFT
  static constexpr double centerMFT[3] = {0, 0, -61.4}; // Field at center of MFT

  o2::base::Propagator::MatCorrType matCorr =
    o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  o2::parameters::GRPMagField* grpmag = nullptr;

  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> mVtxPath{"mVtxPath", "GLO/Calib/MeanVertex", "Path of the mean vertex file"};

  Configurable<bool> produceExtra{"produceExtra", false, "Produce table with refitted track parameters"};

  using ExtBCs = soa::Join<aod::BCs, aod::Timestamps, aod::MatchedBCCollisionsSparseMulti>;

  void init(o2::framework::InitContext& initContext)
  {
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
  }

  void initCCDB(ExtBCs::iterator const& bc)
  {
    if (runNumber == bc.runNumber()) {
      return;
    }
    grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, bc.timestamp());
    LOG(info) << "Setting magnetic field to current " << grpmag->getL3Current()
              << " A for run " << bc.runNumber()
              << " from its GRPMagField CCDB object";
    o2::base::Propagator::initFieldFromGRP(grpmag);
    runNumber = bc.runNumber();

    if (doprocessMFT) {
      o2::field::MagneticField* field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
      Bz = field->getBz(centerMFT);
      LOG(info) << "The field at the center of the MFT is Bz = " << Bz;
    }
  }

  static constexpr TrackSelectionFlags::flagtype trackSelectionITS =
    TrackSelectionFlags::kITSNCls | TrackSelectionFlags::kITSChi2NDF |
    TrackSelectionFlags::kITSHits;

  static constexpr TrackSelectionFlags::flagtype trackSelectionTPC =
    TrackSelectionFlags::kTPCNCls |
    TrackSelectionFlags::kTPCCrossedRowsOverNCls |
    TrackSelectionFlags::kTPCChi2NDF;

  using ExTracksSel = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection>;

  void processCentral(ExTracksSel const&,
                      aod::Collisions const&, ExtBCs const& bcs,
                      aod::AmbiguousTracks const& atracks)
  {
    if (bcs.size() == 0) {
      return;
    }
    initCCDB(bcs.begin());

    gpu::gpustd::array<float, 2> dcaInfo;
    float bestDCA[2];
    o2::track::TrackParametrization<float> bestTrackPar;
    for (auto& atrack : atracks) {
      dcaInfo[0] = 999; // DCAxy
      dcaInfo[1] = 999; // DCAz
      bestDCA[0] = 999;
      bestDCA[1] = 999;

      auto track = atrack.track_as<ExTracksSel>();
      auto bestCol = track.has_collision() ? track.collisionId() : -1;
      if ((track.trackCutFlag() & trackSelectionITS) != trackSelectionITS) {
        continue;
      }
      if ((track.detectorMap() & (uint8_t)o2::aod::track::TPC) == (uint8_t)o2::aod::track::TPC) {
        if ((track.trackCutFlag() & trackSelectionTPC) != trackSelectionTPC) {
          continue;
        }
      }
      // Only re-propagate tracks which have passed the innermost wall of the
      // TPC (e.g. skipping loopers etc).
      auto trackPar = getTrackPar(track);
      if (track.x() < o2::constants::geom::XTPCInnerRef + 0.1) {
        auto compatibleBCs = atrack.bc_as<ExtBCs>();
        for (auto& bc : compatibleBCs) {
          if (!bc.has_collisions()) {
            continue;
          }
          auto collisions = bc.collisions();
          for (auto const& collision : collisions) {
            o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackPar, 2.f, matCorr, &dcaInfo);
            if ((std::abs(dcaInfo[0]) < std::abs(bestDCA[0])) && (std::abs(dcaInfo[1]) < std::abs(bestDCA[1]))) {
              bestCol = collision.globalIndex();
              bestDCA[0] = dcaInfo[0];
              bestDCA[1] = dcaInfo[1];
              bestTrackPar = trackPar;
            }
          }
        }
      }
      tracksReassignedCore(bestCol, track.globalIndex(), bestDCA[0], bestDCA[1]);
      if (produceExtra) {
        tracksReassignedExtra(bestTrackPar.getX(), bestTrackPar.getAlpha(),
                              bestTrackPar.getY(), bestTrackPar.getZ(), bestTrackPar.getSnp(),
                              bestTrackPar.getTgl(), bestTrackPar.getQ2Pt(), bestTrackPar.getPt(),
                              bestTrackPar.getP(), bestTrackPar.getEta(), bestTrackPar.getPhi());
      }
    }
  }
  PROCESS_SWITCH(AmbiguousTrackPropagation, processCentral, "Fill ReassignedTracks for central ambiguous tracks", true);

  void processMFT(aod::MFTTracks const&,
                  aod::Collisions const&, ExtBCs const& bcs,
                  aod::AmbiguousMFTTracks const& atracks)
  {

    if (bcs.size() == 0) {
      return;
    }
    if (atracks.size() == 0) {
      return;
    }
    initCCDB(bcs.begin());

    // Only on DCAxy
    float dcaInfo;
    float bestDCA;
    o2::track::TrackParCovFwd bestTrackPar;

    for (auto& atrack : atracks) {
      dcaInfo = 999; // DCAxy
      bestDCA = 999;

      auto track = atrack.mfttrack();
      auto bestCol = track.has_collision() ? track.collisionId() : -1;

      std::vector<double> v1; // Temporary null vector for the computation of the covariance matrix
      SMatrix55 tcovs(v1.begin(), v1.end());
      SMatrix5 tpars(track.x(), track.y(), track.phi(), track.tgl(), track.signed1Pt());
      o2::track::TrackParCovFwd trackPar{track.z(), tpars, tcovs, track.chi2()};

      auto compatibleBCs = atrack.bc_as<ExtBCs>();
      for (auto& bc : compatibleBCs) {
        if (!bc.has_collisions()) {
          continue;
        }
        auto collisions = bc.collisions();
        for (auto const& collision : collisions) {

          trackPar.propagateToZhelix(collision.posZ(), Bz); // track parameters propagation to the position of the z vertex

          const auto dcaX(trackPar.getX() - collision.posX());
          const auto dcaY(trackPar.getY() - collision.posY());
          dcaInfo = std::sqrt(dcaX * dcaX + dcaY * dcaY);

          if ((dcaInfo < bestDCA)) {
            bestCol = collision.globalIndex();
            bestDCA = dcaInfo;
            bestTrackPar = trackPar;
          }
        }
      }

      fwdtracksBestCollisions(
        bestCol, bestDCA, bestTrackPar.getX(),
        bestTrackPar.getY(), bestTrackPar.getZ(),
        bestTrackPar.getTgl(), bestTrackPar.getInvQPt(), bestTrackPar.getPt(),
        bestTrackPar.getP(), bestTrackPar.getEta(), bestTrackPar.getPhi());
    }
  }
  PROCESS_SWITCH(AmbiguousTrackPropagation, processMFT, "Fill BestCollisionsFwd for MFT ambiguous tracks", false);
};

//****************************************************************************************
/**
 * Workflow definition.
 */
//****************************************************************************************
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<AmbiguousTrackPropagation>(cfgc)};
  return workflow;
}
