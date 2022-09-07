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
#include "CCDB/CcdbApi.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "CommonConstants/GeomConstants.h"
#include "CommonUtils/NameConf.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/TrackFwd.h"
#include "Math/MatrixFunctions.h"
#include "Math/SMatrix.h"

#include "DetectorsBase/Propagator.h"
#include "Field/MagneticField.h"
#include "TGeoGlobalMagField.h"

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
// using namespace o2::framework::expressions;

namespace o2::aod
{
namespace track
{
DECLARE_SOA_INDEX_COLUMN_FULL(BestCollision, bestCollision, int32_t, Collisions, "");
DECLARE_SOA_COLUMN(BestDCAXY, bestDCAXY, float);
DECLARE_SOA_COLUMN(BestDCAZ, bestDCAZ, float);
DECLARE_SOA_COLUMN(PtStatic, pts, float);
DECLARE_SOA_COLUMN(PStatic, ps, float);
DECLARE_SOA_COLUMN(EtaStatic, etas, float);
DECLARE_SOA_COLUMN(PhiStatic, phis, float);
} // namespace track
namespace fwdtrack
{
DECLARE_SOA_INDEX_COLUMN_FULL(BestCollision, bestCollision, int32_t, Collisions, "");
DECLARE_SOA_COLUMN(BestDCAXY, bestDCAXY, float);
DECLARE_SOA_COLUMN(PtStatic, pts, float);
DECLARE_SOA_COLUMN(PStatic, ps, float);
DECLARE_SOA_COLUMN(EtaStatic, etas, float);
DECLARE_SOA_COLUMN(PhiStatic, phis, float);
} // namespace fwdtrack
DECLARE_SOA_TABLE(BestCollisions, "AOD", "BESTCOLL",
                  aod::track::BestCollisionId, aod::track::BestDCAXY,
                  aod::track::BestDCAZ, track::X, track::Alpha, track::Y,
                  track::Z, track::Snp, track::Tgl, track::Signed1Pt,
                  track::PtStatic, track::PStatic, track::EtaStatic,
                  track::PhiStatic);

DECLARE_SOA_TABLE(BestCollisionsFwd, "AOD", "BESTCOLLFWD",
                  aod::fwdtrack::BestCollisionId, aod::fwdtrack::BestDCAXY,
                  fwdtrack::X, fwdtrack::Y,
                  fwdtrack::Z, fwdtrack::Tgl, fwdtrack::Signed1Pt,
                  fwdtrack::PtStatic, fwdtrack::PStatic, fwdtrack::EtaStatic,
                  fwdtrack::PhiStatic); // Snp does not exist
namespace indices
{
DECLARE_SOA_ARRAY_INDEX_COLUMN(Collision, collisions);
}

DECLARE_SOA_TABLE(MatchedMulti, "AOD", "MAMU",
                  indices::BCId, indices::CollisionIds);
} // namespace o2::aod

struct MultiCollisionAssociation {
  Produces<aod::MatchedMulti> mm;
  struct {
    std::vector<int> colids;
  } filler;
  void process(aod::BCs::iterator const& bc, soa::SmallGroups<aod::Collisions> const& collisions)
  {
    filler.colids.clear();
    for (auto const& collision : collisions) {
      filler.colids.emplace_back(collision.globalIndex());
    }
    mm(bc.globalIndex(), filler.colids);
  }
};

struct AmbiguousTrackPropagation {
  Produces<aod::BestCollisions> tracksBestCollisions;
  Produces<aod::BestCollisionsFwd> fwdtracksBestCollisions;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  int runNumber = -1;
  float Bz = 0; // Magnetic field for MFT

  o2::base::Propagator::MatCorrType matCorr =
    o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  o2::parameters::GRPMagField* grpmag = nullptr;

  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> mVtxPath{"mVtxPath", "GLO/Calib/MeanVertex", "Path of the mean vertex file"};

  using ExtBCs = soa::Join<aod::BCs, aod::Timestamps, aod::MatchedMulti>;

  void init(o2::framework::InitContext& initContext)
  {
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>(geoPath);
    }
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

    o2::field::MagneticField* field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
    double centerMFT[3] = {0, 0, -61.4}; // Field at center of MFT
    Bz = field->getBz(centerMFT);
    LOG(info) << "The field at the center of the MFT is Bz = " << Bz;
  }

  void processCentral(soa::Join<aod::Tracks, aod::TracksExtra> const&,
                      aod::Collisions const&, ExtBCs const& bcs,
                      aod::AmbiguousTracks const& atracks)
  {
    if (bcs.size() == 0) {
      return;
    }
    initCCDB(bcs.begin());

    gpu::gpustd::array<float, 2> dcaInfo;
    float bestDCA[2];

    for (auto& atrack : atracks) {
      dcaInfo[0] = 999; // DCAxy
      dcaInfo[1] = 999; // DCAz
      bestDCA[0] = 999;
      bestDCA[1] = 999;

      auto track = atrack.track_as<soa::Join<aod::Tracks, aod::TracksExtra>>();
      auto bestCol = track.has_collision() ? track.collisionId() : -1;

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
            if ((dcaInfo[0] < bestDCA[0]) && (dcaInfo[1] < bestDCA[1])) {
              bestCol = collision.globalIndex();
              bestDCA[0] = dcaInfo[0];
              bestDCA[1] = dcaInfo[1];
            }
          }
        }
      }
      tracksBestCollisions(
        bestCol, dcaInfo[0], dcaInfo[1], trackPar.getX(), trackPar.getAlpha(),
        trackPar.getY(), trackPar.getZ(), trackPar.getSnp(),
        trackPar.getTgl(), trackPar.getQ2Pt(), trackPar.getPt(),
        trackPar.getP(), trackPar.getEta(), trackPar.getPhi());
    }
  }
  PROCESS_SWITCH(AmbiguousTrackPropagation, processCentral, "Fill BestCollisions for central ambiguous tracks", true);

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

    // Only on DCAxy
    float dcaInfo;
    float bestDCA;

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
          auto dcaInfo = std::sqrt(dcaX * dcaX + dcaY * dcaY);
          if ((dcaInfo < bestDCA)) {
            bestCol = collision.globalIndex();
            bestDCA = dcaInfo;
          }
        }
      }

      fwdtracksBestCollisions(
        bestCol, dcaInfo, trackPar.getX(),
        trackPar.getY(), trackPar.getZ(),
        trackPar.getTgl(), trackPar.getInvQPt(), trackPar.getPt(),
        trackPar.getP(), trackPar.getEta(), trackPar.getPhi());
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
  WorkflowSpec workflow{adaptAnalysisTask<AmbiguousTrackPropagation>(cfgc),
                        adaptAnalysisTask<MultiCollisionAssociation>(cfgc)};
  return workflow;
}
