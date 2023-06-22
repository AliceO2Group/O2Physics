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

// \file   trackPropagation.cxx
// \author Anton Alkin <anton.alkin@cern.ch>
// \author Sarah Herrmann <sarah.herrmann@cern.ch>
//
// \brief This code loops over central and MFT tracks and among the compatible
// collisions to this track, picks the one with the smallest DCAxy and puts it
// in a table

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

#include "Common/DataModel/CollisionAssociationTables.h"
#include "bestCollisionTable.h"

using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<Double_t, 5>;

// This is a special version of the propagation task chosing the closest vertex
// among the compatible, which is defined by track-to-collision-associator

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::track;

AxisSpec DCAxyAxis = {500, -1, 50};

struct AmbiguousTrackPropagation {
  //  Produces<aod::BestCollisions> tracksBestCollisions;
  Produces<aod::BestCollisionsFwd> fwdtracksBestCollisions;
  Produces<aod::BestCollFwdExtra> fwdtracksBestCollExtra;
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
  Configurable<bool> produceHistos{"produceHistos", false, "Produce control histograms"};

  HistogramRegistry registry{
    "registry",
    {

    } //
  };

  using ExtBCs = soa::Join<aod::BCs, aod::Timestamps, aod::MatchedBCCollisionsSparseMulti>;

  void init(o2::framework::InitContext& initContext)
  {
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    if (produceHistos) {
      registry.add({"DeltaZ", " ; #Delta#it{z}", {HistType::kTH1F, {{201, -10.1, 10.1}}}});
      registry.add({"TracksDCAXY", " ; DCA_{XY} (cm)", {HistType::kTH1F, {DCAxyAxis}}});
      registry.add({"ReassignedDCAXY", " ; DCA_{XY} (cm)", {HistType::kTH1F, {DCAxyAxis}}});
      registry.add({"TracksOrigDCAXY", " ; DCA_{XY} (wrt orig coll) (cm)", {HistType::kTH1F, {DCAxyAxis}}});
      registry.add({"TracksAmbDegree", " ; N_{coll}^{comp}", {HistType::kTH1I, {{41, -0.5, 40.5}}}});
      registry.add({"TrackIsAmb", " ; isAmbiguous", {HistType::kTH1I, {{2, -0.5, 1.5}}}});
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

    if (doprocessMFT || doprocessMFTReassoc) {
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

  using ExTracksSel = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::TrackCompColls>;

  void processCentral(ExTracksSel const& tracks,
                      aod::Collisions const& collisions,
                      ExtBCs const&)
  {
    auto bc = collisions.begin().bc_as<ExtBCs>();
    initCCDB(bc);

    gpu::gpustd::array<float, 2> dcaInfo;
    float bestDCA[2];
    o2::track::TrackParametrization<float> bestTrackPar;
    for (auto& track : tracks) {
      dcaInfo[0] = track.dcaXY(); // DCAxy
      dcaInfo[1] = track.dcaZ();  // DCAz
      bestDCA[0] = dcaInfo[0];
      bestDCA[1] = dcaInfo[1];

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
        auto ids = track.compatibleCollIds();
        if (ids.empty() || (ids.size() == 1 && bestCol == ids[0])) {
          continue;
        }
        auto compatibleColls = track.compatibleColl();
        for (auto& collision : compatibleColls) {
          o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackPar, 2.f, matCorr, &dcaInfo);
          if ((std::abs(dcaInfo[0]) < std::abs(bestDCA[0])) && (std::abs(dcaInfo[1]) < std::abs(bestDCA[1]))) {
            bestCol = collision.globalIndex();
            bestDCA[0] = dcaInfo[0];
            bestDCA[1] = dcaInfo[1];
            bestTrackPar = trackPar;
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

    // Minimum only on DCAxy
    float dcaInfo;
    float bestDCA, bestDCAx, bestDCAy;
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

      int degree = 0; // degree of ambiguity of the track

      auto compatibleBCs = atrack.bc_as<ExtBCs>();
      for (auto& bc : compatibleBCs) {
        if (!bc.has_collisions()) {
          continue;
        }
        auto collisions = bc.collisions();
        for (auto const& collision : collisions) {
          degree++;
          trackPar.propagateToZhelix(collision.posZ(), Bz); // track parameters propagation to the position of the z vertex

          const auto dcaX(trackPar.getX() - collision.posX());
          const auto dcaY(trackPar.getY() - collision.posY());
          dcaInfo = std::sqrt(dcaX * dcaX + dcaY * dcaY);

          if ((dcaInfo < bestDCA)) {
            bestCol = collision.globalIndex();
            bestDCA = dcaInfo;
            bestDCAx = dcaX;
            bestDCAy = dcaY;
            bestTrackPar = trackPar;
          }

          if (produceHistos) {
            registry.fill(HIST("TracksDCAXY"), dcaInfo);
          }
          if ((track.collisionId() != collision.globalIndex()) && produceHistos) {
            registry.fill(HIST("DeltaZ"), track.collision().posZ() - collision.posZ()); // deltaZ between the 1st coll zvtx and the other compatible ones
          }

          if ((collision.globalIndex() == track.collisionId()) && produceHistos) {
            registry.fill(HIST("TracksOrigDCAXY"), dcaInfo);
          }
        }
      }

      if ((bestCol != track.collisionId()) && produceHistos) {
        // reassigned
        registry.fill(HIST("ReassignedDCAXY"), bestDCA);
      }
      if (produceHistos) {
        registry.fill(HIST("TracksAmbDegree"), degree);
      }

      fwdtracksBestCollisions(-1, degree, bestCol, bestDCA, bestDCAx, bestDCAy);
      if (produceExtra) {
        fwdtracksBestCollExtra(bestTrackPar.getX(),
                               bestTrackPar.getY(), bestTrackPar.getZ(),
                               bestTrackPar.getTgl(), bestTrackPar.getInvQPt(), bestTrackPar.getPt(),
                               bestTrackPar.getP(), bestTrackPar.getEta(), bestTrackPar.getPhi());
      }
    }
  }
  PROCESS_SWITCH(AmbiguousTrackPropagation, processMFT, "Fill BestCollisionsFwd for MFT ambiguous tracks", false);

  using MFTTracksWColls = soa::Join<o2::aod::MFTTracks, aod::MFTTrkCompColls>;

  void processMFTReassoc(MFTTracksWColls const& tracks,
                         aod::Collisions const&, ExtBCs const& bcs)
  {

    if (bcs.size() == 0) {
      return;
    }
    if (tracks.size() == 0) {
      return;
    }
    initCCDB(bcs.begin());

    float dcaInfo;
    float bestDCA, bestDCAx, bestDCAy;
    o2::track::TrackParCovFwd bestTrackPar;

    for (auto& track : tracks) {
      dcaInfo = 999; // DCAxy
      bestDCA = 999;

      auto bestCol = track.has_collision() ? track.collisionId() : -1;

      // auto ids = track.compatibleCollIds();

      // if (ids.empty() || (ids.size() == 1 && bestCol == ids[0]))
      // {
      //   continue;
      // }

      auto compatibleColls = track.compatibleColl();

      std::vector<double> v1; // Temporary null vector for the computation of the covariance matrix
      SMatrix55 tcovs(v1.begin(), v1.end());
      SMatrix5 tpars(track.x(), track.y(), track.phi(), track.tgl(), track.signed1Pt());
      o2::track::TrackParCovFwd trackPar{track.z(), tpars, tcovs, track.chi2()};

      for (auto& collision : compatibleColls) {

        trackPar.propagateToZhelix(collision.posZ(), Bz); // track parameters propagation to the position of the z vertex

        const auto dcaX(trackPar.getX() - collision.posX());
        const auto dcaY(trackPar.getY() - collision.posY());
        dcaInfo = std::sqrt(dcaX * dcaX + dcaY * dcaY);

        if ((dcaInfo < bestDCA)) {
          bestCol = collision.globalIndex();
          bestDCA = dcaInfo;
          bestDCAx = dcaX;
          bestDCAy = dcaY;
          bestTrackPar = trackPar;
        }
        if ((track.collisionId() != collision.globalIndex()) && produceHistos) {
          registry.fill(HIST("DeltaZ"), track.collision().posZ() - collision.posZ()); // deltaZ between the 1st coll zvtx and the other compatible ones
        }
        if (produceHistos) {
          registry.fill(HIST("TracksDCAXY"), dcaInfo);
        }

        if ((collision.globalIndex() == track.collisionId()) && produceHistos) {
          registry.fill(HIST("TracksOrigDCAXY"), dcaInfo);
        }
      }
      if ((bestCol != track.collisionId()) && produceHistos) {
        // reassigned
        registry.fill(HIST("ReassignedDCAXY"), bestDCA);
      }
      if (produceHistos) {
        int isAmbiguous = 0;
        registry.fill(HIST("TracksAmbDegree"), compatibleColls.size());
        if (compatibleColls.size() > 1) {
          isAmbiguous = 1;
        }
        registry.fill(HIST("TrackIsAmb"), isAmbiguous);
      }

      fwdtracksBestCollisions(track.globalIndex(), compatibleColls.size(), bestCol, bestDCA, bestDCAx, bestDCAy);
      if (produceExtra) {
        fwdtracksBestCollExtra(bestTrackPar.getX(),
                               bestTrackPar.getY(), bestTrackPar.getZ(),
                               bestTrackPar.getTgl(), bestTrackPar.getInvQPt(), bestTrackPar.getPt(),
                               bestTrackPar.getP(), bestTrackPar.getEta(), bestTrackPar.getPhi());
      }
    }
  }
  PROCESS_SWITCH(AmbiguousTrackPropagation, processMFTReassoc, "Fill BestCollisionsFwd for MFT ambiguous tracks with the new data model", false);
};

//****************************************************************************************
/**
 * Workflow definition.
 */
//****************************************************************************************
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return {adaptAnalysisTask<AmbiguousTrackPropagation>(cfgc)};
}
