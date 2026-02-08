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
///
/// \file   ambiguousTrackPropagation.cxx
/// \brief This code loops over central and MFT tracks and among the compatible
/// collisions to this track, picks the one with the smallest DCAxy/DCAz and puts it
/// in a table
/// \author Anton Alkin <anton.alkin@cern.ch>
/// \author Sarah Herrmann <sarah.herrmann@cern.ch>
/// \author Gyula Bencedi <gyula.bencedi@cern.ch>
/// \author Tulika Tripathy <tulika.tripathy@cern.ch>

#include "bestCollisionTable.h"

#include "Common/Core/fwdtrackUtilities.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/GeomConstants.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/Propagator.h"
#include "Field/MagneticField.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/TrackFwd.h"

#include "Math/MatrixFunctions.h"
#include "Math/SMatrix.h"
#include "TGeoGlobalMagField.h"

#include <string>
#include <vector>

// This is a special version of the propagation task chosing the closest vertex
// among the compatible, which is defined by track-to-collision-associator

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::track;

struct AmbiguousTrackPropagation {
  Produces<aod::BestCollisionsFwd> fwdtracksBestCollisions;
  Produces<aod::BestCollFwdExtra> fwdtracksBestCollExtra;
  Produces<aod::BestCollisionsFwd3d> fwdtracksBestCollisions3d;
  Produces<aod::BestCollisionsFwd3dExtra> fwdtracksBestCollisions3dExtra;
  Produces<aod::ReassignedTracksCore> tracksReassignedCore;
  Produces<aod::ReassignedTracksExtra> tracksReassignedExtra;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  int runNumber = -1;
  float bZ = 0;                                          // Magnetic field for MFT
  static constexpr double CcenterMFT[3] = {0, 0, -61.4}; // Field at center of MFT
  float mZShift = 0;                                     // z-vertex shift

  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  o2::parameters::GRPMagField* grpmag = nullptr;

  Configurable<std::string> ccdburl{"ccdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> mVtxPath{"mVtxPath", "GLO/Calib/MeanVertex", "Path of the mean vertex file"};

  Configurable<bool> produceExtra{"produceExtra", false, "Produce table with refitted track parameters"};
  Configurable<bool> produceHistos{"produceHistos", false, "Produce control histograms"};
  Configurable<bool> removeTrivialAssoc{"removeTrivialAssoc", false, "Skip trivial associations"};
  Configurable<uint> cfgDCAtype{"cfgDCAtype", 2, "DCA coordinate type [0: DCA-X, 1: DCA-Y, 2: DCA-XY]"};

  Configurable<bool> cfgApplyZShiftFromCCDB{"cfgApplyZShiftFromCCDB", false, "flag to apply z shift from CCDB"};
  Configurable<std::string> cfgZShiftPath{"cfgZShiftPath", "Users/m/mcoquet/ZShift", "CCDB path for z shift to apply to forward tracks"};
  Configurable<float> cfgManualZShift{"cfgManualZShift", 0.0f, "manual z-shift for propagation of global muon to PV"};

  ConfigurableAxis binsDCAxy{"binsDCAxy", {200, -1., 1.}, ""};
  ConfigurableAxis binsDCAz{"binsDCAz", {200, -1., 1.}, ""};

  HistogramRegistry registry{
    "registry",
    {

    } //
  };

  using ExtBCs = soa::Join<aod::BCs, aod::Timestamps, aod::MatchedBCCollisionsSparseMulti>;

  void init(o2::framework::InitContext& /*initContext*/)
  {

    AxisSpec dcaXYAxis = {binsDCAxy, "dcaXYAxis", "dcaXYAxis"};
    AxisSpec dcaZAxis = {binsDCAz, "dcaZAxis", "dcaZAxis"};

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    if (produceHistos) {
      if (doprocessMFT || doprocessMFTReassoc || doprocessMFTReassoc3D) {
        registry.add({"DeltaX", " ; #Delta#it{x}", {HistType::kTH1F, {{201, -10.1, 10.1}}}});
        registry.add({"DeltaY", " ; #Delta#it{y}", {HistType::kTH1F, {{201, -10.1, 10.1}}}});
        registry.add({"DeltaZ", " ; #Delta#it{z}", {HistType::kTH1F, {{201, -10.1, 10.1}}}});
        registry.add({"TracksDCAXY", " ; DCA_{XY} (cm)", {HistType::kTH1F, {dcaXYAxis}}});
        registry.add({"ReassignedDCAXY", " ; DCA_{XY} (cm)", {HistType::kTH1F, {dcaXYAxis}}});
        registry.add({"TracksOrigDCAXY", " ; DCA_{XY} (wrt orig coll) (cm)", {HistType::kTH1F, {dcaXYAxis}}});
        registry.add({"TracksAmbDegree", " ; N_{coll}^{comp}", {HistType::kTH1D, {{41, -0.5, 40.5}}}});
        registry.add({"TrackIsAmb", " ; isAmbiguous", {HistType::kTH1D, {{2, -0.5, 1.5}}}});
        if (doprocessMFTReassoc3D) {
          registry.add({"TracksAmbDegreeWoTrivial", " ; N_{coll}^{comp}", {HistType::kTH1F, {{41, -0.5, 40.5}}}});
          registry.add({"TracksFirstDCAXY", " ; DCA_{XY} (cm)", {HistType::kTH1F, {dcaXYAxis}}});
          registry.add({"TracksFirstDCAZ", " ; DCA_{Z} (cm)", {HistType::kTH1F, {dcaZAxis}}});
          registry.add({"TracksDCAZ", " ; DCA_{Z} (cm)", {HistType::kTH1F, {dcaZAxis}}});
          registry.add({"ReassignedDCAZ", " ; DCA_{Z} (cm)", {HistType::kTH1F, {dcaZAxis}}});
          registry.add({"AssignedDCAXY", " ; DCA_{XY} (cm)", {HistType::kTH1F, {dcaXYAxis}}});
          registry.add({"AssignedDCAZ", " ; DCA_{Z} (cm)", {HistType::kTH1F, {dcaZAxis}}});
          registry.add({"TracksOrigDCAZ", " ; DCA_{Z} (wrt orig coll) (cm)", {HistType::kTH1F, {dcaZAxis}}});
        }
      }
      if (doprocessCentral) {
        registry.add({"PropagationFailures", "", {HistType::kTH1F, {{5, 0.5, 5.5}}}});
        auto h = registry.get<TH1>(HIST("PropagationFailures"));
        auto* x = h->GetXaxis();
        x->SetBinLabel(1, "Total");
        x->SetBinLabel(2, "Propagated");
        x->SetBinLabel(3, "Failed 1");
        x->SetBinLabel(4, "Failed 2");
        x->SetBinLabel(5, "Failed 3+");
      }
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

    if (doprocessMFT || doprocessMFTReassoc || doprocessMFTReassoc3D) {
      o2::field::MagneticField* field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
      bZ = field->getBz(CcenterMFT);
      LOG(info) << "The field at the center of the MFT is bZ = " << bZ;
    }

    if (cfgApplyZShiftFromCCDB) {
      auto* zShift = ccdb->getForTimeStamp<std::vector<float>>(cfgZShiftPath, bc.timestamp());
      if (zShift != nullptr && !zShift->empty()) {
        LOGF(info, "reading z shift %f from %s", (*zShift)[0], cfgZShiftPath.value);
        mZShift = (*zShift)[0];
      } else {
        LOGF(info, "z shift is not found in ccdb path %s. set to 0 cm", cfgZShiftPath.value);
        mZShift = 0;
      }
    } else {
      LOGF(info, "z shift is manually set to %f cm", cfgManualZShift.value);
      mZShift = cfgManualZShift;
    }
  }

  static constexpr TrackSelectionFlags::flagtype CtrackSelectionITS =
    TrackSelectionFlags::kITSNCls | TrackSelectionFlags::kITSChi2NDF |
    TrackSelectionFlags::kITSHits;

  static constexpr TrackSelectionFlags::flagtype CtrackSelectionTPC =
    TrackSelectionFlags::kTPCNCls |
    TrackSelectionFlags::kTPCCrossedRowsOverNCls |
    TrackSelectionFlags::kTPCChi2NDF;

  using ExTracksSel = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::TrackCompColls>;

  void processCentral(ExTracksSel const& tracks,
                      aod::Collisions const&,
                      ExtBCs const& bcs)
  {
    if (bcs.size() == 0) {
      return;
    }
    auto bc = bcs.begin();
    initCCDB(bc);

    std::array<float, 2> dcaInfo;
    float bestDCA[2];
    o2::track::TrackParametrization<float> bestTrackPar;
    for (auto const& track : tracks) {
      dcaInfo[0] = track.dcaXY(); // DCAxy
      dcaInfo[1] = track.dcaZ();  // DCAz
      bestDCA[0] = dcaInfo[0];
      bestDCA[1] = dcaInfo[1];

      auto bestCol = track.has_collision() ? track.collisionId() : -1;
      if ((track.trackCutFlag() & CtrackSelectionITS) != CtrackSelectionITS) {
        continue;
      }
      if ((track.detectorMap() & (uint8_t)o2::aod::track::TPC) == (uint8_t)o2::aod::track::TPC) {
        if ((track.trackCutFlag() & CtrackSelectionTPC) != CtrackSelectionTPC) {
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
        if (produceHistos) {
          registry.fill(HIST("PropagationFailures"), 1);
        }
        auto compatibleColls = track.compatibleColl();
        int failures = 0;
        for (auto const& collision : compatibleColls) {
          auto propagated = o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackPar, 2.f, matCorr, &dcaInfo);
          if (!propagated) {
            ++failures;
          }
          if (propagated && ((std::abs(dcaInfo[0]) < std::abs(bestDCA[0])) && (std::abs(dcaInfo[1]) < std::abs(bestDCA[1])))) {
            bestCol = collision.globalIndex();
            bestDCA[0] = dcaInfo[0];
            bestDCA[1] = dcaInfo[1];
            bestTrackPar = trackPar;
          }
        }
        if (produceHistos) {
          switch (failures) {
            case 0:
              registry.fill(HIST("PropagationFailures"), 2);
              break;
            case 1:
              registry.fill(HIST("PropagationFailures"), 3);
              break;
            case 2:
              registry.fill(HIST("PropagationFailures"), 4);
              break;
            default:
              registry.fill(HIST("PropagationFailures"), 5);
              break;
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
    float dcaInfo = 0.f;
    float bestDCA = 0.f, bestDCAx = 0.f, bestDCAy = 0.f;
    o2::track::TrackParCovFwd bestTrackPar;

    for (auto const& atrack : atracks) {
      dcaInfo = 999; // DCAxy
      bestDCA = 999;

      auto track = atrack.mfttrack();
      auto bestCol = track.has_collision() ? track.collisionId() : -1;

      o2::track::TrackParCovFwd trackPar = o2::aod::fwdtrackutils::getTrackParCovFwdShift(track, mZShift);

      int degree = 0; // degree of ambiguity of the track

      auto compatibleBCs = atrack.bc_as<ExtBCs>();
      for (auto const& bc : compatibleBCs) {
        if (!bc.has_collisions()) {
          continue;
        }
        auto collisions = bc.collisions();
        for (auto const& collision : collisions) {
          degree++;
          trackPar.propagateToZhelix(collision.posZ(), bZ); // track parameters propagation to the position of the z vertex

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

    float dcaInfo = 0.f;
    float bestDCA = 0.f, bestDCAx = 0.f, bestDCAy = 0.f;
    o2::track::TrackParCovFwd bestTrackPar;

    for (auto const& track : tracks) {
      dcaInfo = 999; // DCAxy
      bestDCA = 999;

      auto bestCol = track.has_collision() ? track.collisionId() : -1;

      // auto ids = track.compatibleCollIds();
      // if (ids.empty() || (ids.size() == 1 && bestCol == ids[0]))
      // {
      //   continue;
      // }

      auto compatibleColls = track.compatibleColl();

      o2::track::TrackParCovFwd trackPar = o2::aod::fwdtrackutils::getTrackParCovFwdShift(track, mZShift);

      for (auto const& collision : compatibleColls) {

        trackPar.propagateToZhelix(collision.posZ(), bZ); // track parameters propagation to the position of the z vertex

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

  void processMFTReassoc3D(MFTTracksWColls const& tracks, aod::Collisions const&, ExtBCs const& bcs)
  {
    if (bcs.size() == 0) {
      return;
    }
    if (tracks.size() == 0) {
      return;
    }
    auto bc = bcs.begin();
    initCCDB(bc);

    std::array<double, 3> dcaInfOrig;
    std::array<double, 2> dcaInfo;
    double bestDCA[2];
    o2::track::TrackParCovFwd bestTrackPar;

    for (auto const& track : tracks) {
      dcaInfOrig[0] = 999.f; // original DCAx from propagation
      dcaInfOrig[1] = 999.f; // original DCAy from propagation
      dcaInfOrig[2] = 999.f; // original DCAz from propagation
      dcaInfo[0] = 999.f;    // calcualted DCAxy
      dcaInfo[1] = 999.f;    // calculated DCAz - same as original
      bestDCA[0] = 999.f;    // minimal DCAxy
      bestDCA[1] = 999.f;    // minimal DCAz

      auto bestCol = track.has_collision() ? track.collisionId() : -1;

      if (removeTrivialAssoc) {
        if (track.compatibleCollIds().empty() || (track.compatibleCollIds().size() == 1 && bestCol == track.compatibleCollIds()[0])) {
          if (produceHistos) {
            registry.fill(HIST("TracksAmbDegreeWoTrivial"), track.compatibleCollIds().size());
          }
          continue;
        }
      }

      auto compatibleColls = track.compatibleColl();

      o2::track::TrackParCovFwd trackPar = o2::aod::fwdtrackutils::getTrackParCovFwdShift(track, mZShift);

      for (auto const& collision : compatibleColls) {

        trackPar.propagateToDCAhelix(bZ, {collision.posX(), collision.posY(), collision.posZ()}, dcaInfOrig);

        if (cfgDCAtype == 0) {
          dcaInfo[0] = dcaInfOrig[0];
        } else if (cfgDCAtype == 1) {
          dcaInfo[0] = dcaInfOrig[1];
        } else if (cfgDCAtype == 2) {
          dcaInfo[0] = std::sqrt(dcaInfOrig[0] * dcaInfOrig[0] + dcaInfOrig[1] * dcaInfOrig[1]);
        }
        dcaInfo[1] = dcaInfOrig[2];

        if ((std::abs(dcaInfo[0]) < std::abs(bestDCA[0])) && (std::abs(dcaInfo[1]) < std::abs(bestDCA[1]))) {
          bestCol = collision.globalIndex();
          bestDCA[0] = dcaInfo[0];
          bestDCA[1] = dcaInfo[1];
          bestTrackPar = trackPar;
        }

        if ((track.collisionId() != collision.globalIndex()) && produceHistos) {
          registry.fill(HIST("DeltaZ"), track.collision().posZ() - collision.posZ()); // deltaZ between the 1st coll zvtx and the other compatible ones
          registry.fill(HIST("DeltaX"), track.collision().posX() - collision.posX());
          registry.fill(HIST("DeltaY"), track.collision().posY() - collision.posY());
          registry.fill(HIST("TracksFirstDCAXY"), dcaInfo[0]);
          registry.fill(HIST("TracksFirstDCAZ"), dcaInfo[1]);
        }
        if (produceHistos) {
          registry.fill(HIST("TracksDCAXY"), dcaInfo[0]);
          registry.fill(HIST("TracksDCAZ"), dcaInfo[1]);
        }
        if ((collision.globalIndex() == track.collisionId()) && produceHistos) {
          registry.fill(HIST("TracksOrigDCAXY"), dcaInfo[0]);
          registry.fill(HIST("TracksOrigDCAZ"), dcaInfo[1]);
        }
      }
      if ((bestCol == track.collisionId()) && produceHistos) {
        registry.fill(HIST("AssignedDCAXY"), bestDCA[0]);
        registry.fill(HIST("AssignedDCAZ"), bestDCA[1]);
      }
      if ((bestCol != track.collisionId()) && produceHistos) {
        // reassigned
        registry.fill(HIST("ReassignedDCAXY"), bestDCA[0]);
        registry.fill(HIST("ReassignedDCAZ"), bestDCA[1]);
      }
      if (produceHistos) {
        int isAmbiguous = 0;
        registry.fill(HIST("TracksAmbDegree"), compatibleColls.size());
        if (compatibleColls.size() > 1) {
          isAmbiguous = 1;
        }
        registry.fill(HIST("TrackIsAmb"), isAmbiguous);
      }

      fwdtracksBestCollisions3d(track.globalIndex(), compatibleColls.size(), bestCol, bestDCA[0], bestDCA[1]);
      // LOGP(info, "track {}: {} {} {} {}", track.globalIndex(), compatibleColls.size(), bestCol, bestDCA[0], bestDCA[1]);
      if (produceExtra) {
        // LOGP(info, "track {}: {} {} {} {} {}", track.globalIndex(), bestTrackPar.getX(), bestTrackPar.getY(), bestTrackPar.getZ(), bestTrackPar.getTgl(), bestTrackPar.getInvQPt());
        fwdtracksBestCollisions3dExtra(bestTrackPar.getX(), bestTrackPar.getY(), bestTrackPar.getZ(),
                                       bestTrackPar.getTgl(), bestTrackPar.getInvQPt(), bestTrackPar.getPt(),
                                       bestTrackPar.getP(), bestTrackPar.getEta(), bestTrackPar.getPhi());
      }
    }
  }
  PROCESS_SWITCH(AmbiguousTrackPropagation, processMFTReassoc3D, "Fill ReassignedTracks for MFT ambiguous tracks", false);
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
