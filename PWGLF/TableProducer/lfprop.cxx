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
/// \file lfprop.h
/// \since 20-11-2023
/// \author Carolina Anna Reetz <carolina.reetz@cern.ch>
/// \author Nicolò Jacazio <nicolo.jacazio@cern.ch>
/// \brief Task to study the propagation of tracks from the TPC inner wall to the collision vertex
///

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
#include "Common/DataModel/PIDResponse.h"
#include "TPDGCode.h"

using namespace o2;
using namespace o2::framework;

struct LfProp {

  Produces<aod::StoredTracks> tracksParPropagated;
  Produces<aod::TracksExtension> tracksParExtensionPropagated;

  Produces<aod::StoredTracksCov> tracksParCovPropagated;
  Produces<aod::TracksCovExtension> tracksParCovExtensionPropagated;

  Produces<aod::TracksDCA> tracksDCA;
  Produces<aod::TracksDCACov> tracksDCACov;

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;

  o2::parameters::GRPMagField* grpmag = nullptr;
  o2::base::MatLayerCylSet* lut = nullptr;

  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> mVtxPath{"mVtxPath", "GLO/Calib/MeanVertex", "Path of the mean vertex file"};
  Configurable<float> minPropagationRadius{"minPropagationDistance", o2::constants::geom::XTPCInnerRef + 0.1, "Only tracks which are at a smaller radius will be propagated, defaults to TPC inner wall"};

  int runNumber = 0;
  void init(o2::framework::InitContext& initContext)
  {
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
    // mVtx = ccdb->getForTimeStamp<o2::dataformats::MeanVertexObject>(mVtxPath, bc.timestamp());
    runNumber = bc.runNumber();
  }

  template <typename particleType>
  o2::track::PID::ID pdgToIndex(const particleType& p)
  {
    switch (std::abs(p.pdgCode())) {
      case kElectron:
        return o2::track::PID::Electron;
      case kMuonMinus:
        return o2::track::PID::Muon;
      case kPiPlus:
        return o2::track::PID::Pion;
      case kKPlus:
        return o2::track::PID::Kaon;
      case kProton:
        return o2::track::PID::Proton;
      case 1000010020:
        return o2::track::PID::Deuteron;
      case 1000010030:
        return o2::track::PID::Triton;
      case 1000020030:
        return o2::track::PID::Helium3;
      case 1000020040:
        return o2::track::PID::Alpha;
      default:
        return o2::track::PID::NIDs;
    }
  }

  void process(soa::Join<aod::StoredTracksIU, aod::TracksCovIU, aod::McTrackLabels> const& tracks,
               aod::Collisions const&,
               aod::McParticles const&,
               aod::BCsWithTimestamps const& bcs)
  {
    if (bcs.size() == 0) {
      return;
    }
    initCCDB(bcs.begin());

    o2::dataformats::DCA dcaInfoCov;
    o2::dataformats::VertexBase vtx;

    tracksParPropagated.reserve(tracks.size());
    tracksParExtensionPropagated.reserve(tracks.size());
    tracksParCovPropagated.reserve(tracks.size());
    tracksParCovExtensionPropagated.reserve(tracks.size());
    tracksDCA.reserve(tracks.size());
    tracksDCACov.reserve(tracks.size());

    for (auto& track : tracks) {
      dcaInfoCov.set(999, 999, 999, 999, 999);
      auto trackParCov = getTrackParCov(track);
      if (track.has_mcParticle()) {
        trackParCov.setPID(pdgToIndex(track.mcParticle()));
      }
      aod::track::TrackTypeEnum trackType = (aod::track::TrackTypeEnum)track.trackType();
      // Only propagate tracks which have passed the innermost wall of the TPC (e.g. skipping loopers etc). Others fill unpropagated.
      if (track.trackType() == aod::track::TrackIU && track.x() < minPropagationRadius) {
        if (track.has_collision()) {
          auto const& collision = track.collision();
          vtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
          vtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
          trackType = aod::track::Track;
          o2::base::Propagator::Instance()->propagateToDCABxByBz(vtx, trackParCov, 2.f, matCorr, &dcaInfoCov);
        } else {
          //   vtx.setPos({mVtx->getX(), mVtx->getY(), mVtx->getZ()});
          //   vtx.setCov(mVtx->getSigmaX() * mVtx->getSigmaX(), 0.0f, mVtx->getSigmaY() * mVtx->getSigmaY(), 0.0f, 0.0f, mVtx->getSigmaZ() * mVtx->getSigmaZ());
        }
      }
      tracksParPropagated(track.collisionId(), trackType, trackParCov.getX(), trackParCov.getAlpha(), trackParCov.getY(), trackParCov.getZ(), trackParCov.getSnp(), trackParCov.getTgl(), trackParCov.getQ2Pt());
      tracksParExtensionPropagated(trackParCov.getPt(), trackParCov.getP(), trackParCov.getEta(), trackParCov.getPhi());
      tracksDCA(dcaInfoCov.getY(), dcaInfoCov.getZ());
      tracksDCACov(dcaInfoCov.getSigmaY2(), dcaInfoCov.getSigmaZ2());
      // TODO do we keep the rho as 0? Also the sigma's are duplicated information
      tracksParCovPropagated(std::sqrt(trackParCov.getSigmaY2()), std::sqrt(trackParCov.getSigmaZ2()), std::sqrt(trackParCov.getSigmaSnp2()),
                             std::sqrt(trackParCov.getSigmaTgl2()), std::sqrt(trackParCov.getSigma1Pt2()), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
      tracksParCovExtensionPropagated(trackParCov.getSigmaY2(), trackParCov.getSigmaZY(), trackParCov.getSigmaZ2(), trackParCov.getSigmaSnpY(),
                                      trackParCov.getSigmaSnpZ(), trackParCov.getSigmaSnp2(), trackParCov.getSigmaTglY(), trackParCov.getSigmaTglZ(), trackParCov.getSigmaTglSnp(),
                                      trackParCov.getSigmaTgl2(), trackParCov.getSigma1PtY(), trackParCov.getSigma1PtZ(), trackParCov.getSigma1PtSnp(), trackParCov.getSigma1PtTgl(),
                                      trackParCov.getSigma1Pt2());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<LfProp>(cfgc)}; }
