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

/// \file muonRealignment.cxx
/// \brief Task for muon re-alignment at analysis level
/// \author Chi Zhang <chi.zhang@cern.ch>, CEA-Saclay

#include <filesystem>
#include <string>
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CCDBTimeStampUtils.h"
#include "CommonUtils/NameConf.h"
#include "CommonUtils/ConfigurableParam.h"
#include "DataFormatsMCH/Cluster.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "DetectorsBase/Propagator.h"
#include "MathUtils/Cartesian.h"
#include "MCHGeometryTransformer/Transformations.h"
#include "MCHTracking/Track.h"
#include "MCHTracking/TrackExtrap.h"
#include "MCHTracking/TrackParam.h"
#include "MCHTracking/TrackFitter.h"
#include "MCHBase/TrackerParam.h"
#include "GlobalTracking/MatchGlobalFwd.h"
#include "ReconstructionDataFormats/TrackFwd.h"
#include "Common/DataModel/FwdTrackReAlignTables.h"

using namespace std;
using namespace o2;
using namespace o2::framework;
using namespace o2::mch;
using namespace o2::framework::expressions;

const int fgNDetElemCh[10] = {4, 4, 4, 4, 18, 18, 26, 26, 26, 26};
const int fgSNDetElemCh[11] = {0, 4, 8, 12, 16, 34, 52, 78, 104, 130, 156};

struct MuonRealignment {
  Produces<aod::StoredFwdTracksReAlign> realignFwdTrks;
  Produces<aod::StoredFwdTrksCovReAlign> realignFwdTrksCov;

  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> geoRefPath{"geoRefPath", "GLO/Config/GeometryAligned", "Path of the reference geometry file"};
  Configurable<std::string> geoNewPath{"geoNewPath", "GLO/Config/GeometryAligned", "Path of the new geometry file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> grpPathLocal{"grpPathLocal", "", "Local path of the GRP object if not using CCDB"};
  Configurable<std::string> geoNewPathLocal{"geoNewPathLocal", "", "Local path of the GRP object if not using CCDB"};
  Configurable<int64_t> nolaterthanRef{"ccdb-no-later-than-ref", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object of reference basis"};
  Configurable<int64_t> nolaterthanNew{"ccdb-no-later-than-new", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object of new basis"};
  Configurable<double> cfgChamberResolutionX{"cfgChamberResolutionX", 0.04, "Chamber resolution along X configuration for refit"}; // 0.4cm pp, 0.2cm PbPb
  Configurable<double> cfgChamberResolutionY{"cfgChamberResolutionY", 0.04, "Chamber resolution along Y configuration for refit"}; // 0.4cm pp, 0.2cm PbPb
  Configurable<double> cfgSigmaCutImprove{"cfgSigmaCutImprove", 6., "Sigma cut for track improvement"};                            // 6 for pp, 4 for PbPb
  Configurable<bool> fUseRemoteField{"cfgUseRemoteField", true, "Chose whether to fetch the magnetic field from ccdb or set it manually"};
  Configurable<bool> fUseRemoteGeometry{"cfgUseRemoteGeometry", false, "Chose whether to fetch new geometry from ccdb or set it manually"};

  parameters::GRPMagField* grpmag = nullptr;
  base::MatLayerCylSet* lut = nullptr;
  TrackFitter trackFitter; // Track fitter from MCH tracking library
  geo::TransformationCreator transformation;
  map<int, math_utils::Transform3D> transformRef; // reference geometry w.r.t track data
  map<int, math_utils::Transform3D> transformNew; // new geometry
  globaltracking::MatchGlobalFwd mMatching;
  int fCurrentRun;        // needed to detect if the run changed and trigger update of calibrations etc.
  double mImproveCutChi2; // Chi2 cut for track improvement.
  Service<ccdb::BasicCCDBManager> ccdb;
  TGeoManager* geoNew = nullptr;
  TGeoManager* geoRef = nullptr;

  Preslice<aod::FwdTrkCl> perMuon = aod::fwdtrkcl::fwdtrackId;

  int GetDetElemId(int iDetElemNumber)
  {
    // make sure detector number is valid
    if (!(iDetElemNumber >= fgSNDetElemCh[0] &&
          iDetElemNumber < fgSNDetElemCh[10])) {
      LOGF(fatal, "Invalid detector element number: %d", iDetElemNumber);
    }
    /// get det element number from ID
    // get chamber and element number in chamber
    int iCh = 0;
    int iDet = 0;
    for (int i = 1; i <= 10; i++) {
      if (iDetElemNumber < fgSNDetElemCh[i]) {
        iCh = i;
        iDet = iDetElemNumber - fgSNDetElemCh[i - 1];
        break;
      }
    }

    // make sure detector index is valid
    if (!(iCh > 0 && iCh <= 10 && iDet < fgNDetElemCh[iCh - 1])) {
      LOGF(fatal, "Invalid detector element id: %d", 100 * iCh + iDet);
    }

    // add number of detectors up to this chamber
    return 100 * iCh + iDet;
  }

  bool RemoveTrack(mch::Track& track)
  {
    // Refit track with re-aligned clusters
    bool removeTrack = false;
    try {
      trackFitter.fit(track, false);
    } catch (exception const& e) {
      removeTrack = true;
      return removeTrack;
    }

    auto itStartingParam = std::prev(track.rend());

    while (true) {

      try {
        trackFitter.fit(track, true, false, (itStartingParam == track.rbegin()) ? nullptr : &itStartingParam);
      } catch (exception const&) {
        removeTrack = true;
        break;
      }

      double worstLocalChi2 = -1.0;

      track.tagRemovableClusters(0x1F, false);

      auto itWorstParam = track.end();

      for (auto itParam = track.begin(); itParam != track.end(); ++itParam) {
        if (itParam->getLocalChi2() > worstLocalChi2) {
          worstLocalChi2 = itParam->getLocalChi2();
          itWorstParam = itParam;
        }
      }

      if (worstLocalChi2 < mImproveCutChi2) {
        break;
      }

      if (!itWorstParam->isRemovable()) {
        removeTrack = true;
        track.removable();
        break;
      }

      auto itNextParam = track.removeParamAtCluster(itWorstParam);
      auto itNextToNextParam = (itNextParam == track.end()) ? itNextParam : std::next(itNextParam);
      itStartingParam = track.rbegin();

      if (track.getNClusters() < 10) {
        removeTrack = true;
        break;
      } else {
        while (itNextToNextParam != track.end()) {
          if (itNextToNextParam->getClusterPtr()->getChamberId() != itNextParam->getClusterPtr()->getChamberId()) {
            itStartingParam = std::make_reverse_iterator(++itNextParam);
            break;
          }
          ++itNextToNextParam;
        }
      }
    }

    if (!removeTrack) {
      for (auto& param : track) {
        param.setParameters(param.getSmoothParameters());
        param.setCovariances(param.getSmoothCovariances());
      }
    }

    return removeTrack;
  }

  void init(InitContext const&)
  {
    fCurrentRun = 0;

    // Configuration for CCDB server
    ccdb->setURL(ccdburl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    // Configuration for track fitter
    const auto& trackerParam = TrackerParam::Instance();
    trackFitter.setBendingVertexDispersion(trackerParam.bendingVertexDispersion);
    trackFitter.setChamberResolution(cfgChamberResolutionX.value, cfgChamberResolutionY.value);
    trackFitter.smoothTracks(true);
    trackFitter.useChamberResolution();
    mImproveCutChi2 = 2. * cfgSigmaCutImprove.value * cfgSigmaCutImprove.value;
  }

  void process(aod::Collision const& collision, aod::FwdTracks const& tracks, aod::FwdTrkCls const& clusters, aod::BCsWithTimestamps const&)
  {
    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    if (fCurrentRun != bc.runNumber()) {
      // Load magnetic field information from CCDB/local
      if (fUseRemoteField) {
        ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
        grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, bc.timestamp());
        if (grpmag != nullptr) {
          base::Propagator::initFieldFromGRP(grpmag);
          TrackExtrap::setField();
          TrackExtrap::useExtrapV2();
        } else {
          LOGF(fatal, "GRP object is not available in CCDB at timestamp=%llu", bc.timestamp());
        }
      } else {
        if (std::filesystem::exists(grpPathLocal.value)) {
          const auto grp = parameters::GRPObject::loadFrom(grpPathLocal.value);
          base::Propagator::initFieldFromGRP(grp);
          TrackExtrap::setField();
          TrackExtrap::useExtrapV2();
        } else {
          LOGF(fatal, "GRP object is not available in local path: %s", grpPathLocal.value);
        }
      }

      // Load geometry information from CCDB/local
      LOGF(info, "Loading reference aligned geometry from CCDB no later than %d", nolaterthanRef.value);
      ccdb->setCreatedNotAfter(nolaterthanRef.value); // this timestamp has to be consistent with what has been used in reco
      geoRef = ccdb->getForTimeStamp<TGeoManager>(geoRefPath, bc.timestamp());
      ccdb->clearCache(geoRefPath);
      if (geoRef != nullptr) {
        transformation = geo::transformationFromTGeoManager(*geoRef);
      } else {
        LOGF(fatal, "Reference aligned geometry object is not available in CCDB at timestamp=%llu", bc.timestamp());
      }
      for (int i = 0; i < 156; i++) {
        int iDEN = GetDetElemId(i);
        transformRef[iDEN] = transformation(iDEN);
      }

      if (fUseRemoteGeometry) {
        LOGF(info, "Loading new aligned geometry from CCDB no later than %d", nolaterthanNew.value);
        ccdb->setCreatedNotAfter(nolaterthanNew.value); // make sure this timestamp can be resolved regarding the reference one
        geoNew = ccdb->getForTimeStamp<TGeoManager>(geoNewPath, bc.timestamp());
        ccdb->clearCache(geoNewPath);
        if (geoNew != nullptr) {
          transformation = geo::transformationFromTGeoManager(*geoNew);
        } else {
          LOGF(fatal, "New aligned geometry object is not available in CCDB at timestamp=%llu", bc.timestamp());
        }
        for (int i = 0; i < 156; i++) {
          int iDEN = GetDetElemId(i);
          transformNew[iDEN] = transformation(iDEN);
        }
      } else {
        LOGF(info, "Loading new aligned geometry from local path: %s", geoNewPathLocal.value);
        if (std::filesystem::exists(geoNewPathLocal.value)) {
          base::GeometryManager::loadGeometry(geoNewPathLocal.value);
          transformation = geo::transformationFromTGeoManager(*gGeoManager);
          for (int i = 0; i < 156; i++) {
            int iDEN = GetDetElemId(i);
            transformNew[iDEN] = transformation(iDEN);
          }
        } else {
          LOGF(fatal, "New geometry file is not available in local path: %s", geoNewPathLocal.value);
        }
      }

      fCurrentRun = bc.runNumber();
    }

    // Reserve storage for output table
    realignFwdTrks.reserve(tracks.size());
    realignFwdTrksCov.reserve(tracks.size());

    // Loop over forward tracks
    for (auto const& track : tracks) {
      if (track.has_collision()) {
        if (track.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack || track.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::MCHStandaloneTrack) {

          auto clustersSliced = clusters.sliceBy(perMuon, track.globalIndex()); // Slice clusters by muon id
          mch::Track convertedTrack = mch::Track();                             // Temporary variable to store re-aligned clusters
          int clIndex = -1;
          // Get re-aligned clusters associated to current track
          for (auto const& cluster : clustersSliced) {
            clIndex += 1;

            mch::Cluster* clusterMCH = new mch::Cluster();

            math_utils::Point3D<double> local;
            math_utils::Point3D<double> master;
            master.SetXYZ(cluster.x(), cluster.y(), cluster.z());

            // Transformation from reference geometry frame to new geometry frame
            transformRef[cluster.deId()].MasterToLocal(master, local);
            transformNew[cluster.deId()].LocalToMaster(local, master);

            clusterMCH->x = master.x();
            clusterMCH->y = master.y();
            clusterMCH->z = master.z();

            uint32_t ClUId = mch::Cluster::buildUniqueId(static_cast<int>(cluster.deId() / 100) - 1, cluster.deId(), clIndex);
            clusterMCH->uid = ClUId;
            clusterMCH->ex = cluster.isGoodX() ? 0.2 : 10.0;
            clusterMCH->ey = cluster.isGoodY() ? 0.2 : 10.0;

            // Add transformed cluster into temporary variable
            convertedTrack.createParamAtCluster(*clusterMCH);
            LOGF(debug, "Track %d, cluster DE%d:  x:%g  y:%g  z:%g", track.globalIndex(), cluster.deId(), cluster.x(), cluster.y(), cluster.z());
            LOGF(debug, "Track %d, re-aligned cluster DE%d:  x:%g  y:%g  z:%g", track.globalIndex(), cluster.deId(), clusterMCH->getX(), clusterMCH->getY(), clusterMCH->getZ());
          }

          // Refit the re-aligned track
          int removable = 0;
          if (convertedTrack.getNClusters() != 0) {
            removable = RemoveTrack(convertedTrack);
          } else {
            LOGF(fatal, "Muon track %d has no associated clusters.", track.globalIndex());
          }

          // Get the re-aligned track parameter: track param at the first cluster
          mch::TrackParam trackParam = mch::TrackParam(convertedTrack.first());

          // Convert MCH track to FWD track and get new parameters
          auto fwdtrack = mMatching.MCHtoFwd(trackParam);
          fwdtrack.setTrackChi2(trackParam.getTrackChi2() / convertedTrack.getNDF());
          float sigX = TMath::Sqrt(fwdtrack.getCovariances()(0, 0));
          float sigY = TMath::Sqrt(fwdtrack.getCovariances()(1, 1));
          float sigPhi = TMath::Sqrt(fwdtrack.getCovariances()(2, 2));
          float sigTgl = TMath::Sqrt(fwdtrack.getCovariances()(3, 3));
          float sig1Pt = TMath::Sqrt(fwdtrack.getCovariances()(4, 4));
          float rhoXY = (Char_t)(128. * fwdtrack.getCovariances()(0, 1) / (sigX * sigY));
          float rhoPhiX = (Char_t)(128. * fwdtrack.getCovariances()(0, 2) / (sigPhi * sigX));
          float rhoPhiY = (Char_t)(128. * fwdtrack.getCovariances()(1, 2) / (sigPhi * sigY));
          float rhoTglX = (Char_t)(128. * fwdtrack.getCovariances()(0, 3) / (sigTgl * sigX));
          float rhoTglY = (Char_t)(128. * fwdtrack.getCovariances()(1, 3) / (sigTgl * sigY));
          float rhoTglPhi = (Char_t)(128. * fwdtrack.getCovariances()(2, 3) / (sigTgl * sigPhi));
          float rho1PtX = (Char_t)(128. * fwdtrack.getCovariances()(0, 4) / (sig1Pt * sigX));
          float rho1PtY = (Char_t)(128. * fwdtrack.getCovariances()(1, 4) / (sig1Pt * sigY));
          float rho1PtPhi = (Char_t)(128. * fwdtrack.getCovariances()(2, 4) / (sig1Pt * sigPhi));
          float rho1PtTgl = (Char_t)(128. * fwdtrack.getCovariances()(3, 4) / (sig1Pt * sigTgl));

          LOGF(debug, "TrackParm %d, x:%g  y:%g  z:%g  phi:%g  tgl:%g  InvQPt:%g  chi2:%g  nClusters:%d", track.globalIndex(), track.x(), track.y(), track.z(), track.phi(), track.tgl(), track.signed1Pt(), track.chi2(), track.nClusters());
          LOGF(debug, "Re-aligned trackParm %d, x:%g  y:%g  z:%g  phi:%g  tgl:%g  InvQPt:%g  chi2:%g  nClusters:%d  removable:%d", track.globalIndex(), fwdtrack.getX(), fwdtrack.getY(), fwdtrack.getZ(), fwdtrack.getPhi(), fwdtrack.getTgl(), fwdtrack.getInvQPt(), fwdtrack.getTrackChi2(), convertedTrack.getNClusters(), removable);
          // Fill refitted track info
          realignFwdTrks(fwdtrack.getX(), fwdtrack.getY(), fwdtrack.getZ(), fwdtrack.getPhi(), fwdtrack.getTgl(), fwdtrack.getInvQPt(), fwdtrack.getTrackChi2(), removable);
          realignFwdTrksCov(sigX, sigY, sigPhi, sigTgl, sig1Pt, rhoXY, rhoPhiX, rhoPhiY, rhoTglX, rhoTglY, rhoTglPhi, rho1PtX, rho1PtY, rho1PtPhi, rho1PtTgl);
        } else {
          // Fill nothing for global muons
          realignFwdTrks(-999., -999., -999., -999., -999., -999., -999., -999.);
          realignFwdTrksCov(-999., -999., -999., -999., -999., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        }
      } else {
        // Fill nothing for tracks having no associated collision
        realignFwdTrks(-999., -999., -999., -999., -999., -999., -999., -999.);
        realignFwdTrksCov(-999., -999., -999., -999., -999., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<MuonRealignment>(cfgc)};
}
