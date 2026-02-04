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

/// \file createPCM.cxx
/// \brief This code produces photon data tables.
/// \author Daiki Sekihata <daiki.sekihata@cern.ch>, Tokyo

#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/PCMUtilities.h"
#include "PWGEM/PhotonMeson/Utils/TrackSelection.h"
//
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <DCAFitter/DCAFitterN.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DataFormatsParameters/GRPObject.h>
#include <DetectorsBase/GeometryManager.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/runDataProcessing.h>

#include <algorithm>
#include <array>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::pwgem::photonmeson;

using FullTracksExtIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA>;
// using FullTracksExt = soa::Join<aod::Tracks, aod::TracksCov, aod::TracksExtra, aod::TracksDCA>;
using FullTrackExtIU = FullTracksExtIU::iterator;

struct createPCM {
  SliceCache cache;
  Preslice<aod::TracksIU> perCol = o2::aod::track::collisionId;
  Produces<aod::V0Indices> v0indices;
  Produces<aod::V0CoresBase> v0cores;
  Produces<aod::V0TrackXs> v0trackXs;

  // Basic checks
  HistogramRegistry registry{
    "createPCM",
    {
      {"hEventCounter", "hEventCounter", {HistType::kTH1F, {{5, 0.5f, 5.5f}}}},
      {"hV0xy", "hV0xy;X (cm);Y(cm)", {HistType::kTH2F, {{400, -100, +100}, {400, -100, +100}}}},
      {"hV0xy_recalculated", "hV0xy_recalculated;X (cm);Y(cm)", {HistType::kTH2F, {{400, -100, +100}, {400, -100, +100}}}},
    },
  };

  // Configurables
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  Configurable<float> d_bz_input{"d_bz_input", -999, "bz field in kG, -999 is automatic"};
  Configurable<bool> d_UseAbsDCA{"d_UseAbsDCA", true, "Use Abs DCAs"};
  Configurable<bool> d_UseWeightedPCA{"d_UseWeightedPCA", false, "Vertices use cov matrices"};
  Configurable<int> useMatCorrType{"useMatCorrType", 0, "0: none, 1: TGeo, 2: LUT"};

  Configurable<float> minv0cospa{"minv0cospa", 0.90, "minimum V0 CosPA"};
  Configurable<float> maxdcav0dau{"maxdcav0dau", 2.0, "max DCA between V0 Daughters"};
  Configurable<float> v0Rmin{"v0Rmin", 0.0, "v0Rmin"};
  Configurable<float> v0Rmax{"v0Rmax", 180.0, "v0Rmax"};
  Configurable<float> dcamin{"dcamin", 0.1, "dcamin"};
  Configurable<float> dcamax{"dcamax", 1e+10, "dcamax"};
  Configurable<int> nsw{"nsw", 1, "number of searching window in collisions"};
  Configurable<float> maxX{"maxX", 83.1, "maximum X (starting point X of track iu)"};
  Configurable<float> maxY{"maxY", 20.0, "maximum Y (starting point Y of track iu)"};
  Configurable<float> minpt{"minpt", 0.01, "min pT for single track in GeV/c"};
  Configurable<float> maxeta{"maxeta", 0.9, "eta acceptance for single track"};
  Configurable<int> mincrossedrows{"mincrossedrows", 10, "min crossed rows"};
  Configurable<float> maxchi2tpc{"maxchi2tpc", 4.0, "max chi2/NclsTPC"};
  Configurable<float> maxchi2its{"maxchi2its", 5.0, "max chi2/NclsITS"};
  Configurable<float> maxpt_itsonly{"maxpt_itsonly", 0.5, "max pT for ITSonly tracks"};
  Configurable<float> min_tpcdEdx{"min_tpcdEdx", 30.0, "min TPC dE/dx"};
  Configurable<float> max_tpcdEdx{"max_tpcdEdx", 110.0, "max TPC dE/dx"};
  Configurable<float> margin_r{"margin_r", 7.0, "margin for r cut"};
  Configurable<float> max_qt_arm{"max_qt_arm", 0.03, "max qt for AP cut in GeV/c"};
  Configurable<float> max_r_req_its{"max_r_req_its", 16.0, "min Rxy for V0 with ITS hits"};
  Configurable<float> min_r_tpconly{"min_r_tpconly", 32.0, "min Rxy for V0 with TPConly tracks"};

  int mRunNumber;
  float d_bz;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut = nullptr;
  o2::vertexing::DCAFitterN<2> fitter;
  // Material correction in the DCA fitter
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  float calculateDCAStraightToPV(float X, float Y, float Z, float Px, float Py, float Pz, float pvX, float pvY, float pvZ)
  {
    return std::sqrt((std::pow((pvY - Y) * Pz - (pvZ - Z) * Py, 2) + std::pow((pvX - X) * Pz - (pvZ - Z) * Px, 2) + std::pow((pvX - X) * Py - (pvY - Y) * Px, 2)) / (Px * Px + Py * Py + Pz * Pz));
  }

  void init(InitContext&)
  {
    mRunNumber = 0;
    d_bz = 0;

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    if (useMatCorrType == 1) {
      LOGF(info, "TGeo correction requested, loading geometry");
      if (!o2::base::GeometryManager::isGeometryLoaded()) {
        ccdb->get<TGeoManager>(geoPath);
      }
    }
    if (useMatCorrType == 2) {
      LOGF(info, "LUT correction requested, loading LUT");
      lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(lutPath));
    }

    fitter.setPropagateToPCA(true);
    fitter.setMaxR(200.);
    fitter.setMinParamChange(1e-3);
    fitter.setMinRelChi2Change(0.9);
    fitter.setMaxDZIni(1e9);
    fitter.setMaxChi2(1e9);
    fitter.setUseAbsDCA(d_UseAbsDCA);
    fitter.setWeightedFinalPCA(d_UseWeightedPCA);

    if (useMatCorrType == 1) {
      matCorr = o2::base::Propagator::MatCorrType::USEMatCorrTGeo;
    }
    if (useMatCorrType == 2) {
      matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
    }
    fitter.setMatCorrType(matCorr);
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    // In case override, don't proceed, please - no CCDB access required
    if (d_bz_input > -990) {
      d_bz = d_bz_input;
      fitter.setBz(d_bz);
      o2::parameters::GRPMagField grpmag;
      if (std::fabs(d_bz) > 1e-5) {
        grpmag.setL3Current(30000.f / (d_bz / 5.0f));
      }
      o2::base::Propagator::initFieldFromGRP(&grpmag);
      mRunNumber = bc.runNumber();
      return;
    }

    auto run3grp_timestamp = bc.timestamp();
    auto* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    o2::parameters::GRPMagField* grpmag = nullptr;
    if (grpo != nullptr) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      // Fetch magnetic field from ccdb for current collision
      d_bz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (grpmag == nullptr) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      // Fetch magnetic field from ccdb for current collision
      d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    }
    mRunNumber = bc.runNumber();
    // Set magnetic field value once known
    fitter.setBz(d_bz);

    if (useMatCorrType == 2) {
      // setMatLUT only after magfield has been initalized
      // (setMatLUT has implicit and problematic init field call if not)
      o2::base::Propagator::Instance()->setMatLUT(lut);
    }
  }

  template <typename TTrack>
  bool reconstructV0(TTrack const& ele, TTrack const& pos)
  {
    bool isITSonly_pos = pos.hasITS() && !pos.hasTPC();
    bool isITSonly_ele = ele.hasITS() && !ele.hasTPC();
    bool isTPConly_pos = !pos.hasITS() && pos.hasTPC();
    bool isTPConly_ele = !ele.hasITS() && ele.hasTPC();

    if ((isITSonly_pos && isTPConly_ele) || (isITSonly_ele && isTPConly_pos)) {
      return false;
    }

    // fitter is memeber variable.
    auto pTrack = getTrackParCov(pos); // positive
    auto nTrack = getTrackParCov(ele); // negative
    std::array<float, 3> svpos = {0.}; // secondary vertex position
    std::array<float, 3> pvec0 = {0.};
    std::array<float, 3> pvec1 = {0.};

    int nCand = fitter.process(pTrack, nTrack);
    if (nCand != 0) {
      fitter.propagateTracksToVertex();
      const auto& vtx = fitter.getPCACandidate();
      for (int i = 0; i < 3; i++) {
        svpos[i] = vtx[i];
      }
      fitter.getTrack(0).getPxPyPzGlo(pvec0); // positive
      fitter.getTrack(1).getPxPyPzGlo(pvec1); // negative
    } else {
      return false;
    }

    float v0dca = std::sqrt(fitter.getChi2AtPCACandidate()); // distance between 2 legs.
    if (v0dca > maxdcav0dau) {
      return false;
    }

    if (!checkAP(v0_alpha(pvec0[0], pvec0[1], pvec0[2], pvec1[0], pvec1[1], pvec1[2]), v0_qt(pvec0[0], pvec0[1], pvec0[2], pvec1[0], pvec1[1], pvec1[2]), 0.95, max_qt_arm)) { // store only photon conversions
      return false;
    }
    if (ele.hasITS() && pos.hasITS() && !checkAP(v0_alpha(pvec0[0], pvec0[1], pvec0[2], pvec1[0], pvec1[1], pvec1[2]), v0_qt(pvec0[0], pvec0[1], pvec0[2], pvec1[0], pvec1[1], pvec1[2]), 0.95, 0.02)) { // store only photon conversions
      return false;
    }

    float xyz[3] = {0.f, 0.f, 0.f};
    Vtx_recalculation(o2::base::Propagator::Instance(), pos, ele, xyz, matCorr);
    float recalculatedVtxR = std::sqrt(std::pow(xyz[0], 2) + std::pow(xyz[1], 2));
    // LOGF(info, "recalculated vtx : x = %f , y = %f , z = %f", xyz[0], xyz[1], xyz[2]);
    if (recalculatedVtxR > std::min(pos.x(), ele.x()) + margin_r && (pos.x() > 1.f && ele.x() > 1.f)) {
      return false;
    }

    if (recalculatedVtxR < max_r_req_its && (!pos.hasITS() || !ele.hasITS())) {
      return false;
    }
    if (recalculatedVtxR < min_r_tpconly && (!pos.hasITS() && !ele.hasITS())) {
      return false;
    }

    return true;
  }

  template <typename TCollision, typename TTrack>
  void fillV0Table(TCollision const& collision, TTrack const& ele, TTrack const& pos, const bool filltable)
  {
    std::array<float, 3> pVtx = {collision.posX(), collision.posY(), collision.posZ()};
    std::array<float, 3> svpos = {0.}; // secondary vertex position
    std::array<float, 3> pvec0 = {0.};
    std::array<float, 3> pvec1 = {0.};

    auto pTrack = getTrackParCov(pos); // positive
    auto nTrack = getTrackParCov(ele); // negative

    int nCand = fitter.process(pTrack, nTrack);
    if (nCand != 0) {
      fitter.propagateTracksToVertex();
      const auto& vtx = fitter.getPCACandidate();
      for (int i = 0; i < 3; i++) {
        svpos[i] = vtx[i];
      }
      fitter.getTrack(0).getPxPyPzGlo(pvec0); // positive
      fitter.getTrack(1).getPxPyPzGlo(pvec1); // negative
    } else {
      return;
    }

    std::array<float, 3> pvxyz{pvec0[0] + pvec1[0], pvec0[1] + pvec1[1], pvec0[2] + pvec1[2]};

    float v0dca = std::sqrt(fitter.getChi2AtPCACandidate()); // distance between 2 legs.
    float v0CosinePA = RecoDecay::cpa(pVtx, svpos, pvxyz);
    float v0radius = RecoDecay::sqrtSumOfSquares(svpos[0], svpos[1]);
    float dcaV0toPV = calculateDCAStraightToPV(svpos[0], svpos[1], svpos[2], pvxyz[0], pvxyz[1], pvxyz[2], pVtx[0], pVtx[1], pVtx[2]);

    if (v0dca > maxdcav0dau) {
      return;
    }
    if (v0CosinePA < minv0cospa) {
      return;
    }

    if (!checkAP(v0_alpha(pvec0[0], pvec0[1], pvec0[2], pvec1[0], pvec1[1], pvec1[2]), v0_qt(pvec0[0], pvec0[1], pvec0[2], pvec1[0], pvec1[1], pvec1[2]), 0.95, max_qt_arm)) { // store only photon conversions
      return;
    }

    if (filltable) {
      if (v0radius < v0Rmin || v0Rmax < v0radius) {
        return;
      }

      registry.fill(HIST("hV0xy"), svpos[0], svpos[1]); // this should have worst resolution
      float xyz_tmp[3] = {0.f, 0.f, 0.f};
      Vtx_recalculation(o2::base::Propagator::Instance(), pos, ele, xyz_tmp, matCorr);
      registry.fill(HIST("hV0xy_recalculated"), xyz_tmp[0], xyz_tmp[1]); // this should have good resolution

      // populates the various tables that comprise V0Datas
      v0indices(pos.globalIndex(), ele.globalIndex(), collision.globalIndex(), -1);
      v0trackXs(fitter.getTrack(0).getX(), fitter.getTrack(1).getX());
      v0cores(svpos[0], svpos[1], svpos[2],
              pvec0[0], pvec0[1], pvec0[2],
              pvec1[0], pvec1[1], pvec1[2],
              v0dca, pos.dcaXY(), ele.dcaXY(),
              v0CosinePA, dcaV0toPV, 3); // v0 type: photon-exclusive
    } else {
      // LOGF(info, "storing: collision.globalIndex() = %d , pos.globalIndex() = %d , ele.globalIndex() = %d, cospa = %f", collision.globalIndex(), pos.globalIndex(), ele.globalIndex(), v0CosinePA);
      pca_map[std::make_tuple(pos.globalIndex(), ele.globalIndex(), collision.globalIndex())] = v0dca;
      cospa_map[std::make_tuple(pos.globalIndex(), ele.globalIndex(), collision.globalIndex())] = v0CosinePA;
    } // store indices
  }

  std::pair<int8_t, std::set<uint8_t>> its_ib_Requirement = {0, {0, 1, 2}}; // no hit on 3 ITS ib layers.
  template <typename TTrack>
  bool isSelected(TTrack const& track)
  {
    if (track.pt() < minpt || std::abs(track.eta()) > maxeta) {
      return false;
    }
    if (std::abs(track.dcaXY()) < dcamin || dcamax < std::abs(track.dcaXY())) {
      return false;
    }
    if (!track.hasITS() && !track.hasTPC()) {
      return false;
    }

    if (track.hasITS() && !track.hasTPC() && (track.hasTRD() || track.hasTOF())) { // remove unrealistic track. this should not happen.
      return false;
    }

    if (track.hasTPC()) {
      if (track.tpcNClsCrossedRows() < mincrossedrows || track.tpcChi2NCl() > maxchi2tpc) {
        return false;
      }
      if (track.tpcSignal() < min_tpcdEdx || max_tpcdEdx < track.tpcSignal()) {
        return false;
      }
    }

    if (track.hasITS()) {
      if (track.itsChi2NCl() > maxchi2its) {
        return false;
      }

      if (std::abs(track.z() / track.x() - track.tgl()) > 0.5) {
        return false;
      }

      auto hits_ib = std::count_if(its_ib_Requirement.second.begin(), its_ib_Requirement.second.end(), [&](auto&& requiredLayer) { return track.itsClusterMap() & (1 << requiredLayer); });
      bool its_ob_only = hits_ib <= its_ib_Requirement.first;
      if (!its_ob_only) {
        return false;
      }

      bool isITSonly = isITSonlyTrack(track);
      if (isITSonly) {
        if (track.pt() > maxpt_itsonly) {
          return false;
        }
      }
    }

    return true;
  }

  Filter trackFilter = o2::aod::track::x < maxX && nabs(o2::aod::track::y) < maxY && o2::aod::track::pt > minpt&& nabs(o2::aod::track::eta) < maxeta&& dcamin < nabs(o2::aod::track::dcaXY) && nabs(o2::aod::track::dcaXY) < dcamax && ((min_tpcdEdx < o2::aod::track::tpcSignal && o2::aod::track::tpcSignal < max_tpcdEdx) || o2::aod::track::tpcSignal < -10.f);
  using MyFilteredTracks = soa::Filtered<FullTracksExtIU>;

  std::map<std::tuple<int32_t, int32_t, int32_t>, float> pca_map;
  std::map<std::tuple<int32_t, int32_t, int32_t>, float> cospa_map;

  // Partition<MyFilteredTracks> orphan_posTracks = o2::aod::track::signed1Pt > 0.f && o2::aod::track::collisionId < int32_t(0);
  // Partition<MyFilteredTracks> orphan_negTracks = o2::aod::track::signed1Pt < 0.f && o2::aod::track::collisionId < int32_t(0);
  Partition<MyFilteredTracks> posTracks = o2::aod::track::signed1Pt > 0.f;
  Partition<MyFilteredTracks> negTracks = o2::aod::track::signed1Pt < 0.f;
  std::vector<decltype(negTracks->sliceByCached(o2::aod::track::collisionId, 0, cache))> negTracks_sw;
  std::vector<decltype(posTracks->sliceByCached(o2::aod::track::collisionId, 0, cache))> posTracks_sw;

  void processSA(MyFilteredTracks const& tracks, aod::Collisions const& collisions, aod::BCsWithTimestamps const&)
  {
    // LOGF(info, "collisions.size() = %d, tracks.size() = %d", collisions.size(), tracks.size());
    for (int64_t icoll = 0; icoll < collisions.size(); icoll += nsw) { // don't repeat the same collision
      auto collision = collisions.rawIteratorAt(icoll);
      // LOGF(info, "collision.globalIndex() = %d", collision.globalIndex());

      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      // registry.fill(HIST("hEventCounter"), 1);

      int32_t min_sw = std::max(static_cast<int64_t>(0), collision.globalIndex());
      int32_t max_sw = std::min(static_cast<int64_t>(min_sw + nsw), static_cast<int64_t>(collisions.size()));

      // LOGF(info, "orphan_posTracks.size() = %d, orphan_negTracks.size() = %d", orphan_posTracks.size(), orphan_negTracks.size());
      negTracks_sw.reserve(max_sw - min_sw);
      posTracks_sw.reserve(max_sw - min_sw);

      // int npos = 0, nneg = 0;
      for (int32_t isw = min_sw; isw < max_sw; isw++) {
        negTracks_sw.emplace_back(negTracks->sliceByCached(o2::aod::track::collisionId, isw, cache));
        posTracks_sw.emplace_back(posTracks->sliceByCached(o2::aod::track::collisionId, isw, cache));
        // npos += posTracks_sw.back().size();
        // nneg += negTracks_sw.back().size();
        // LOGF(info, "collision.globalIndex() = %d , posTracks_sw.back().size() = %d , negTracks_sw.back().size() = %d", collision.globalIndex(), posTracks_sw.back().size(), negTracks_sw.back().size());
      }
      // LOGF(info, "min_sw = %d , max_sw = %d , collision.globalIndex() = %d , n posTracks_sw = %d , n negTracks_sw = %d", min_sw, max_sw, collision.globalIndex(), npos, nneg);

      for (const auto& negTracks_coll : negTracks_sw) {
        for (const auto& posTracks_coll : posTracks_sw) {
          for (const auto& [ele, pos] : combinations(CombinationsFullIndexPolicy(negTracks_coll, posTracks_coll))) {
            if (!isSelected(ele) || !isSelected(pos)) {
              continue;
            }
            if (!reconstructV0(ele, pos)) { // this is needed for speed-up.
              continue;
            }

            for (int32_t isw = min_sw; isw < max_sw; isw++) {
              auto collision_in_sw = collisions.rawIteratorAt(isw);

              if (ele.isPVContributor() && isw != ele.collisionId()) {
                continue;
              }
              if (pos.isPVContributor() && isw != pos.collisionId()) {
                continue;
              }

              // LOGF(info, "pairing: collision_in_sw.globalIndex() = %d , ele.collisionId() = %d , pos.collisionId() = %d ele.globalIndex() = %d , pos.globalIndex() = %d",
              //     collision_in_sw.globalIndex(), ele.collisionId(), pos.collisionId(), ele.globalIndex(), pos.globalIndex());
              fillV0Table(collision_in_sw, ele, pos, false);
            } // end of searching window loop
          } // end of pairing loop
        } // end of pos track loop in sw
      } // end of pos track loop in sw

      // LOGF(info, "possible number of V0 = %d", cospa_map.size());
      std::map<std::pair<uint32_t, uint32_t>, bool> used_pair_map;

      for (const auto& [key, value] : cospa_map) {
        auto pos = tracks.rawIteratorAt(std::get<0>(key));
        auto ele = tracks.rawIteratorAt(std::get<1>(key));

        // LOGF(info, "candidate : pos.globalIndex() = %d , ele.globalIndex() = %d , collision.globalIndex() = %d , cospa = %f , pca = %f", std::get<0>(key), std::get<1>(key), std::get<2>(key), value, pca_map[key]);

        std::vector<float> vec_cospa; // vector for each searching window
        vec_cospa.reserve(max_sw - min_sw);
        for (int32_t isw = min_sw; isw < max_sw; isw++) {
          auto collision_in_sw = collisions.rawIteratorAt(isw);
          if (cospa_map.find(std::make_tuple(pos.globalIndex(), ele.globalIndex(), collision_in_sw.globalIndex())) != cospa_map.end()) {
            vec_cospa.emplace_back(cospa_map[std::make_tuple(pos.globalIndex(), ele.globalIndex(), collision_in_sw.globalIndex())]);
          } else {
            vec_cospa.emplace_back(-999.f);
          }
        } // end of searching window loop

        // search for the most probable collision where V0 belongs by maximal cospa.
        int32_t collision_id_most_prob = std::distance(vec_cospa.begin(), std::max_element(vec_cospa.begin(), vec_cospa.end())) + min_sw;
        auto collision_most_prob = collisions.rawIteratorAt(collision_id_most_prob);
        // float max_cospa = *std::max_element(vec_cospa.begin(), vec_cospa.end());
        // LOGF(info, "max cospa is found! collision_most_prob.globalIndex() = %d , pos.collisionId() = %d , ele.collisionId() = %d, max_cospa = %f", collision_most_prob.globalIndex(), pos.collisionId(), ele.collisionId(), max_cospa);
        vec_cospa.clear();
        vec_cospa.shrink_to_fit();

        // next, check pca between 2 legs in this searching window and select V0s that have the smallest pca to avoid double counting of legs.
        float v0pca = pca_map[std::make_tuple(pos.globalIndex(), ele.globalIndex(), collision_most_prob.globalIndex())];
        bool is_closest_v0 = true;
        for (const auto& [key_tmp, value_tmp] : pca_map) {
          auto pos_tmp = tracks.rawIteratorAt(std::get<0>(key_tmp));
          auto ele_tmp = tracks.rawIteratorAt(std::get<1>(key_tmp));

          float v0pca_tmp = value_tmp;
          // float v0pca_tmp = 999.f;
          // if(pca_map.find(std::make_tuple(pos_tmp.globalIndex(), ele_tmp.globalIndex(), collision_most_prob.globalIndex())) != pca_map.end()){
          //   v0pca_tmp = pca_map[std::make_tuple(pos_tmp.globalIndex(), ele_tmp.globalIndex(), collision_most_prob.globalIndex())];
          // }

          if (ele.globalIndex() == ele_tmp.globalIndex() && pos.globalIndex() == pos_tmp.globalIndex()) { // skip exactly the same V0
            continue;
          }
          if ((ele.globalIndex() == ele_tmp.globalIndex() || pos.globalIndex() == pos_tmp.globalIndex()) && v0pca > v0pca_tmp) {
            // LOGF(info, "!reject! | collision id = %d | posid1 = %d , eleid1 = %d , posid2 = %d , eleid2 = %d , pca1 = %f , pca2 = %f",
            // collision.globalIndex(), pos.globalIndex(), ele.globalIndex(), pos_tmp.globalIndex(), ele_tmp.globalIndex(), v0pca, v0pca_tmp);
            is_closest_v0 = false;
            break;
          }
        } // end of pca_map loop

        if (is_closest_v0 && used_pair_map.find(std::make_pair(pos.globalIndex(), ele.globalIndex())) == used_pair_map.end()) {
          // LOGF(info, "store : pos.globalIndex() = %d , ele.globalIndex() = %d , collision.globalIndex() = %d , cospa = %f , pca = %f", std::get<0>(key), std::get<1>(key), std::get<2>(key), value, pca_map[key]);
          fillV0Table(collision_most_prob, ele, pos, true);
          used_pair_map[std::make_pair(pos.globalIndex(), ele.globalIndex())] = true;
        }
      } // end of pca_map loop
      used_pair_map.clear();

      pca_map.clear();
      cospa_map.clear();

      negTracks_sw.clear();
      posTracks_sw.clear();
      negTracks_sw.shrink_to_fit();
      posTracks_sw.shrink_to_fit();
    } // end of collision loop

  } // end of process
  PROCESS_SWITCH(createPCM, processSA, "create V0s with stand-alone way", true);

  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;
  void processTrkCollAsso(aod::TrackAssoc const& trackIndices, FullTracksExtIU const&, aod::Collisions const& collisions, aod::BCsWithTimestamps const&)
  {
    for (const auto& collision : collisions) {
      registry.fill(HIST("hEventCounter"), 1);

      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, collision.globalIndex());

      // LOGF(info,"%d tracks in collision %d", trackIdsThisCollision.size(), collision.globalIndex());
      for (const auto& [eleId, posId] : combinations(CombinationsStrictlyUpperIndexPolicy(trackIdsThisCollision, trackIdsThisCollision))) {
        auto ele = eleId.track_as<FullTracksExtIU>();
        auto pos = posId.track_as<FullTracksExtIU>();
        // LOGF(info,"eleId = %d , posId = %d", ele.globalIndex(), pos.globalIndex());

        if (ele.sign() * pos.sign() > 0) { // reject same sign combination
          continue;
        }

        if (!isSelected(ele) || !isSelected(pos)) {
          continue;
        }

        if (ele.sign() < 0) {
          fillV0Table(collision, ele, pos, true);
        } else {
          fillV0Table(collision, pos, ele, true);
        }
      }
    } // end of collision loop
  } // end of process
  PROCESS_SWITCH(createPCM, processTrkCollAsso, "create V0s with track-to-collision associator", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<createPCM>(cfgc, TaskName{"v0-finder"})};
}
