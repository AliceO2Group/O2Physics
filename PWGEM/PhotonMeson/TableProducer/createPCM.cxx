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
// ========================
//
// This code produces photon data tables.
//    Please write to: daiki.sekihata@cern.ch

#include <array>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/CollisionAssociation.h"
#include "DCAFitter/DCAFitterN.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using std::array;

using FullTracksExtIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA>;
// using FullTracksExt = soa::Join<aod::Tracks, aod::TracksCov, aod::TracksExtra, aod::TracksDCA>;
using FullTrackExtIU = FullTracksExtIU::iterator;

struct createPCM {
  SliceCache cache;
  Preslice<aod::TracksIU> perCol = o2::aod::track::collisionId;
  Produces<aod::StoredV0Datas> v0data;

  // Basic checks
  HistogramRegistry registry{
    "createPCM",
    {
      {"hEventCounter", "hEventCounter", {HistType::kTH1F, {{5, 0.5f, 5.5f}}}},
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

  Configurable<float> minv0cospa{"minv0cospa", 0.95, "minimum V0 CosPA"};
  Configurable<float> maxdcav0dau{"maxdcav0dau", 1.5, "max DCA between V0 Daughters"};
  Configurable<float> v0Rmin{"v0Rmin", 0.0, "v0Rmin"};
  Configurable<float> v0Rmax{"v0Rmax", 180.0, "v0Rmax"};
  Configurable<float> dcamin{"dcamin", 0.1, "dcamin"};
  Configurable<float> dcamax{"dcamax", 1e+10, "dcamax"};
  Configurable<float> minpt{"minpt", 0.01, "min pT for single track in GeV/c"};
  Configurable<float> maxeta{"maxeta", 0.9, "eta acceptance for single track"};
  Configurable<int> mincrossedrows{"mincrossedrows", 10, "min crossed rows"};
  Configurable<float> maxchi2tpc{"maxchi2tpc", 4.0, "max chi2/NclsTPC"};
  Configurable<float> min_tpcdEdx{"min_tpcdEdx", 30.0, "min TPC dE/dx"};
  Configurable<float> max_tpcdEdx{"max_tpcdEdx", 110.0, "max TPC dE/dx"};
  Configurable<bool> useTPConly{"useTPConly", false, "Use truly TPC only tracks for V0 finder"};
  Configurable<bool> rejectTPConly{"rejectTPConly", false, "Reject truly TPC only tracks for V0 finder"};

  int mRunNumber;
  float d_bz;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut = nullptr;
  o2::vertexing::DCAFitterN<2> fitter;

  void init(InitContext& context)
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

    // Material correction in the DCA fitter
    o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
    if (useMatCorrType == 1)
      matCorr = o2::base::Propagator::MatCorrType::USEMatCorrTGeo;
    if (useMatCorrType == 2)
      matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
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
      if (fabs(d_bz) > 1e-5) {
        grpmag.setL3Current(30000.f / (d_bz / 5.0f));
      }
      o2::base::Propagator::initFieldFromGRP(&grpmag);
      mRunNumber = bc.runNumber();
      return;
    }

    auto run3grp_timestamp = bc.timestamp();
    o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      // Fetch magnetic field from ccdb for current collision
      d_bz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (!grpmag) {
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
  bool IsTPConlyTrack(TTrack const& track)
  {
    if (track.hasTPC() && (!track.hasITS() && !track.hasTOF() && !track.hasTRD())) {
      return true;
    } else {
      return false;
    }
  }

  template <typename TCollision, typename TTrack>
  void fillV0Table(TCollision const& collision, TTrack const& ele, TTrack const& pos)
  {
    array<float, 3> pVtx = {collision.posX(), collision.posY(), collision.posZ()};
    array<float, 3> svpos = {0.}; // secondary vertex position
    array<float, 3> pvec0 = {0.};
    array<float, 3> pvec1 = {0.};

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

    float px = pvec0[0] + pvec1[0];
    float py = pvec0[1] + pvec1[1];
    float pz = pvec0[2] + pvec1[2];

    float v0dca = fitter.getChi2AtPCACandidate(); // distance between 2 legs.
    float v0CosinePA = RecoDecay::cpa(pVtx, array{svpos[0], svpos[1], svpos[2]}, array{px, py, pz});
    float v0radius = RecoDecay::sqrtSumOfSquares(svpos[0], svpos[1]);

    if (v0dca > maxdcav0dau) {
      return;
    }
    if (v0radius < v0Rmin || v0Rmax < v0radius) {
      return;
    }
    if (v0CosinePA < minv0cospa) {
      return;
    }

    v0data(pos.globalIndex(), ele.globalIndex(), collision.globalIndex(), -1,
           fitter.getTrack(0).getX(), fitter.getTrack(1).getX(),
           svpos[0], svpos[1], svpos[2],
           pvec0[0], pvec0[1], pvec0[2],
           pvec1[0], pvec1[1], pvec1[2],
           v0dca, pos.dcaXY(), ele.dcaXY());
  }

  Filter trackFilter = o2::aod::track::pt > minpt&& nabs(o2::aod::track::eta) < maxeta&& dcamin < nabs(o2::aod::track::dcaXY) && nabs(o2::aod::track::dcaXY) < dcamax&& o2::aod::track::tpcChi2NCl < maxchi2tpc&& min_tpcdEdx < o2::aod::track::tpcSignal&& o2::aod::track::tpcSignal < max_tpcdEdx;
  using MyFilteredTracks = soa::Filtered<FullTracksExtIU>;
  Partition<MyFilteredTracks> posTracks = o2::aod::track::signed1Pt > 0.f;
  Partition<MyFilteredTracks> negTracks = o2::aod::track::signed1Pt < 0.f;

  void processSA(MyFilteredTracks const& tracks, aod::Collisions const& collisions, aod::BCsWithTimestamps const&)
  {
    for (auto& collision : collisions) {
      registry.fill(HIST("hEventCounter"), 1);

      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      auto negTracks_coll = negTracks->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);
      auto posTracks_coll = posTracks->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);

      // LOGF(info, "collision.globalIndex() = %d , negTracks_coll.size() = %d , posTracks_coll.size() = %d", collision.globalIndex(), negTracks_coll.size(), posTracks_coll.size());

      for (auto& [ele, pos] : combinations(CombinationsFullIndexPolicy(negTracks_coll, posTracks_coll))) {
        if (!ele.hasTPC() || !pos.hasTPC()) {
          continue;
        }
        if (ele.tpcNClsCrossedRows() < mincrossedrows || pos.tpcNClsCrossedRows() < mincrossedrows) {
          continue;
        }

        if (useTPConly && (!IsTPConlyTrack(ele) || !IsTPConlyTrack(pos))) {
          continue;
        }
        if (rejectTPConly && (IsTPConlyTrack(ele) || IsTPConlyTrack(pos))) {
          continue;
        }
        fillV0Table(collision, ele, pos);
      }
    } // end of collision loop
  }   // end of process
  PROCESS_SWITCH(createPCM, processSA, "create V0s with stand-alone way", true);

  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;
  void processTrkCollAsso(aod::TrackAssoc const& trackIndices, FullTracksExtIU const& tracks, aod::Collisions const& collisions, aod::BCsWithTimestamps const&)
  {
    for (auto& collision : collisions) {
      registry.fill(HIST("hEventCounter"), 1);

      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, collision.globalIndex());

      // LOGF(info,"%d tracks in collision %d", trackIdsThisCollision.size(), collision.globalIndex());
      for (auto& [eleId, posId] : combinations(CombinationsStrictlyUpperIndexPolicy(trackIdsThisCollision, trackIdsThisCollision))) {
        auto ele = eleId.track_as<FullTracksExtIU>();
        auto pos = posId.track_as<FullTracksExtIU>();
        // LOGF(info,"eleId = %d , posId = %d", ele.globalIndex(), pos.globalIndex());

        if (ele.sign() * pos.sign() > 0) { // reject same sign combination
          continue;
        }
        if ((abs(ele.dcaXY()) < dcamin || dcamax < abs(ele.dcaXY())) || (abs(pos.dcaXY()) < dcamin || dcamax < abs(pos.dcaXY()))) {
          continue;
        }
        if (!ele.hasTPC() || !pos.hasTPC()) {
          continue;
        }
        if (ele.tpcNClsCrossedRows() < mincrossedrows || pos.tpcNClsCrossedRows() < mincrossedrows) {
          continue;
        }
        if (ele.tpcChi2NCl() > maxchi2tpc || pos.tpcChi2NCl() > maxchi2tpc) {
          continue;
        }
        if (abs(ele.eta()) > maxeta || abs(pos.eta()) > maxeta) {
          continue;
        }
        if (ele.pt() < minpt || pos.pt() < minpt) {
          continue;
        }
        if (ele.tpcSignal() < min_tpcdEdx || max_tpcdEdx < ele.tpcSignal()) {
          continue;
        }
        if (pos.tpcSignal() < min_tpcdEdx || max_tpcdEdx < pos.tpcSignal()) {
          continue;
        }

        if (useTPConly && (!IsTPConlyTrack(ele) || !IsTPConlyTrack(pos))) {
          continue;
        }
        if (rejectTPConly && (IsTPConlyTrack(ele) || IsTPConlyTrack(pos))) {
          continue;
        }

        if (ele.sign() < 0) {
          fillV0Table(collision, ele, pos);
        } else {
          fillV0Table(collision, pos, ele);
        }
      }
    } // end of collision loop
  }   // end of process
  PROCESS_SWITCH(createPCM, processTrkCollAsso, "create V0s with track-to-collision associator", false);
};

// Extends the v0data table with expression columns
struct v0Initializer {
  Spawns<aod::V0Datas> v0datas;
  void init(InitContext const&) {}
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<createPCM>(cfgc, TaskName{"v0-finder"}),
    adaptAnalysisTask<v0Initializer>(cfgc, TaskName{"v0-initializer"})};
}
