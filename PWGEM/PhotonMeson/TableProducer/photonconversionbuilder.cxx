// Copyright 2019-2023 CERN and copyright holders of ALICE O2.
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
// \file photonconversionbuilder.cxx
// \brief this task produces photon data table with KFParticle.
//
// \author Daiki Sekihata <daiki.sekihata@cern.ch>, Tokyo

#include <cmath>
#include <array>
#include <cstdlib>
#include <map>
#include <iterator>
#include <utility>

#include "Math/Vector4D.h"

#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"

#include "Tools/KFparticle/KFUtilities.h"

#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/PCMUtilities.h"
#include "PWGEM/PhotonMeson/Utils/TrackSelection.h"

using namespace o2;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace o2::pwgem::photonmeson;
using std::array;

using MyCollisions = aod::Collisions;
using MyCollisionsMC = soa::Join<MyCollisions, aod::McCollisionLabels>;

using MyTracksIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCFullEl, aod::pidTPCFullPi>;
using MyTracksIUMC = soa::Join<MyTracksIU, aod::McTrackLabels>;

struct PhotonConversionBuilder {
  Produces<aod::V0PhotonsKF> v0photonskf;
  Produces<aod::V0Legs> v0legs;
  Produces<aod::V0Recalculation> fFuncTableV0Recalculated;

  // CCDB options
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};

  // Operation and minimisation criteria
  Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};
  Configurable<int> useMatCorrType{"useMatCorrType", 0, "0: none, 1: TGeo, 2: LUT"};

  // single track cuts
  Configurable<int> min_ncluster_tpc{"min_ncluster_tpc", 10, "min ncluster tpc"};
  Configurable<int> mincrossedrows{"mincrossedrows", 10, "min crossed rows"};
  Configurable<float> maxchi2tpc{"maxchi2tpc", 5.0, "max chi2/NclsTPC"}; // default 4.0 + 1.0
  Configurable<float> maxchi2its{"maxchi2its", 6.0, "max chi2/NclsITS"}; // default 5.0 + 1.0
  Configurable<float> maxpt_itsonly{"maxpt_itsonly", 0.15, "max pT for ITSonly tracks at SV"};
  Configurable<float> maxTPCNsigmaEl{"maxTPCNsigmaEl", 4.0, "max. TPC n sigma for electron"};
  Configurable<float> dcanegtopv{"dcanegtopv", 0.1, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", 0.1, "DCA Pos To PV"};
  Configurable<float> min_pt_leg{"min_pt_leg", 0.04, "min pT for v0 legs at SV"};
  Configurable<float> max_mean_its_cluster_size{"max_mean_its_cluster_size", 16.f, "max. <ITS cluster size> x cos(lambda) for ITSonly tracks"}; // this is to suppress random combination for V0s with ITSonly tracks. default 3 + 1 for skimming.
  Configurable<float> maxX{"maxX", 83.1, "max X for track IU"};

  // v0 cuts
  Configurable<float> min_v0cospa_tpconly{"min_v0cospa_tpconly", 0.99, "min V0 CosPA to V0s with TPConly tracks"}; // double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<float> min_v0cospa_its{"min_v0cospa_its", 0.99, "min V0 CosPA to V0s with ITs hits"};               // double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<float> max_dcav0dau_tpconly{"max_dcav0dau_tpconly", 3.0, "max distance btween 2 legs to V0s with TPConly tracks"};
  Configurable<float> max_dcav0dau_its{"max_dcav0dau_its", 0.5, "max distance btween 2 legs to V0s with ITS hits"};
  Configurable<float> max_dcav0dau_itsibss{"max_dcav0dau_itsibss", 1.0, "max distance btween 2 legs to V0s with ITS hits on ITSib SS"};
  Configurable<float> max_dcav0dau_tpc_inner_fc{"max_dcav0dau_tpc_inner_fc", 1.5, "max distance btween 2 legs to V0s with ITS hits on TPC inner FC"};
  Configurable<float> min_v0radius{"min_v0radius", 1.0, "min v0 radius"};
  Configurable<float> margin_r_its{"margin_r_its", 3.0, "margin for r cut in cm"};
  Configurable<float> margin_r_tpconly{"margin_r_tpconly", 7.0, "margin for r cut in cm"};
  Configurable<float> margin_z{"margin_z", 7.0, "margin for z cut in cm"};
  Configurable<float> max_alpha_ap{"max_alpha_ap", 0.95, "max alpha for AP cut"};
  Configurable<float> max_qt_ap{"max_qt_ap", 0.01, "max qT for AP cut"};
  Configurable<float> min_pt_v0{"min_pt_v0", 0.05, "min pT for v0 photons at PV"};
  Configurable<float> max_pt_v0_itsonly{"max_pt_v0_itsonly", 0.3, "max pT for v0 photons wth 2 ITSonly tracks at PV"};
  Configurable<float> max_eta_v0{"max_eta_v0", 0.9, "max eta for v0 photons at PV"};
  Configurable<float> kfMassConstrain{"kfMassConstrain", -1.f, "mass constrain for the KFParticle mother particle"};
  Configurable<float> max_r_req_its{"max_r_req_its", 16.0, "max Rxy for V0 with ITS hits"};
  Configurable<float> min_r_tpconly{"min_r_tpconly", 36.0, "min Rxy for V0 with TPConly tracks"};
  Configurable<float> max_r_itsmft_ss{"max_r_itsmft_ss", 66.0, "max Rxy for ITS/MFT SS"};
  Configurable<float> max_dcatopv_xy_v0{"max_dcatopv_xy_v0", +1e+10, "max. DCAxy to PV for V0"};
  Configurable<float> max_dcatopv_z_v0{"max_dcatopv_z_v0", +1e+10, "max. DCAz to PV for V0"};

  int mRunNumber;
  float d_bz;
  float maxSnp;  // max sine phi for propagation
  float maxStep; // max step size (cm) for propagation
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut = nullptr;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  HistogramRegistry registry{
    "registry",
    {
      {"hCollisionCounter", "hCollisionCounter", {HistType::kTH1F, {{1, 0.5f, 1.5f}}}},
      {"V0/hAP", "Armenteros Podolanski;#alpha;q_{T} (GeV/c)", {HistType::kTH2F, {{200, -1.0f, 1.0f}, {250, 0, 0.25}}}},
      {"V0/hConversionPointXY", "conversion point in XY;X (cm);Y (cm)", {HistType::kTH2F, {{400, -100.0f, 100.0f}, {400, -100.f, 100.f}}}},
      {"V0/hConversionPointRZ", "conversion point in RZ;Z (cm);R_{xy} (cm)", {HistType::kTH2F, {{200, -100.0f, 100.0f}, {200, 0.f, 100.f}}}},
      {"V0/hPt", "pT of V0 at PV;p_{T,#gamma} (GeV/c)", {HistType::kTH1F, {{1000, 0.0f, 10.0f}}}},
      {"V0/hEtaPhi", "#eta vs. #varphi of V0 at PV;#varphi (rad.);#eta", {HistType::kTH2F, {{72, 0.0f, 2 * M_PI}, {400, -2, +2}}}},
      {"V0/hCosPA", "cosine of pointing angle;cosine of pointing angle", {HistType::kTH1F, {{100, 0.9f, 1.f}}}},
      {"V0/hCosPA_Rxy", "cosine of pointing angle;r_{xy} (cm);cosine of pointing angle", {HistType::kTH2F, {{200, 0, 100}, {100, 0.9f, 1.f}}}},
      {"V0/hCosPAXY_Rxy", "cosine of pointing angle;r_{xy} (cm);cosine of pointing angle", {HistType::kTH2F, {{200, 0, 100}, {100, 0.9f, 1.f}}}},
      {"V0/hCosPARZ_Rxy", "cosine of pointing angle;r_{xy} (cm);cosine of pointing angle", {HistType::kTH2F, {{200, 0, 100}, {100, 0.9f, 1.f}}}},
      {"V0/hPCA", "distance between 2 legs at SV;PCA (cm)", {HistType::kTH1F, {{500, 0.0f, 5.f}}}},
      {"V0/hPCA_Rxy", "distance between 2 legs at SV;R_{xy} (cm);PCA (cm)", {HistType::kTH2F, {{200, 0, 100}, {500, 0.0f, 5.f}}}},
      {"V0/hDCAxyz", "DCA to PV;DCA_{xy} (cm);DCA_{z} (cm)", {HistType::kTH2F, {{200, -5.f, +5.f}, {200, -5.f, +5.f}}}},
      {"V0/hMeeSV_Rxy", "mee at SV vs. R_{xy};R_{xy} (cm);m_{ee} at SV (GeV/c^{2})", {HistType::kTH2F, {{200, 0.0f, 100.f}, {100, 0, 0.1f}}}},
      {"V0Leg/hPt", "pT of leg at SV;p_{T,e} (GeV/c)", {HistType::kTH1F, {{1000, 0.0f, 10.0f}}}},
      {"V0Leg/hEtaPhi", "#eta vs. #varphi of leg at SV;#varphi (rad.);#eta", {HistType::kTH2F, {{72, 0.0f, 2 * M_PI}, {400, -2, +2}}}},
      {"V0Leg/hDCAxyz", "DCA xy vs. z to PV;DCA_{xy} (cm);DCA_{z} (cm)", {HistType::kTH2F, {{200, -10.f, 10.f}, {200, -10.f, +10.f}}}},
      {"V0Leg/hdEdx_Pin", "TPC dE/dx vs. p_{in};p_{in} (GeV/c);TPC dE/dx", {HistType::kTH2F, {{1000, 0.f, 10.f}, {200, 0.f, 200.f}}}},
      {"V0Leg/hTPCNsigmaEl", "TPC dE/dx vs. p_{in};p_{in} (GeV/c);n #sigma_{e}^{TPC}", {HistType::kTH2F, {{1000, 0.f, 10.f}, {100, -5.f, +5.f}}}},
    }};

  void init(InitContext& context)
  {
    mRunNumber = 0;
    d_bz = 0;
    maxSnp = 0.85f;  // could be changed later
    maxStep = 2.00f; // could be changed later

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

    if (useMatCorrType == 1)
      matCorr = o2::base::Propagator::MatCorrType::USEMatCorrTGeo;
    if (useMatCorrType == 2)
      matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    // In case override, don't proceed, please - no CCDB access required
    if (d_bz_input > -990) {
      d_bz = d_bz_input;
      o2::parameters::GRPMagField grpmag;
      if (fabs(d_bz) > 1e-5) {
        grpmag.setL3Current(30000.f / (d_bz / 5.0f));
      }
      o2::base::Propagator::initFieldFromGRP(&grpmag);
      mRunNumber = bc.runNumber();
      return;
    }

    auto run3grp_timestamp = bc.timestamp();
    o2::parameters::GRPObject* grpo = 0x0;
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (!skipGRPOquery)
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
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

    if (useMatCorrType == 2) {
      // setMatLUT only after magfield has been initalized (setMatLUT has implicit and problematic init field call if not)
      o2::base::Propagator::Instance()->setMatLUT(lut);
    }
    /// Set magnetic field for KF vertexing
    float magneticField = o2::base::Propagator::Instance()->getNominalBz();
    KFParticle::SetField(magneticField);
  }

  std::pair<int8_t, std::set<uint8_t>> its_ib_Requirement = {0, {0, 1, 2}}; // no hit on 3 ITS ib layers.
  template <bool isMC, typename TTrack>
  bool checkV0leg(TTrack const& track)
  {
    if constexpr (isMC) {
      if (!track.has_mcParticle()) {
        return false;
      }
    }

    if (track.x() > maxX) {
      return false;
    }

    if (!track.hasITS() && !track.hasTPC()) {
      return false;
    }

    if (track.hasITS() && !track.hasTPC() && (track.hasTRD() || track.hasTOF())) { // remove unrealistic track. this should not happen.
      return false;
    }

    if (track.hasTPC()) {
      if (track.tpcNClsFound() < min_ncluster_tpc) {
        return false;
      }
      if (track.tpcNClsCrossedRows() < mincrossedrows || track.tpcChi2NCl() > maxchi2tpc) {
        return false;
      }
      if (abs(track.tpcNSigmaEl()) > maxTPCNsigmaEl) {
        return false;
      }
    }

    if (track.hasITS()) {
      if (track.itsChi2NCl() > maxchi2its) {
        return false;
      }

      auto hits_ib = std::count_if(its_ib_Requirement.second.begin(), its_ib_Requirement.second.end(), [&](auto&& requiredLayer) { return track.itsClusterMap() & (1 << requiredLayer); });
      bool its_ob_only = hits_ib <= its_ib_Requirement.first;
      if (!its_ob_only) {
        return false;
      }

      if (isITSonlyTrack(track)) {
        uint32_t itsClusterSizes = track.itsClusterSizes();
        int total_cluster_size = 0, nl = 0;
        for (unsigned int layer = 0; layer < 7; layer++) {
          int cluster_size_per_layer = (itsClusterSizes >> (layer * 4)) & 0xf;
          if (cluster_size_per_layer > 0) {
            nl++;
          }
          total_cluster_size += cluster_size_per_layer;
        }
        if (static_cast<float>(total_cluster_size) / static_cast<float>(nl) * std::cos(std::atan(track.tgl())) > max_mean_its_cluster_size) {
          return false;
        }
      }
    }

    return true;
  }

  float cospaXY_KF(KFParticle kfp, KFParticle PV)
  {
    float lx = kfp.GetX() - PV.GetX(); // flight length X
    float ly = kfp.GetY() - PV.GetY(); // flight length Y

    float px = kfp.GetPx();
    float py = kfp.GetPy();
    float cospaXY = RecoDecay::dotProd(std::array{lx, ly}, std::array{px, py}) / (RecoDecay::sqrtSumOfSquares(lx, ly) * RecoDecay::sqrtSumOfSquares(px, py));
    if (cospaXY < -1.) {
      return -1.;
    } else if (cospaXY > 1.) {
      return 1.;
    }
    return cospaXY;
  }

  float cospaRZ_KF(KFParticle kfp, KFParticle PV)
  {
    float lx = kfp.GetX() - PV.GetX();              // flight length X
    float ly = kfp.GetY() - PV.GetY();              // flight length Y
    float lz = kfp.GetZ() - PV.GetZ();              // flight length Z
    float lt = RecoDecay::sqrtSumOfSquares(lx, ly); // flight length R, i.e. transverse plane.

    float pt = RecoDecay::sqrtSumOfSquares(kfp.GetPx(), kfp.GetPy());
    float pz = kfp.GetPz();

    float cospaRZ = RecoDecay::dotProd(std::array{lt, lz}, std::array{pt, pz}) / (RecoDecay::sqrtSumOfSquares(lt, lz) * RecoDecay::sqrtSumOfSquares(pt, pz));
    if (cospaRZ < -1.) {
      return -1.;
    } else if (cospaRZ > 1.) {
      return 1.;
    }
    return cospaRZ;
  }

  template <typename TTrack, typename TKFParticle>
  void fillTrackTable(TTrack const& track, TKFParticle const& kfp, float dcaXY, float dcaZ)
  {
    v0legs(track.collisionId(), track.globalIndex(), track.sign(),
           kfp.GetPx(), kfp.GetPy(), kfp.GetPz(), dcaXY, dcaZ,
           track.tpcNClsFindable(), track.tpcNClsFindableMinusFound(), track.tpcNClsFindableMinusCrossedRows(),
           track.tpcChi2NCl(), track.tpcInnerParam(), track.tpcSignal(),
           track.tpcNSigmaEl(), track.tpcNSigmaPi(),
           track.itsClusterSizes(), track.itsChi2NCl(), track.detectorMap(),
           track.x(), track.y(), track.z(), track.tgl());
  }

  template <bool isMC, class TCollision, class TTrack, typename TV0>
  void fillV0Table(TV0 const& v0, const bool filltable)
  {
    // Get tracks
    auto pos = v0.template posTrack_as<TTrack>();
    auto ele = v0.template negTrack_as<TTrack>();
    auto collision = v0.template collision_as<TCollision>(); // collision where this v0 belongs to.

    if (pos.sign() * ele.sign() > 0) { // reject same sign pair
      return;
    }

    if (pos.globalIndex() == ele.globalIndex()) {
      return;
    }

    if (isITSonlyTrack(pos) && !ele.hasITS()) {
      return;
    }
    if (isITSonlyTrack(ele) && !pos.hasITS()) {
      return;
    }

    if (!checkV0leg<isMC>(pos) || !checkV0leg<isMC>(ele)) {
      return;
    }
    // LOGF(info, "v0.collisionId() = %d , v0.posTrackId() = %d , v0.negTrackId() = %d", v0.collisionId(), v0.posTrackId(), v0.negTrackId());

    // Calculate DCA with respect to the collision associated to the v0, not individual tracks
    gpu::gpustd::array<float, 2> dcaInfo;

    auto pTrack = getTrackPar(pos);
    o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, pTrack, 2.f, matCorr, &dcaInfo);
    auto posdcaXY = dcaInfo[0];
    auto posdcaZ = dcaInfo[1];

    auto nTrack = getTrackPar(ele);
    o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, nTrack, 2.f, matCorr, &dcaInfo);
    auto eledcaXY = dcaInfo[0];
    auto eledcaZ = dcaInfo[1];

    if (fabs(posdcaXY) < dcapostopv || fabs(eledcaXY) < dcanegtopv) {
      return;
    }

    float xyz[3] = {0.f, 0.f, 0.f};
    Vtx_recalculation(o2::base::Propagator::Instance(), pos, ele, xyz, matCorr);
    float rxy_tmp = RecoDecay::sqrtSumOfSquares(xyz[0], xyz[1]);
    if (rxy_tmp > maxX + margin_r_tpconly) {
      return;
    }
    if (rxy_tmp < abs(xyz[2]) * TMath::Tan(2 * TMath::ATan(TMath::Exp(-max_eta_v0))) - margin_z) {
      return; // RZ line cut
    }

    KFPTrack kfp_track_pos = createKFPTrackFromTrack(pos);
    KFPTrack kfp_track_ele = createKFPTrackFromTrack(ele);
    KFParticle kfp_pos(kfp_track_pos, -11);
    KFParticle kfp_ele(kfp_track_ele, 11);
    const KFParticle* GammaDaughters[2] = {&kfp_pos, &kfp_ele};

    KFParticle gammaKF;
    gammaKF.SetConstructMethod(2);
    gammaKF.Construct(GammaDaughters, 2);
    if (kfMassConstrain > -0.1) {
      gammaKF.SetNonlinearMassConstraint(kfMassConstrain);
    }
    KFPVertex kfpVertex = createKFPVertexFromCollision(collision);
    KFParticle KFPV(kfpVertex);

    // Transport the gamma to the recalculated decay vertex
    KFParticle gammaKF_DecayVtx = gammaKF; // with respect to (0,0,0)
    gammaKF_DecayVtx.TransportToPoint(xyz);

    float cospa_kf = cpaFromKF(gammaKF_DecayVtx, KFPV);
    if (!ele.hasITS() && !pos.hasITS()) {
      if (cospa_kf < min_v0cospa_tpconly) {
        return;
      }
    } else {
      if (cospa_kf < min_v0cospa_its) {
        return;
      }
    }

    float rxy = RecoDecay::sqrtSumOfSquares(gammaKF_DecayVtx.GetX(), gammaKF_DecayVtx.GetY());
    if (rxy < min_v0radius) {
      return;
    }

    if (pos.hasITS() && ele.hasITS()) {
      if (rxy > std::min(pos.x(), ele.x()) + margin_r_its) {
        return;
      }
    } else if (!pos.hasITS() && ele.hasITS()) {
      if (rxy > std::min(83.f, ele.x()) + margin_r_its) {
        return;
      }
    } else if (pos.hasITS() && !ele.hasITS()) {
      if (rxy > std::min(pos.x(), 83.f) + margin_r_its) {
        return;
      }
    } else if (!pos.hasITS() && !ele.hasITS()) {
      if (rxy > std::min(83.f, 83.f) + margin_r_tpconly) {
        return;
      }
    }

    if (rxy > maxX + margin_r_tpconly) {
      return;
    }
    if (rxy < abs(gammaKF_DecayVtx.GetZ()) * TMath::Tan(2 * TMath::ATan(TMath::Exp(-max_eta_v0))) - margin_z) {
      return; // RZ line cut
    }

    if ((!pos.hasITS() || !ele.hasITS()) && rxy < max_r_req_its) { // conversion points smaller than max_r_req_its have to be detected with ITS hits.
      return;
    }

    if ((!pos.hasITS() && !ele.hasITS()) && rxy < min_r_tpconly) { // TPConly tracks can detect conversion points larger than min_r_tpconly.
      return;
    }

    // Apply a topological constraint of the gamma to the PV. Parameters will be given at the primary vertex.
    KFParticle gammaKF_PV = gammaKF;
    gammaKF_PV.SetProductionVertex(KFPV);
    float v0pt = RecoDecay::sqrtSumOfSquares(gammaKF_PV.GetPx(), gammaKF_PV.GetPy());
    float v0eta = RecoDecay::eta(std::array{gammaKF_PV.GetPx(), gammaKF_PV.GetPy(), gammaKF_PV.GetPz()});
    float v0phi = RecoDecay::phi(gammaKF_PV.GetPx(), gammaKF_PV.GetPy()) > 0.f ? RecoDecay::phi(gammaKF_PV.GetPx(), gammaKF_PV.GetPy()) : RecoDecay::phi(gammaKF_PV.GetPx(), gammaKF_PV.GetPy()) + TMath::TwoPi();

    // KFParticle gammaKF_DecayVtx2 = gammaKF;
    // gammaKF_DecayVtx2.SetProductionVertex(KFPV);
    // gammaKF_DecayVtx2.TransportToPoint(xyz);
    // LOGF(info, "gammaKF_PV.GetPx() = %f, gammaKF_DecayVtx.GetPx() = %f, gammaKF_DecayVtx2.GetPx() = %f", gammaKF_PV.GetPx(), gammaKF_DecayVtx.GetPx(), gammaKF_DecayVtx2.GetPx());
    // LOGF(info, "gammaKF_PV.GetPy() = %f, gammaKF_DecayVtx.GetPy() = %f, gammaKF_DecayVtx2.GetPy() = %f", gammaKF_PV.GetPy(), gammaKF_DecayVtx.GetPy(), gammaKF_DecayVtx2.GetPy());
    // LOGF(info, "gammaKF_PV.GetPz() = %f, gammaKF_DecayVtx.GetPz() = %f, gammaKF_DecayVtx2.GetPz() = %f", gammaKF_PV.GetPz(), gammaKF_DecayVtx.GetPz(), gammaKF_DecayVtx2.GetPz());

    if (fabs(v0eta) > max_eta_v0 || v0pt < min_pt_v0) {
      return;
    }

    if (isITSonlyTrack(ele) && isITSonlyTrack(pos) && v0pt > max_pt_v0_itsonly) {
      return;
    }

    KFParticle kfp_pos_DecayVtx = kfp_pos;  // Don't set Primary Vertex
    KFParticle kfp_ele_DecayVtx = kfp_ele;  // Don't set Primary Vertex
    kfp_pos_DecayVtx.TransportToPoint(xyz); // Don't set Primary Vertex
    kfp_ele_DecayVtx.TransportToPoint(xyz); // Don't set Primary Vertex

    float pca_kf = kfp_pos_DecayVtx.GetDistanceFromParticle(kfp_ele_DecayVtx);
    if (!ele.hasITS() && !pos.hasITS()) { // V0s with TPConly-TPConly
      if (max_r_itsmft_ss < rxy && rxy < maxX + margin_r_tpconly) {
        if (pca_kf > max_dcav0dau_tpc_inner_fc) {
          return;
        }
      } else {
        if (pca_kf > max_dcav0dau_tpconly) {
          return;
        }
      }
    } else { // V0s with ITS hits
      if (rxy < max_r_req_its) {
        if (pca_kf > max_dcav0dau_itsibss) {
          return;
        }
      } else {
        if (pca_kf > max_dcav0dau_its) {
          return;
        }
      }
    }

    float pos_pt = RecoDecay::sqrtSumOfSquares(kfp_pos_DecayVtx.GetPx(), kfp_pos_DecayVtx.GetPy());
    float ele_pt = RecoDecay::sqrtSumOfSquares(kfp_ele_DecayVtx.GetPx(), kfp_ele_DecayVtx.GetPy());
    if (pos_pt < min_pt_leg || ele_pt < min_pt_leg) {
      return;
    }

    if (isITSonlyTrack(pos) && pos_pt > maxpt_itsonly) {
      return;
    }

    if (isITSonlyTrack(ele) && ele_pt > maxpt_itsonly) {
      return;
    }

    // calculate DCAxy,z to PV
    float v0mom = RecoDecay::sqrtSumOfSquares(gammaKF_DecayVtx.GetPx(), gammaKF_DecayVtx.GetPy(), gammaKF_DecayVtx.GetPz());
    float length = RecoDecay::sqrtSumOfSquares(gammaKF_DecayVtx.GetX() - collision.posX(), gammaKF_DecayVtx.GetY() - collision.posY(), gammaKF_DecayVtx.GetZ() - collision.posZ());
    float dca_x_v0_to_pv = (gammaKF_DecayVtx.GetX() - gammaKF_DecayVtx.GetPx() * cospa_kf * length / v0mom) - collision.posX();
    float dca_y_v0_to_pv = (gammaKF_DecayVtx.GetY() - gammaKF_DecayVtx.GetPy() * cospa_kf * length / v0mom) - collision.posY();
    float dca_z_v0_to_pv = (gammaKF_DecayVtx.GetZ() - gammaKF_DecayVtx.GetPz() * cospa_kf * length / v0mom) - collision.posZ();
    float sign_tmp = dca_x_v0_to_pv * dca_y_v0_to_pv > 0 ? +1.f : -1.f;
    float dca_xy_v0_to_pv = RecoDecay::sqrtSumOfSquares(dca_x_v0_to_pv, dca_y_v0_to_pv) * sign_tmp;
    if (abs(dca_xy_v0_to_pv) > max_dcatopv_xy_v0 || abs(dca_z_v0_to_pv) > max_dcatopv_z_v0) {
      return;
    }

    float alpha = v0_alpha(kfp_pos_DecayVtx.GetPx(), kfp_pos_DecayVtx.GetPy(), kfp_pos_DecayVtx.GetPz(), kfp_ele_DecayVtx.GetPx(), kfp_ele_DecayVtx.GetPy(), kfp_ele_DecayVtx.GetPz());
    float qt = v0_qt(kfp_pos_DecayVtx.GetPx(), kfp_pos_DecayVtx.GetPy(), kfp_pos_DecayVtx.GetPz(), kfp_ele_DecayVtx.GetPx(), kfp_ele_DecayVtx.GetPy(), kfp_ele_DecayVtx.GetPz());
    if (!checkAP(alpha, qt, max_alpha_ap, max_qt_ap)) { // store only photon conversions
      return;
    }
    pca_map[std::make_tuple(v0.globalIndex(), collision.globalIndex(), pos.globalIndex(), ele.globalIndex())] = pca_kf;
    cospa_map[std::make_tuple(v0.globalIndex(), collision.globalIndex(), pos.globalIndex(), ele.globalIndex())] = cospa_kf;

    if (filltable) {
      registry.fill(HIST("V0/hAP"), alpha, qt);
      registry.fill(HIST("V0/hConversionPointXY"), gammaKF_DecayVtx.GetX(), gammaKF_DecayVtx.GetY());
      registry.fill(HIST("V0/hConversionPointRZ"), gammaKF_DecayVtx.GetZ(), rxy);
      registry.fill(HIST("V0/hPt"), v0pt);
      registry.fill(HIST("V0/hEtaPhi"), v0phi, v0eta);
      registry.fill(HIST("V0/hCosPA"), cospa_kf);
      registry.fill(HIST("V0/hCosPA_Rxy"), rxy, cospa_kf);
      registry.fill(HIST("V0/hPCA"), pca_kf);
      registry.fill(HIST("V0/hPCA_Rxy"), rxy, pca_kf);
      registry.fill(HIST("V0/hDCAxyz"), dca_xy_v0_to_pv, dca_z_v0_to_pv);

      float cospaXY_kf = cospaXY_KF(gammaKF_DecayVtx, KFPV);
      float cospaRZ_kf = cospaRZ_KF(gammaKF_DecayVtx, KFPV);
      // LOGF(info, "cospa_kf = %f, cospaXY_kf = %f, cospaRZ_kf = %f", cospa_kf, cospaXY_kf, cospaRZ_kf);
      registry.fill(HIST("V0/hCosPAXY_Rxy"), rxy, cospaXY_kf);
      registry.fill(HIST("V0/hCosPARZ_Rxy"), rxy, cospaRZ_kf);

      float chi2kf = gammaKF_DecayVtx.GetChi2() / gammaKF_DecayVtx.GetNDF();

      for (auto& leg : {kfp_pos_DecayVtx, kfp_ele_DecayVtx}) {
        float legpt = RecoDecay::sqrtSumOfSquares(leg.GetPx(), leg.GetPy());
        float legeta = RecoDecay::eta(std::array{leg.GetPx(), leg.GetPy(), leg.GetPz()});
        float legphi = RecoDecay::phi(leg.GetPx(), leg.GetPy()) > 0.f ? RecoDecay::phi(leg.GetPx(), leg.GetPy()) : RecoDecay::phi(leg.GetPx(), leg.GetPy()) + TMath::TwoPi();
        registry.fill(HIST("V0Leg/hPt"), legpt);
        registry.fill(HIST("V0Leg/hEtaPhi"), legphi, legeta);
      } // end of leg loop
      for (auto& leg : {pos, ele}) {
        registry.fill(HIST("V0Leg/hdEdx_Pin"), leg.tpcInnerParam(), leg.tpcSignal());
        registry.fill(HIST("V0Leg/hTPCNsigmaEl"), leg.tpcInnerParam(), leg.tpcNSigmaEl());
      } // end of leg loop
      registry.fill(HIST("V0Leg/hDCAxyz"), posdcaXY, posdcaZ);
      registry.fill(HIST("V0Leg/hDCAxyz"), eledcaXY, eledcaZ);

      ROOT::Math::PxPyPzMVector vpos_sv(kfp_pos_DecayVtx.GetPx(), kfp_pos_DecayVtx.GetPy(), kfp_pos_DecayVtx.GetPz(), o2::constants::physics::MassElectron);
      ROOT::Math::PxPyPzMVector vele_sv(kfp_ele_DecayVtx.GetPx(), kfp_ele_DecayVtx.GetPy(), kfp_ele_DecayVtx.GetPz(), o2::constants::physics::MassElectron);
      ROOT::Math::PxPyPzMVector v0_sv = vpos_sv + vele_sv;
      registry.fill(HIST("V0/hMeeSV_Rxy"), rxy, v0_sv.M());

      v0photonskf(collision.globalIndex(), v0legs.lastIndex() + 1, v0legs.lastIndex() + 2,
                  gammaKF_DecayVtx.GetX(), gammaKF_DecayVtx.GetY(), gammaKF_DecayVtx.GetZ(),
                  gammaKF_PV.GetPx(), gammaKF_PV.GetPy(), gammaKF_PV.GetPz(),
                  v0_sv.M(), dca_xy_v0_to_pv, dca_z_v0_to_pv,
                  cospa_kf, pca_kf, alpha, qt, chi2kf);

      fFuncTableV0Recalculated(xyz[0], xyz[1], xyz[2]);
      fillTrackTable(pos, kfp_pos_DecayVtx, posdcaXY, posdcaZ); // positive leg first
      fillTrackTable(ele, kfp_ele_DecayVtx, eledcaXY, eledcaZ); // negative leg second
    }                                                           // end of fill table
  }

  Preslice<aod::V0s> perCollision = o2::aod::v0::collisionId;
  std::map<std::tuple<int64_t, int64_t, int64_t, int64_t>, float> pca_map;   //(v0.globalIndex(), collision.globalIndex(), pos.globalIndex(), ele.globalIndex()) -> pca
  std::map<std::tuple<int64_t, int64_t, int64_t, int64_t>, float> cospa_map; //(v0.globalIndex(), collision.globalIndex(), pos.globalIndex(), ele.globalIndex()) -> cospa
  std::vector<std::pair<int64_t, int64_t>> stored_v0Ids;                     //(pos.globalIndex(), ele.globalIndex())

  template <bool isMC, typename TCollisions, typename TV0s, typename TTracks, typename TBCs>
  void build(TCollisions const& collisions, TV0s const& v0s, TTracks const& tracks, TBCs const&)
  {
    for (auto& collision : collisions) {
      if constexpr (isMC) {
        if (!collision.has_mcCollision()) {
          continue;
        }
      }

      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      registry.fill(HIST("hCollisionCounter"), 1);

      auto v0s_per_coll = v0s.sliceBy(perCollision, collision.globalIndex());
      // LOGF(info, "n v0 = %d", v0s_per_coll.size());
      for (auto& v0 : v0s_per_coll) {
        // LOGF(info, "collision.globalIndex() = %d, v0.globalIndex() = %d, v0.posTrackId() = %d, v0.negTrackId() = %d", collision.globalIndex(), v0.globalIndex(), v0.posTrackId() , v0.negTrackId());
        fillV0Table<isMC, TCollisions, TTracks>(v0, false);
      } // end of v0 loop
    }   // end of collision loop

    stored_v0Ids.reserve(pca_map.size()); // number of photon candidates per DF

    // find minimal pca
    for (const auto& [key, value] : pca_map) {
      auto v0Id = std::get<0>(key);
      auto collisionId = std::get<1>(key);
      auto posId = std::get<2>(key);
      auto eleId = std::get<3>(key);
      float v0pca = value;
      float cospa = cospa_map[key];
      bool is_closest_v0 = true;
      bool is_most_aligned_v0 = true;

      for (const auto& [key_tmp, value_tmp] : pca_map) {
        auto v0Id_tmp = std::get<0>(key_tmp);
        auto collisionId_tmp = std::get<1>(key_tmp);
        auto posId_tmp = std::get<2>(key_tmp);
        auto eleId_tmp = std::get<3>(key_tmp);
        float v0pca_tmp = value_tmp;
        float cospa_tmp = cospa_map[key_tmp];

        if (v0Id == v0Id_tmp) { // skip exactly the same v0
          continue;
        }

        if (collisionId != collisionId_tmp && eleId == eleId_tmp && posId == posId_tmp && cospa < cospa_tmp) { // same ele and pos, but attached to different collision
          // LOGF(info, "!reject! | collision id = %d | posid1 = %d , eleid1 = %d , posid2 = %d , eleid2 = %d , cospa1 = %f , cospa2 = %f", collisionId, posId, eleId, posId_tmp, eleId_tmp, cospa, cospa_tmp);
          is_most_aligned_v0 = false;
          break;
        }

        if ((eleId == eleId_tmp || posId == posId_tmp) && v0pca > v0pca_tmp) {
          // LOGF(info, "!reject! | collision id = %d | posid1 = %d , eleid1 = %d , posid2 = %d , eleid2 = %d , pca1 = %f , pca2 = %f", collisionId, posId, eleId, posId_tmp, eleId_tmp, v0pca, v0pca_tmp);
          is_closest_v0 = false;
          break;
        }
      } // end of pca_map tmp loop

      bool is_stored = std::find(stored_v0Ids.begin(), stored_v0Ids.end(), std::make_pair(posId, eleId)) != stored_v0Ids.end();
      if (is_closest_v0 && is_most_aligned_v0 && !is_stored) {
        auto v0 = v0s.rawIteratorAt(v0Id);
        // auto collision = collisions.rawIteratorAt(collisionId);
        // auto pos = tracks.rawIteratorAt(posId);
        // auto ele = tracks.rawIteratorAt(eleId);
        // LOGF(info, "!accept! | collision id = %d | v0id1 = %d , posid1 = %d , eleid1 = %d , pca1 = %f , cospa = %f", collisionId, v0Id, posId, eleId, v0pca, cospa);
        fillV0Table<isMC, TCollisions, TTracks>(v0, true);
        stored_v0Ids.emplace_back(std::make_pair(posId, eleId));
      }
    } // end of pca_map loop
    // LOGF(info, "pca_map.size() = %d", pca_map.size());
    pca_map.clear();
    cospa_map.clear();
    stored_v0Ids.clear();
    stored_v0Ids.shrink_to_fit();
  } // end of build

  //! type of V0. 0: built solely for cascades (does not pass standard V0 cuts), 1: standard 2, 3: photon-like with TPC-only use. Regular analysis should always use type 1 or 3.
  Filter v0Filter = o2::aod::v0::v0Type > (uint8_t)0;
  using filteredV0s = soa::Filtered<aod::V0s>;

  void processRec(MyCollisions const& collisions, filteredV0s const& v0s, MyTracksIU const& tracks, aod::BCsWithTimestamps const& bcs)
  {
    build<false>(collisions, v0s, tracks, bcs);
  }
  PROCESS_SWITCH(PhotonConversionBuilder, processRec, "process reconstructed info for data", true);

  void processMC(MyCollisionsMC const& collisions, filteredV0s const& v0s, MyTracksIUMC const& tracks, aod::BCsWithTimestamps const& bcs)
  {
    build<true>(collisions, v0s, tracks, bcs);
  }
  PROCESS_SWITCH(PhotonConversionBuilder, processMC, "process reconstructed info for MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<PhotonConversionBuilder>(cfgc, TaskName{"photon-conversion-builder"})};
}
