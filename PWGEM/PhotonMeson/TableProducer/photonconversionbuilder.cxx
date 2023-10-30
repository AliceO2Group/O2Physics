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
#include "DCAFitter/DCAFitterN.h"
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

using MyTracksIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCFullEl, aod::pidTPCFullPi>;

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
  Configurable<bool> d_UseAbsDCA{"d_UseAbsDCA", true, "Use Abs DCAs"};
  Configurable<bool> d_UseWeightedPCA{"d_UseWeightedPCA", false, "Vertices use cov matrices"};
  Configurable<int> useMatCorrType{"useMatCorrType", 0, "0: none, 1: TGeo, 2: LUT"};

  Configurable<float> dcanegtopv{"dcanegtopv", 0.1, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", 0.1, "DCA Pos To PV"};
  Configurable<float> min_v0cospa{"min_v0cospa", 0.99, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<float> max_dcav0dau{"max_dcav0dau", 1.5, "max distance btween 2 legs"};
  Configurable<float> min_v0radius{"min_v0radius", 1.0, "min v0 radius"};
  Configurable<float> margin_r{"margin_r", 3.0, "margin for r cut in cm"};
  Configurable<float> margin_z{"margin_z", 7.0, "margin for z cut in cm"};

  Configurable<float> min_pt_leg{"min_pt_leg", 0.05, "min pT for v0 legs at SV"};
  Configurable<float> min_pt_v0{"min_pt_v0", 0.05, "min pT for v0 photons at PV"};
  Configurable<float> max_eta_v0{"max_eta_v0", 0.9, "max eta for v0 photons at PV"};
  Configurable<int> mincrossedrows{"mincrossedrows", 10, "min crossed rows"};
  Configurable<float> maxchi2tpc{"maxchi2tpc", 4.0, "max chi2/NclsTPC"};
  Configurable<float> maxchi2its{"maxchi2its", 5.0, "max chi2/NclsITS"};
  Configurable<float> maxpt_itsonly{"maxpt_itsonly", 0.15, "max pT for ITSonly tracks at SV"};
  Configurable<float> maxTPCNsigmaEl{"maxTPCNsigmaEl", 4.0, "max. TPC n sigma for electron"};
  Configurable<float> kfMassConstrain{"kfMassConstrain", -1.f, "mass constrain for the KFParticle mother particle"};

  int mRunNumber;
  float d_bz;
  float maxSnp;  // max sine phi for propagation
  float maxStep; // max step size (cm) for propagation
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut = nullptr;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  //// Define o2 fitter, 2-prong, active memory (no need to redefine per event)
  // o2::vertexing::DCAFitterN<2> fitter;

  HistogramRegistry registry{
    "registry",
    {
      {"hCollisionCounter", "hCollisionCounter", {HistType::kTH1F, {{1, 0.5f, 1.5f}}}},
      {"V0/hAP", "Armenteros Podolanski;#alpha;q_{T} (GeV/c)", {HistType::kTH2F, {{200, -1.0f, 1.0f}, {250, 0, 0.25}}}},
      {"V0/hConversionPointXY", "conversion point in XY;X (cm);Y (cm)", {HistType::kTH2F, {{400, -100.0f, 100.0f}, {400, -100.f, 100.f}}}},
      {"V0/hConversionPointRZ", "conversion point in RZ;Z (cm);R_{xy} (cm)", {HistType::kTH2F, {{200, -100.0f, 100.0f}, {200, 0.f, 100.f}}}},
      {"V0/hPt", "pT of V0 at PV;p_{T,#gamma} (GeV/c)", {HistType::kTH1F, {{1000, 0.0f, 10.0f}}}},
      {"V0/hEtaPhi", "#eta vs. #varphi of V0 at PV;#varphi (rad.);#eta", {HistType::kTH2F, {{72, 0.0f, TMath::TwoPi()}, {400, -2, +2}}}},
      {"V0/hCosPA", "cosine of pointing angle;cosine of pointing angle", {HistType::kTH1F, {{100, 0.9f, 1.f}}}},
      {"V0/hPCA", "distance between 2 legs at SV;PCA (cm)", {HistType::kTH1F, {{500, 0.0f, 5.f}}}},
      {"V0/hMee_SVPV", "mee at PV and SV;m_{ee} at PV (GeV/c^{2});m_{ee} at SV (GeV/c^{2})", {HistType::kTH2F, {{100, 0.0f, 0.1f}, {100, 0, 0.1f}}}},
      {"V0/hMeeSV_Rxy", "mee at SV vs. R_{xy};R_{xy} (cm);m_{ee} at SV (GeV/c^{2})", {HistType::kTH2F, {{200, 0.0f, 100.f}, {100, 0, 0.1f}}}},
      {"V0/hKFChi2", "V0 KF chi2;p_{T,#gamma} (GeV/c);KF #chi^{2}/ndf", {HistType::kTH2F, {{100, 0.0f, 10.f}, {100, 0, 100}}}},
      {"V0/hKPtDiff", "V0 KF pt leg sum vs. gamma pt;p_{T,#gamma} (GeV/c);p_{T,ee} (GeV/c)", {HistType::kTH2F, {{1000, 0.0f, 10.f}, {1000, 0, 10}}}},
      {"V0Leg/hPt", "pT of leg at SV;p_{T,e} (GeV/c)", {HistType::kTH1F, {{1000, 0.0f, 10.0f}}}},
      {"V0Leg/hEtaPhi", "#eta vs. #varphi of leg at SV;#varphi (rad.);#eta", {HistType::kTH2F, {{72, 0.0f, TMath::TwoPi()}, {400, -2, +2}}}},
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

    //// initialize O2 2-prong fitter (only once)
    // fitter.setPropagateToPCA(true);
    // fitter.setMaxR(200.);
    // fitter.setMinParamChange(1e-3);
    // fitter.setMinRelChi2Change(0.9);
    // fitter.setMaxDZIni(1e9);
    // fitter.setMaxChi2(1e9);
    // fitter.setUseAbsDCA(d_UseAbsDCA);
    // fitter.setWeightedFinalPCA(d_UseWeightedPCA);

    if (useMatCorrType == 1)
      matCorr = o2::base::Propagator::MatCorrType::USEMatCorrTGeo;
    if (useMatCorrType == 2)
      matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
    // fitter.setMatCorrType(matCorr);
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    // In case override, don't proceed, please - no CCDB access required
    if (d_bz_input > -990) {
      d_bz = d_bz_input;
      // fitter.setBz(d_bz);
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
    // Set magnetic field value once known
    // fitter.setBz(d_bz);

    if (useMatCorrType == 2) {
      // setMatLUT only after magfield has been initalized (setMatLUT has implicit and problematic init field call if not)
      o2::base::Propagator::Instance()->setMatLUT(lut);
    }
    /// Set magnetic field for KF vertexing
    float magneticField = o2::base::Propagator::Instance()->getNominalBz();
    KFParticle::SetField(magneticField);
  }

  //  static float v0_alpha(float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg)
  //  {
  //    float momTot = RecoDecay::p(pxpos + pxneg, pypos + pyneg, pzpos + pzneg);
  //    float lQlNeg = RecoDecay::dotProd(std::array{pxneg, pyneg, pzneg}, std::array{pxpos + pxneg, pypos + pyneg, pzpos + pzneg}) / momTot;
  //    float lQlPos = RecoDecay::dotProd(std::array{pxpos, pypos, pzpos}, std::array{pxpos + pxneg, pypos + pyneg, pzpos + pzneg}) / momTot;
  //    return (lQlPos - lQlNeg) / (lQlPos + lQlNeg); // longitudinal momentum asymmetry of v0
  //  }
  //
  //  static float v0_qt(float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg)
  //  {
  //    float momTot = RecoDecay::p2(pxpos + pxneg, pypos + pyneg, pzpos + pzneg);
  //    float dp = RecoDecay::dotProd(std::array{pxneg, pyneg, pzneg}, std::array{pxpos + pxneg, pypos + pyneg, pzpos + pzneg});
  //    return std::sqrt(RecoDecay::p2(pxneg, pyneg, pzneg) - dp * dp / momTot); // qt of v0
  //  }

  template <typename TTrack>
  bool checkV0leg(TTrack const& track)
  {
    if (!track.hasITS() && !track.hasTPC()) {
      return false;
    }

    if (track.hasTPC()) {
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
    }

    return true;
  }

  template <typename TTrack, typename TKFParticle>
  void fillTrackTable(TTrack const& track, TKFParticle const& kfp, float dcaXY, float dcaZ)
  {
    v0legs(track.collisionId(), track.globalIndex(), track.sign(),
           kfp.GetPx(), kfp.GetPy(), kfp.GetPz(), dcaXY, dcaZ,
           track.tpcNClsFindable(), track.tpcNClsFindableMinusFound(), track.tpcNClsFindableMinusCrossedRows(),
           track.tpcChi2NCl(), track.tpcInnerParam(), track.tpcSignal(),
           track.tpcNSigmaEl(), track.tpcNSigmaPi(),
           track.itsClusterMap(), track.itsChi2NCl(), track.detectorMap(),
           track.x(), track.y(), track.z(), track.tgl(), track.signed1Pt());
  }

  template <class TCollision, class TTrack, typename TV0>
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

    if (!checkV0leg(pos) || !checkV0leg(ele)) {
      return;
    }
    // LOGF(info, "v0.collisionId() = %d , v0.posTrackId() = %d , v0.negTrackId() = %d", v0.collisionId(), v0.posTrackId(), v0.negTrackId());

    // Calculate DCA with respect to the collision associated to the v0, not individual tracks
    gpu::gpustd::array<float, 2> dcaInfo;

    auto pTrack = getTrackParCov(pos);
    pTrack.setPID(o2::track::PID::Electron);
    o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, pTrack, 2.f, matCorr, &dcaInfo);
    auto posdcaXY = dcaInfo[0];
    auto posdcaZ = dcaInfo[1];

    auto nTrack = getTrackParCov(ele);
    nTrack.setPID(o2::track::PID::Electron);
    o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, nTrack, 2.f, matCorr, &dcaInfo);
    auto eledcaXY = dcaInfo[0];
    auto eledcaZ = dcaInfo[1];

    if (fabs(posdcaXY) < dcapostopv || fabs(eledcaXY) < dcanegtopv) {
      return;
    }

    float xyz[3] = {0.f, 0.f, 0.f};
    Vtx_recalculation(o2::base::Propagator::Instance(), pos, ele, xyz, matCorr);

    float rxy = RecoDecay::sqrtSumOfSquares(xyz[0], xyz[1]);
    if (rxy > std::min(pos.x(), ele.x()) + margin_r || rxy < min_v0radius) {
      return;
    }

    if (rxy < abs(xyz[2]) * TMath::Tan(2 * TMath::ATan(TMath::Exp(-max_eta_v0))) - margin_z) {
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
    if (cospa_kf < min_v0cospa) {
      return;
    }

    //// Apply a topological constraint of the gamma to the PV. Parameters will be given at the primary vertex.
    // KFParticle gammaKF_PV = gammaKF_DecayVtx;
    // gammaKF_PV.SetProductionVertex(KFPV);

    float v0pt = RecoDecay::sqrtSumOfSquares(gammaKF_DecayVtx.GetPx(), gammaKF_DecayVtx.GetPy());
    float v0eta = RecoDecay::eta(std::array{gammaKF_DecayVtx.GetPx(), gammaKF_DecayVtx.GetPy(), gammaKF_DecayVtx.GetPz()});
    float v0phi = RecoDecay::phi(gammaKF_DecayVtx.GetPx(), gammaKF_DecayVtx.GetPy()) > 0.f ? RecoDecay::phi(gammaKF_DecayVtx.GetPx(), gammaKF_DecayVtx.GetPy()) : RecoDecay::phi(gammaKF_DecayVtx.GetPx(), gammaKF_DecayVtx.GetPy()) + TMath::TwoPi();

    if (fabs(v0eta) > max_eta_v0 || v0pt < min_pt_v0) {
      return;
    }

    float chi2kf = -1.f;
    if (gammaKF_DecayVtx.GetNDF() > 0) {
      chi2kf = gammaKF_DecayVtx.GetChi2() / gammaKF_DecayVtx.GetNDF();
    }

    KFParticle kfp_pos_DecayVtx = kfp_pos;  // Don't set Primary Vertex
    KFParticle kfp_ele_DecayVtx = kfp_ele;  // Don't set Primary Vertex
    kfp_pos_DecayVtx.TransportToPoint(xyz); // Don't set Primary Vertex
    kfp_ele_DecayVtx.TransportToPoint(xyz); // Don't set Primary Vertex

    KFParticle kfp_pos_PV = kfp_pos_DecayVtx;
    KFParticle kfp_ele_PV = kfp_ele_DecayVtx;
    kfp_pos_PV.SetProductionVertex(KFPV);
    kfp_ele_PV.SetProductionVertex(KFPV);

    float pca_kf = kfp_pos_DecayVtx.GetDistanceFromParticle(kfp_ele_DecayVtx);
    if (pca_kf > max_dcav0dau) {
      return;
    }

    float pos_pt = RecoDecay::sqrtSumOfSquares(kfp_pos_DecayVtx.GetPx(), kfp_pos_DecayVtx.GetPy());
    float ele_pt = RecoDecay::sqrtSumOfSquares(kfp_ele_DecayVtx.GetPx(), kfp_ele_DecayVtx.GetPy());
    if (pos_pt < min_pt_leg || ele_pt < min_pt_leg) {
      return;
    }

    if (isITSonlyTrack(pos)) {
      float legpt = RecoDecay::sqrtSumOfSquares(kfp_pos_DecayVtx.GetPx(), kfp_pos_DecayVtx.GetPy());
      if (legpt > maxpt_itsonly) {
        return;
      }
    }

    if (isITSonlyTrack(ele)) {
      float legpt = RecoDecay::sqrtSumOfSquares(kfp_ele_DecayVtx.GetPx(), kfp_ele_DecayVtx.GetPy());
      if (legpt > maxpt_itsonly) {
        return;
      }
    }

    float alpha = v0_alpha(kfp_pos_DecayVtx.GetPx(), kfp_pos_DecayVtx.GetPy(), kfp_pos_DecayVtx.GetPz(), kfp_ele_DecayVtx.GetPx(), kfp_ele_DecayVtx.GetPy(), kfp_ele_DecayVtx.GetPz());
    float qt = v0_qt(kfp_pos_DecayVtx.GetPx(), kfp_pos_DecayVtx.GetPy(), kfp_pos_DecayVtx.GetPz(), kfp_ele_DecayVtx.GetPx(), kfp_ele_DecayVtx.GetPy(), kfp_ele_DecayVtx.GetPz());
    if (!checkAP(alpha, qt, 0.95, 0.01)) { // store only photon conversions
      return;
    }
    pca_map[std::make_tuple(v0.globalIndex(), collision.globalIndex(), pos.globalIndex(), ele.globalIndex())] = pca_kf;
    cospa_map[std::make_tuple(v0.globalIndex(), collision.globalIndex(), pos.globalIndex(), ele.globalIndex())] = cospa_kf;

    if (filltable) {
      registry.fill(HIST("V0/hAP"), alpha, qt);
      registry.fill(HIST("V0/hConversionPointXY"), xyz[0], xyz[1]);
      registry.fill(HIST("V0/hConversionPointRZ"), xyz[2], rxy);
      registry.fill(HIST("V0/hPt"), v0pt);
      registry.fill(HIST("V0/hEtaPhi"), v0phi, v0eta);
      registry.fill(HIST("V0/hCosPA"), cospa_kf);
      registry.fill(HIST("V0/hPCA"), pca_kf);
      registry.fill(HIST("V0/hKFChi2"), v0pt, chi2kf);

      float v0pt_sv = RecoDecay::sqrtSumOfSquares(gammaKF_DecayVtx.GetPx(), gammaKF_DecayVtx.GetPy());
      float eept_sv = RecoDecay::sqrtSumOfSquares(kfp_pos_DecayVtx.GetPx() + kfp_ele_DecayVtx.GetPx(), kfp_pos_DecayVtx.GetPy() + kfp_ele_DecayVtx.GetPy());
      registry.fill(HIST("V0/hKPtDiff"), v0pt_sv, eept_sv);

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

      ROOT::Math::PxPyPzMVector vpos_pv(kfp_pos_PV.GetPx(), kfp_pos_PV.GetPy(), kfp_pos_PV.GetPz(), o2::constants::physics::MassElectron);
      ROOT::Math::PxPyPzMVector vele_pv(kfp_ele_PV.GetPx(), kfp_ele_PV.GetPy(), kfp_ele_PV.GetPz(), o2::constants::physics::MassElectron);
      ROOT::Math::PxPyPzMVector v0_pv = vpos_pv + vele_pv;

      ROOT::Math::PxPyPzMVector vpos_sv(kfp_pos_DecayVtx.GetPx(), kfp_pos_DecayVtx.GetPy(), kfp_pos_DecayVtx.GetPz(), o2::constants::physics::MassElectron);
      ROOT::Math::PxPyPzMVector vele_sv(kfp_ele_DecayVtx.GetPx(), kfp_ele_DecayVtx.GetPy(), kfp_ele_DecayVtx.GetPz(), o2::constants::physics::MassElectron);
      ROOT::Math::PxPyPzMVector v0_sv = vpos_sv + vele_sv;
      registry.fill(HIST("V0/hMee_SVPV"), v0_pv.M(), v0_sv.M());
      registry.fill(HIST("V0/hMeeSV_Rxy"), rxy, v0_sv.M());

      v0photonskf(collision.globalIndex(), v0legs.lastIndex() + 1, v0legs.lastIndex() + 2,
                  gammaKF_DecayVtx.GetX(), gammaKF_DecayVtx.GetY(), gammaKF_DecayVtx.GetZ(),
                  gammaKF_DecayVtx.GetPx(), gammaKF_DecayVtx.GetPy(), gammaKF_DecayVtx.GetPz(),
                  v0_sv.M(),
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

  void process(aod::Collisions const& collisions, aod::V0s const& v0s, MyTracksIU const& tracks, aod::BCsWithTimestamps const&)
  {
    for (auto& collision : collisions) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      registry.fill(HIST("hCollisionCounter"), 1);

      auto v0s_per_coll = v0s.sliceBy(perCollision, collision.globalIndex());
      // LOGF(info, "n v0 = %d", v0s_per_coll.size());
      for (auto& v0 : v0s_per_coll) {
        // LOGF(info, "collision.globalIndex() = %d, v0.globalIndex() = %d, v0.posTrackId() = %d, v0.negTrackId() = %d", collision.globalIndex(), v0.globalIndex(), v0.posTrackId() , v0.negTrackId());
        fillV0Table<aod::Collisions, MyTracksIU>(v0, false);
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
        fillV0Table<aod::Collisions, MyTracksIU>(v0, true);
        stored_v0Ids.emplace_back(std::make_pair(posId, eleId));
      }
    } // end of pca_map loop
    // LOGF(info, "pca_map.size() = %d", pca_map.size());
    pca_map.clear();
    cospa_map.clear();
    stored_v0Ids.clear();
    stored_v0Ids.shrink_to_fit();
  } // end of process
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<PhotonConversionBuilder>(cfgc, TaskName{"photon-conversion-builder"})};
}
