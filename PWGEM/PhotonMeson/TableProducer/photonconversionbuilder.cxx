// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
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

/// \file photonconversionbuilder.cxx
/// \brief this task produces photon data table with KFParticle.
/// \author Daiki Sekihata <daiki.sekihata@cern.ch>, Tokyo

#ifndef HomogeneousField
#define HomogeneousField // needed for KFParticle::SetField(magneticField);
#endif

#include "PWGEM/PhotonMeson/Core/EmMlResponsePCM.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCandidate.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/PCMUtilities.h"
#include "PWGEM/PhotonMeson/Utils/TrackSelection.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TPCVDriftManager.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Tools/KFparticle/KFUtilities.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DataFormatsParameters/GRPObject.h>
#include <DetectorsBase/GeometryManager.h>
#include <DetectorsBase/MatLayerCylSet.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/PID.h>

#include <Math/Vector4D.h> // IWYU pragma: keep
#include <Math/Vector4Dfwd.h>
#include <TPDGCode.h>

#include <KFPTrack.h>
#include <KFPVertex.h>
#include <KFParticle.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace o2::pwgem::photonmeson;
using std::array;

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::EMEvSels, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>;
using MyCollisionsWithSWT = soa::Join<MyCollisions, aod::EMSWTriggerBitsTMP>;
using MyCollisionsMC = soa::Join<MyCollisions, aod::McCollisionLabels>;

using MyTracksIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCFullEl, aod::pidTPCFullPi>;
using MyTracksIUMC = soa::Join<MyTracksIU, aod::McTrackLabels, aod::mcTPCTuneOnData>;

enum MatCorrType {
  None = 0,
  TGeo = 1,
  LUT = 2
};

struct PhotonConversionBuilder {
  Produces<aod::V0PhotonsKF> v0photonskf;
  Produces<aod::V0Legs> v0legs;
  Produces<aod::V0LegsXYZ> v0legsXYZ;
  Produces<aod::V0LegsDeDxMC> v0legsDeDxMC;
  Produces<aod::V0PhotonsPhiV> v0photonsphiv;
  // Produces<aod::V0PhotonsKFCov> v0photonskfcov;
  // Produces<aod::EMEventsNgPCM> events_ngpcm;

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
  Configurable<int> min_ncluster_tpc{"min_ncluster_tpc", 0, "min ncluster tpc"};
  Configurable<int> mincrossedrows{"mincrossedrows", 40, "min crossed rows"};
  Configurable<bool> moveTPCTracks{"moveTPCTracks", true, "Move TPC-only tracks under the collision assumption"};
  Configurable<bool> disableITSonlyTracks{"disableITSonlyTracks", false, "disable ITSonly tracks in V0 legs"};
  Configurable<bool> disableTPConlyTracks{"disableTPConlyTracks", false, "disable TPConly tracks in V0 legs"};
  Configurable<bool> requireITShit{"requireITShit", false, "require ITS hit to V0 legs"};

  Configurable<float> maxchi2tpc{"maxchi2tpc", 5.0, "max chi2/NclsTPC"}; // default 4.0 + 1.0
  Configurable<float> maxchi2its{"maxchi2its", 6.0, "max chi2/NclsITS"}; // default 5.0 + 1.0
  Configurable<float> maxpt_itsonly{"maxpt_itsonly", 0.15, "max pT for ITSonly tracks at SV"};
  Configurable<float> maxTPCNsigmaEl{"maxTPCNsigmaEl", 4.0, "max. TPC n sigma for electron"};
  Configurable<float> minTPCNsigmaEl{"minTPCNsigmaEl", -4.0, "min. TPC n sigma for electron"};
  Configurable<float> max_frac_shared_clusters_tpc{"max_frac_shared_clusters_tpc", 999.f, "max fraction of shared clusters in TPC"};
  Configurable<float> dcanegtopv{"dcanegtopv", 0.1, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", 0.1, "DCA Pos To PV"};
  Configurable<float> maxX{"maxX", 83.1, "max X for track IU"};
  Configurable<float> min_pt_trackiu{"min_pt_trackiu", 0.05, "min pT for trackiu"}; // this comes from online processing. pT of track seed is above 50 MeV/c in B = 0.5 T, 20 MeV/c in B = 0.2 T.

  // v0 cuts
  Configurable<float> min_v0cospa_tpconly{"min_v0cospa_tpconly", 0.99, "min V0 CosPA to V0s with TPConly tracks"}; // double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<float> min_v0cospa_its{"min_v0cospa_its", 0.99, "min V0 CosPA to V0s with ITs hits"};               // double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<float> max_dcav0dau_tpconly{"max_dcav0dau_tpconly", 3.0, "max distance btween 2 legs to V0s with TPConly tracks"};
  Configurable<float> max_dcav0dau_its{"max_dcav0dau_its", 0.5, "max distance btween 2 legs to V0s with ITS hits"};
  Configurable<float> max_dcav0dau_itsibss{"max_dcav0dau_itsibss", 1.0, "max distance btween 2 legs to V0s with ITS hits on ITSib SS"};
  Configurable<float> max_dcav0dau_tpc_inner_fc{"max_dcav0dau_tpc_inner_fc", 1.5, "max distance btween 2 legs to V0s with ITS hits on TPC inner FC"};
  Configurable<float> min_v0radius{"min_v0radius", 1.0, "min v0 radius"};
  Configurable<float> max_v0radius{"max_v0radius", 90.0, "max v0 radius"};
  Configurable<float> margin_r_its{"margin_r_its", 3.0, "margin for r cut in cm"};
  Configurable<float> margin_r_tpc{"margin_r_tpc", 7.0, "margin for r cut in cm"};
  Configurable<float> margin_r_itstpc_tpc{"margin_r_itstpc_tpc", 7.0, "margin for r cut in cm"};
  Configurable<float> margin_z{"margin_z", 7.0, "margin for z cut in cm"};
  Configurable<float> max_alpha_ap{"max_alpha_ap", 0.95, "max alpha for AP cut"};
  Configurable<float> max_qt_ap{"max_qt_ap", 0.01, "max qT for AP cut"};
  Configurable<float> min_pt_v0{"min_pt_v0", 0.1, "min pT for v0 photons at PV"};
  Configurable<float> max_pt_v0_itsonly{"max_pt_v0_itsonly", 0.3, "max pT for v0 photons wth 2 ITSonly tracks at PV"};
  Configurable<float> max_eta_v0{"max_eta_v0", 0.9, "max eta for v0 photons at PV"};
  Configurable<float> kfMassConstrain{"kfMassConstrain", -1.f, "mass constrain for the KFParticle mother particle"};
  Configurable<float> max_r_req_its{"max_r_req_its", 16.0, "max Rxy for V0 with ITS hits"};
  Configurable<float> min_r_tpconly{"min_r_tpconly", 32.0, "min Rxy for V0 with TPConly tracks"};
  Configurable<float> max_r_itsmft_ss{"max_r_itsmft_ss", 66.0, "max Rxy for ITS/MFT SS"};
  Configurable<float> max_dcatopv_xy_v0{"max_dcatopv_xy_v0", +1e+10, "max. DCAxy to PV for V0"};
  Configurable<float> max_dcatopv_z_v0{"max_dcatopv_z_v0", +1e+10, "max. DCAz to PV for V0"};
  Configurable<bool> reject_v0_on_itsib{"reject_v0_on_itsib", true, "flag to reject v0s on ITSib"};

  // PCM ML inference
  Configurable<bool> applyPCMMl{"applyPCMMl", false, "Flag to apply ML selections"};
  Configurable<bool> use2DBinning{"use2DBinning", false, "Flag to enable/disable 2D binning for ML application"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};
  Configurable<int> nClassesPCMMl{"nClassesPCMMl", static_cast<int>(o2::analysis::em_cuts_ml::NCutScores), "Number of classes in ML model"};
  Configurable<int> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
  Configurable<std::string> centTypePCMMl{"centTypePCMMl", "CentFT0C", "Centrality type for 2D ML application: CentFT0C, CentFT0M, or CentFT0A"};
  Configurable<std::vector<int>> cutDirPCMMl{"cutDirPCMMl", std::vector<int>{o2::analysis::em_cuts_ml::vecCutDir}, "Whether to reject score values greater or smaller than the threshold"};
  Configurable<std::vector<std::string>> namesInputFeatures{"namesInputFeatures", std::vector<std::string>{"feature1", "feature2"}, "Names of ML model input features"};
  Configurable<std::vector<std::string>> modelPathsCCDB{"modelPathsCCDB", std::vector<std::string>{"path_ccdb/BDT_PCM/"}, "Paths of models on CCDB"};
  Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"ModelHandler_onnx_PCM.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
  Configurable<std::vector<std::string>> labelsBinsPCMMl{"labelsBinsPCMMl", std::vector<std::string>{"bin 0", "bin 1"}, "Labels for bins"};
  Configurable<std::vector<std::string>> labelsCutScoresPCMMl{"labelsCutScoresPCMMl", std::vector<std::string>{o2::analysis::em_cuts_ml::labelsCutScore}, "Labels for cut scores"};
  Configurable<std::vector<double>> binsPtPCMMl{"binsPtPCMMl", std::vector<double>{0.0, +1e+10}, "pT bin limits for ML application"};
  Configurable<std::vector<double>> binsCentPCMMl{"binsCentPCMMl", std::vector<double>{0.0, 100.0}, "Centrality bin limits for ML application"};
  Configurable<std::vector<double>> cutsPCMMlFlat{"cutsPCMMlFlat", {0.5}, "Flattened ML cuts: [bin0_score0, bin0_score1, ..., binN_scoreM]"};

  o2::analysis::EmMlResponsePCM<float> emMlResponse;
  std::vector<float> outputML;
  o2::ccdb::CcdbApi ccdbApi;

  int mRunNumber;
  float d_bz;
  float maxSnp;  // max sine phi for propagation
  float maxStep; // max step size (cm) for propagation
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut = nullptr;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  o2::aod::common::TPCVDriftManager mVDriftMgr;

  HistogramRegistry registry{
    "registry",
    {
      {"hCollisionCounter", "hCollisionCounter", {HistType::kTH1F, {{1, 0.5f, 1.5f}}}},
      {"V0/hAP", "Armenteros Podolanski;#alpha;q_{T} (GeV/c)", {HistType::kTH2F, {{200, -1.0f, 1.0f}, {250, 0, 0.25}}}},
      {"V0/hConversionPointXY", "conversion point in XY;X (cm);Y (cm)", {HistType::kTH2F, {{400, -100.0f, 100.0f}, {400, -100.f, 100.f}}}},
      {"V0/hConversionPointRZ", "conversion point in RZ;Z (cm);R_{xy} (cm)", {HistType::kTH2F, {{200, -100.0f, 100.0f}, {200, 0.f, 100.f}}}},
      {"V0/hPt", "pT of V0 at PV;p_{T,#gamma} (GeV/c)", {HistType::kTH1F, {{1000, 0.0f, 10.0f}}}},
      {"V0/hEtaPhi", "#eta vs. #varphi of V0 at PV;#varphi (rad.);#eta", {HistType::kTH2F, {{72, 0.0f, o2::constants::math::TwoPI}, {200, -1, +1}}}},
      {"V0/hCosPA", "cosine of pointing angle;cosine of pointing angle", {HistType::kTH1F, {{100, 0.99f, 1.f}}}},
      {"V0/hCosPA_Rxy", "cosine of pointing angle;r_{xy} (cm);cosine of pointing angle", {HistType::kTH2F, {{200, 0, 100}, {100, 0.99f, 1.f}}}},
      {"V0/hCosPAXY_Rxy", "cosine of pointing angle;r_{xy} (cm);cosine of pointing angle", {HistType::kTH2F, {{200, 0, 100}, {100, 0.99f, 1.f}}}},
      {"V0/hCosPARZ_Rxy", "cosine of pointing angle;r_{xy} (cm);cosine of pointing angle", {HistType::kTH2F, {{200, 0, 100}, {100, 0.99f, 1.f}}}},
      {"V0/hPCA", "distance between 2 legs at SV;PCA (cm)", {HistType::kTH1F, {{500, 0.0f, 5.f}}}},
      {"V0/hPCA_Rxy", "distance between 2 legs at SV;R_{xy} (cm);PCA (cm)", {HistType::kTH2F, {{200, 0, 100}, {500, 0.0f, 5.f}}}},
      {"V0/hPCA_CosPA", "distance between 2 legs at SV vs. cosPA;cosine of pointing angle;PCA (cm)", {HistType::kTH2F, {{100, 0.99, 1}, {500, 0.0f, 5.f}}}},
      {"V0/hDCAxyz", "DCA to PV;DCA_{xy} (cm);DCA_{z} (cm)", {HistType::kTH2F, {{200, -5.f, +5.f}, {200, -5.f, +5.f}}}},
      {"V0/hMeeSV_Rxy", "mee at SV vs. R_{xy};R_{xy} (cm);m_{ee} at SV (GeV/c^{2})", {HistType::kTH2F, {{200, 0.0f, 100.f}, {100, 0, 0.1f}}}},
      {"V0/hRxy_minX_ITSonly_ITSonly", "min trackiu X vs. R_{xy};trackiu X (cm);min trackiu X - R_{xy} (cm)", {HistType::kTH2F, {{100, 0.0f, 100.f}, {100, -50.0, 50.0f}}}},
      {"V0/hRxy_minX_ITSTPC_ITSTPC", "min trackiu X vs. R_{xy};trackiu X (cm);min trackiu X - R_{xy} (cm)", {HistType::kTH2F, {{100, 0.0f, 100.f}, {100, -50.0, 50.0f}}}},
      {"V0/hRxy_minX_ITSTPC_ITSonly", "min trackiu X vs. R_{xy};trackiu X (cm);min trackiu X - R_{xy} (cm)", {HistType::kTH2F, {{100, 0.0f, 100.f}, {100, -50.0, 50.0f}}}},
      {"V0/hRxy_minX_ITSTPC_TPC", "min trackiu X vs. R_{xy};trackiu X (cm);min trackiu X - R_{xy} (cm)", {HistType::kTH2F, {{100, 0.0f, 100.f}, {100, -50.0, 50.0f}}}},
      {"V0/hRxy_minX_TPC_TPC", "min trackiu X vs. R_{xy};trackiu X (cm);min trackiu X - R_{xy} (cm)", {HistType::kTH2F, {{100, 0.0f, 100.f}, {100, -50.0, 50.0f}}}},
      {"V0/hPCA_diffX", "PCA vs. trackiu X - R_{xy};distance btween 2 legs (cm);min trackiu X - R_{xy} (cm)", {HistType::kTH2F, {{500, 0.0f, 5.f}, {100, -50.0, 50.0f}}}},
      {"V0/hPhiV", "#phi_{V}; #phi_{V} (rad.)", {HistType::kTH1F, {{500, 0.0f, o2::constants::math::TwoPI}}}},
      {"V0/hBDTvalueBeforeCutVsPt", "BDT response before cut vs pT; pT (GeV/c); BDT response", {HistType::kTH2F, {{1000, 0.0f, 20.0f}, {1000, 0.0f, 1.0f}}}},
      {"V0/hBDTvalueAfterCutVsPt", "BDT response after cut vs pT; pT (GeV/c); BDT response", {HistType::kTH2F, {{1000, 0.0f, 20.0f}, {1000, 0.0f, 1.0f}}}},
      {"V0Leg/hPt", "pT of leg at SV;p_{T,e} (GeV/c)", {HistType::kTH1F, {{1000, 0.0f, 10.0f}}}},
      {"V0Leg/hEtaPhi", "#eta vs. #varphi of leg at SV;#varphi (rad.);#eta", {HistType::kTH2F, {{72, 0.0f, o2::constants::math::TwoPI}, {200, -1, +1}}}},
      {"V0Leg/hRelDeltaPt", "pT resolution;p_{T} (GeV/c);#Deltap_{T}/p_{T}", {HistType::kTH2F, {{1000, 0.f, 10.f}, {100, 0, 1}}}},
      {"V0Leg/hDCAxyz", "DCA xy vs. z to PV;DCA_{xy} (cm);DCA_{z} (cm)", {HistType::kTH2F, {{200, -50.f, 50.f}, {200, -50.f, +50.f}}}},
      {"V0Leg/hdEdx_Pin", "TPC dE/dx vs. p_{in};p_{in} (GeV/c);TPC dE/dx", {HistType::kTH2F, {{1000, 0.f, 10.f}, {200, 0.f, 200.f}}}},
      {"V0Leg/hTPCNsigmaEl", "TPC dE/dx vs. p_{in};p_{in} (GeV/c);n #sigma_{e}^{TPC}", {HistType::kTH2F, {{1000, 0.f, 10.f}, {100, -5.f, +5.f}}}},
      {"V0Leg/hXZ", "track iu x vs. z;z (cm);x (cm)", {HistType::kTH2F, {{200, -100.f, 100.f}, {200, 0.f, 100.f}}}},
    }};

  void init(InitContext&)
  {
    mRunNumber = 0;
    d_bz = 0;
    maxSnp = 0.85f;  // could be changed later
    maxStep = 2.00f; // could be changed later

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    if (useMatCorrType == MatCorrType::TGeo) {
      LOGF(info, "TGeo correction requested, loading geometry");
      if (!o2::base::GeometryManager::isGeometryLoaded()) {
        ccdb->get<TGeoManager>(geoPath);
      }
    }
    if (useMatCorrType == MatCorrType::LUT) {
      LOGF(info, "LUT correction requested, loading LUT");
      lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(lutPath));
    }

    if (useMatCorrType == MatCorrType::TGeo) {
      matCorr = o2::base::Propagator::MatCorrType::USEMatCorrTGeo;
    }
    if (useMatCorrType == MatCorrType::LUT) {
      matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
    }

    if (applyPCMMl) {
      if (use2DBinning) {
        int binsNPt = static_cast<int>(binsPtPCMMl->size()) - 1;
        int binsNCent = static_cast<int>(binsCentPCMMl->size()) - 1;
        int binsN = binsNPt * binsNCent;
        if (binsN * static_cast<int>(cutDirPCMMl->size()) != static_cast<int>(cutsPCMMlFlat->size())) {
          LOG(fatal) << "Mismatch in number of bins and cuts provided for 2D ML application: binsN * cutDirPCMMl: " << int(binsN) * int(cutDirPCMMl->size()) << " bins vs. cutsPCMMlFlat: " << cutsPCMMlFlat->size() << " cuts";
        }
        if (binsN != static_cast<int>(onnxFileNames->size())) {
          LOG(fatal) << "Mismatch in number of bins and ONNX files provided for 2D ML application: binsN " << binsN << " bins vs. onnxFileNames: " << onnxFileNames->size() << " ONNX files";
        }
        if (binsN != static_cast<int>(labelsBinsPCMMl->size())) {
          LOG(fatal) << "Mismatch in number of bins and labels provided for 2D ML application: binsN:" << binsN << " bins vs. labelsBinsPCMMl: " << labelsBinsPCMMl->size() << " labels";
        }
        if (static_cast<int>(cutDirPCMMl->size()) != nClassesPCMMl) {
          LOG(fatal) << "Mismatch in number of classes and cut directions provided for 2D ML application: nClassesPCMMl: " << nClassesPCMMl << " classes vs. cutDirPCMMl: " << cutDirPCMMl->size() << " cut directions";
        }
        if (static_cast<int>(labelsCutScoresPCMMl->size()) != nClassesPCMMl) {
          LOG(fatal) << "Mismatch in number of labels for cut scores and number of classes provided for 2D ML application: nClassesPCMMl: " << nClassesPCMMl << " classes vs. labelsCutScoresPCMMl: " << labelsCutScoresPCMMl->size() << " labels";
        }
        LabeledArray<double> cutsPCMMl(cutsPCMMlFlat->data(), binsN, nClassesPCMMl, labelsBinsPCMMl, labelsCutScoresPCMMl);
        emMlResponse.configure2D(binsPtPCMMl, binsCentPCMMl, cutsPCMMl, cutDirPCMMl, nClassesPCMMl);
      } else {
        int binsNPt = static_cast<int>(binsPtPCMMl->size()) - 1;
        if (binsNPt * static_cast<int>(cutDirPCMMl->size()) != static_cast<int>(cutsPCMMlFlat->size())) {
          LOG(fatal) << "Mismatch in number of pT bins and cuts provided for ML application: binsNPt * cutDirPCMMl:" << binsNPt * cutDirPCMMl->size() << " bins vs. cutsPCMMlFlat: " << cutsPCMMlFlat->size() << " cuts";
        }
        if (binsNPt != static_cast<int>(onnxFileNames->size())) {
          LOG(fatal) << "Mismatch in number of pT bins and ONNX files provided for ML application: binsNPt " << binsNPt << " bins vs. onnxFileNames: " << onnxFileNames->size() << " ONNX files";
        }
        if (binsNPt != static_cast<int>(labelsBinsPCMMl->size())) {
          LOG(fatal) << "Mismatch in number of pT bins and labels provided for ML application: binsNPt:" << binsNPt << " bins vs. labelsBinsPCMMl: " << labelsBinsPCMMl->size() << " labels";
        }
        if (nClassesPCMMl != static_cast<int>(cutDirPCMMl->size())) {
          LOG(fatal) << "Mismatch in number of classes and cut directions provided for ML application: nClassesPCMMl: " << nClassesPCMMl << " classes vs. cutDirPCMMl: " << cutDirPCMMl->size() << " cut directions";
        }
        if (static_cast<int>(labelsCutScoresPCMMl->size()) != nClassesPCMMl) {
          LOG(fatal) << "Mismatch in number of labels for cut scores and number of classes provided for ML application: nClassesPCMMl:" << nClassesPCMMl << " classes vs. labelsCutScoresPCMMl: " << labelsCutScoresPCMMl->size() << " labels";
        }
        LabeledArray<double> cutsPCMMl(cutsPCMMlFlat->data(), binsNPt, nClassesPCMMl, labelsBinsPCMMl, labelsCutScoresPCMMl);
        emMlResponse.configure(binsPtPCMMl, cutsPCMMl, cutDirPCMMl, nClassesPCMMl);
      }
      if (loadModelsFromCCDB) {
        ccdbApi.init(ccdburl);
        emMlResponse.setModelPathsCCDB(onnxFileNames, ccdbApi, modelPathsCCDB, timestampCCDB);
      } else {
        emMlResponse.setModelPathsLocal(onnxFileNames);
      }
      emMlResponse.cacheInputFeaturesIndices(namesInputFeatures);
      emMlResponse.init();
    }
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
      if (std::fabs(d_bz) > 1e-5) {
        grpmag.setL3Current(30000.f / (d_bz / 5.0f));
      }
      o2::base::Propagator::initFieldFromGRP(&grpmag);
      mRunNumber = bc.runNumber();
      return;
    }

    auto run3grp_timestamp = bc.timestamp();
    o2::parameters::GRPObject* grpo = nullptr;
    o2::parameters::GRPMagField* grpmag = nullptr;
    if (!skipGRPOquery) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    }
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

    if (useMatCorrType == 2) {
      // setMatLUT only after magfield has been initalized (setMatLUT has implicit and problematic init field call if not)
      o2::base::Propagator::Instance()->setMatLUT(lut);
    }
    /// Set magnetic field for KF vertexing
    const float magneticField = o2::base::Propagator::Instance()->getNominalBz();
    KFParticle::SetField(magneticField);

    mVDriftMgr.init(&ccdb->instance());
  }

  void updateCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    auto timestamp = bc.timestamp();

    mVDriftMgr.update(timestamp);
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

    if (disableITSonlyTracks && isITSonlyTrack(track)) {
      return false;
    }

    if (disableTPConlyTracks && isTPConlyTrack(track)) {
      return false;
    }

    if (requireITShit && !track.hasITS()) {
      return false;
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
      if (track.tpcFractionSharedCls() > max_frac_shared_clusters_tpc) {
        return false;
      }
      if (track.tpcNSigmaEl() < minTPCNsigmaEl || maxTPCNsigmaEl < track.tpcNSigmaEl()) {
        return false;
      }
    }

    if (track.hasITS()) {
      if (track.itsChi2NCl() > maxchi2its) {
        return false;
      }

      if (reject_v0_on_itsib) {
        auto hits_ib = std::count_if(its_ib_Requirement.second.begin(), its_ib_Requirement.second.end(), [&](auto&& requiredLayer) { return track.itsClusterMap() & (1 << requiredLayer); });
        bool its_ob_only = hits_ib <= its_ib_Requirement.first;
        if (!its_ob_only) {
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
    if (cospaXY < -1.f) {
      return -1.f;
    } else if (cospaXY > 1.f) {
      return 1.f;
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
    if (cospaRZ < -1.f) {
      return -1.f;
    } else if (cospaRZ > 1.f) {
      return 1.f;
    }
    return cospaRZ;
  }

  template <bool isMC, typename TTrack, typename TShiftedTrack, typename TKFParticle>
  void fillTrackTable(TTrack const& track, TShiftedTrack const& shiftedtrack, TKFParticle const& kfp, const float dcaXY, const float dcaZ)
  {
    v0legs(track.collisionId(), track.globalIndex(), track.sign(),
           kfp.GetPx(), kfp.GetPy(), kfp.GetPz(), dcaXY, dcaZ,
           track.tpcNClsFindable(), track.tpcNClsFindableMinusFound(), track.tpcNClsFindableMinusCrossedRows(), track.tpcNClsShared(),
           track.tpcChi2NCl(), track.tpcInnerParam(), track.tpcSignal(),
           track.tpcNSigmaEl(), track.tpcNSigmaPi(),
           track.itsClusterSizes(), track.itsChi2NCl(), track.detectorMap());

    v0legsXYZ(shiftedtrack.getX(), shiftedtrack.getY(), shiftedtrack.getZ());

    if constexpr (isMC) {
      v0legsDeDxMC(track.mcTunedTPCSignal());
    }
  }

  template <bool isMC, class TBCs, class TCollisions, class TTracks, typename TV0>
  void fillV0Table(TV0 const& v0, const bool filltable)
  {
    // Get tracks
    const auto& pos = v0.template posTrack_as<TTracks>();
    const auto& ele = v0.template negTrack_as<TTracks>();
    const auto& collision = v0.template collision_as<TCollisions>(); // collision where this v0 belongs to.
    // LOGF(info, "v0.collisionId() = %d, pos.collisionId() = %d, ele.collisionId() = %d", v0.collisionId(), pos.collisionId(), ele.collisionId());

    if (pos.sign() * ele.sign() > 0) { // reject same sign pair
      return;
    }

    if (pos.pt() < min_pt_trackiu || ele.pt() < min_pt_trackiu) {
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

    // if(isTPConlyTrack(ele)){
    //   // LOGF(info, "TPConly: ele.globalIndex() = %d, ele.x() = %f, ele.y() = %f, ele.z() = %f, ele.tgl() = %f, ele.alpha() = %f, ele.snp() = %f, ele.signed1Pt() = %f", ele.globalIndex(), ele.x(), ele.y(), ele.z(), ele.tgl(), ele.alpha(), ele.snp(), ele.signed1Pt());
    //   // LOGF(info, "TPConly: ele.globalIndex() = %d, ele.cYY() = %f, ele.cZY() = %f, ele.cZZ() = %f, ele.cSnpY() = %f, ele.cSnpZ() = %f, ele.cSnpSnp() = %f, ele.cTglY() = %f, ele.cTglZ() = %f, ele.cTglSnp() = %f, ele.cTglTgl() = %f, ele.c1PtY() = %f, ele.c1PtZ() = %f, ele.c1PtSnp() = %f, ele.c1PtTgl() = %f, ele.c1Pt21Pt2() = %f", ele.globalIndex(), ele.cYY(), ele.cZY(), ele.cZZ(), ele.cSnpY(), ele.cSnpZ(), ele.cSnpSnp(), ele.cTglY(), ele.cTglZ(), ele.cTglSnp(), ele.cTglTgl(), ele.c1PtY(), ele.c1PtZ(), ele.c1PtSnp(), ele.c1PtTgl(), ele.c1Pt21Pt2());
    // }

    // Calculate DCA with respect to the collision associated to the v0, not individual tracks
    std::array<float, 2> dcaInfo;

    auto pTrack = getTrackParCov(pos);
    if (moveTPCTracks && isTPConlyTrack(pos) && !mVDriftMgr.moveTPCTrack<TBCs, TCollisions>(collision, pos, pTrack)) {
      LOGP(error, "failed correction for positive tpc track");
      return;
    }
    auto pTrackC = pTrack;
    pTrackC.setPID(o2::track::PID::Electron);
    o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, pTrackC, 2.f, matCorr, &dcaInfo);
    auto posdcaXY = dcaInfo[0];
    auto posdcaZ = dcaInfo[1];

    auto nTrack = getTrackParCov(ele);
    if (moveTPCTracks && isTPConlyTrack(ele) && !mVDriftMgr.moveTPCTrack<TBCs, TCollisions>(collision, ele, nTrack)) {
      LOGP(error, "failed correction for negative tpc track");
      return;
    }
    auto nTrackC = nTrack;
    nTrackC.setPID(o2::track::PID::Electron);
    o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, nTrackC, 2.f, matCorr, &dcaInfo);
    auto eledcaXY = dcaInfo[0];
    auto eledcaZ = dcaInfo[1];

    if (std::fabs(posdcaXY) < dcapostopv || std::fabs(eledcaXY) < dcanegtopv) {
      return;
    }

    float xyz[3] = {0.f, 0.f, 0.f};
    Vtx_recalculationParCov(o2::base::Propagator::Instance(), pTrack, nTrack, xyz, matCorr);
    float rxy_tmp = RecoDecay::sqrtSumOfSquares(xyz[0], xyz[1]);
    if (rxy_tmp > maxX + margin_r_tpc) {
      return;
    }
    if (rxy_tmp < std::fabs(xyz[2]) * std::tan(2 * std::atan(std::exp(-max_eta_v0))) - margin_z) {
      return; // RZ line cut
    }

    KFPTrack kfp_track_pos = createKFPTrackFromTrackParCov(pTrack, pos.sign(), pos.tpcNClsFound(), pos.tpcChi2NCl());
    KFPTrack kfp_track_ele = createKFPTrackFromTrackParCov(nTrack, ele.sign(), ele.tpcNClsFound(), ele.tpcChi2NCl());
    KFParticle kfp_pos(kfp_track_pos, kPositron);
    KFParticle kfp_ele(kfp_track_ele, kElectron);
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
    if (rxy > maxX + margin_r_tpc) {
      return;
    }
    if (rxy < std::fabs(gammaKF_DecayVtx.GetZ()) * std::tan(2 * std::atan(std::exp(-max_eta_v0))) - margin_z) {
      return; // RZ line cut
    }
    if (rxy < min_v0radius || max_v0radius < rxy) {
      return;
    }

    if (!filltable) {
      if (isITSTPCTrack(pos) && isITSTPCTrack(ele)) {
        registry.fill(HIST("V0/hRxy_minX_ITSTPC_ITSTPC"), std::min(pTrack.getX(), nTrack.getX()), std::min(pTrack.getX(), nTrack.getX()) - rxy); // trackiu.x() - rxy should be positive
      } else if (isITSonlyTrack(pos) && isITSonlyTrack(ele)) {
        registry.fill(HIST("V0/hRxy_minX_ITSonly_ITSonly"), std::min(pTrack.getX(), nTrack.getX()), std::min(pTrack.getX(), nTrack.getX()) - rxy); // trackiu.x() - rxy should be positive
      } else if ((isITSTPCTrack(pos) && isITSonlyTrack(ele)) || (isITSTPCTrack(ele) && isITSonlyTrack(pos))) {
        registry.fill(HIST("V0/hRxy_minX_ITSTPC_ITSonly"), std::min(pTrack.getX(), nTrack.getX()), std::min(pTrack.getX(), nTrack.getX()) - rxy); // trackiu.x() - rxy should be positive
      } else if (isITSTPCTrack(pos) && !ele.hasITS()) {
        registry.fill(HIST("V0/hRxy_minX_ITSTPC_TPC"), std::min(pTrack.getX(), 83.f), std::min(pTrack.getX(), 83.f) - rxy); // trackiu.x() - rxy should be positive
      } else if (isITSTPCTrack(ele) && !pos.hasITS()) {
        registry.fill(HIST("V0/hRxy_minX_ITSTPC_TPC"), std::min(nTrack.getX(), 83.f), std::min(nTrack.getX(), 83.f) - rxy); // trackiu.x() - rxy should be positive
      } else {
        registry.fill(HIST("V0/hRxy_minX_TPC_TPC"), std::min(83.f, 83.f), std::min(83.f, 83.f) - rxy); // trackiu.x() - rxy should be positive
      }
    }

    if (pos.hasITS() && ele.hasITS()) { // ITSonly-ITSonly, ITSTPC-ITSTPC, ITSTPC-ITSonly
      if (rxy > std::min(pTrack.getX(), nTrack.getX()) + margin_r_its) {
        return;
      }
    } else if (!pos.hasITS() && ele.hasITS()) { // ITSTPC-TPC
      if (rxy > std::min(83.f, nTrack.getX()) + margin_r_itstpc_tpc) {
        return;
      }
    } else if (pos.hasITS() && !ele.hasITS()) { // ITSTPC-TPC
      if (rxy > std::min(pTrack.getX(), 83.f) + margin_r_itstpc_tpc) {
        return;
      }
    } else if (!pos.hasITS() && !ele.hasITS()) { // TPC-TPC
      if (rxy > std::min(83.f, 83.f) + margin_r_tpc) {
        return;
      }
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
    float v0phi = RecoDecay::constrainAngle(RecoDecay::phi(gammaKF_PV.GetPx(), gammaKF_PV.GetPy()));

    // KFParticle gammaKF_DecayVtx2 = gammaKF;
    // gammaKF_DecayVtx2.SetProductionVertex(KFPV);
    // gammaKF_DecayVtx2.TransportToPoint(xyz);
    // LOGF(info, "gammaKF_PV.GetPx() = %f, gammaKF_DecayVtx.GetPx() = %f, gammaKF_DecayVtx2.GetPx() = %f", gammaKF_PV.GetPx(), gammaKF_DecayVtx.GetPx(), gammaKF_DecayVtx2.GetPx());
    // LOGF(info, "gammaKF_PV.GetPy() = %f, gammaKF_DecayVtx.GetPy() = %f, gammaKF_DecayVtx2.GetPy() = %f", gammaKF_PV.GetPy(), gammaKF_DecayVtx.GetPy(), gammaKF_DecayVtx2.GetPy());
    // LOGF(info, "gammaKF_PV.GetPz() = %f, gammaKF_DecayVtx.GetPz() = %f, gammaKF_DecayVtx2.GetPz() = %f", gammaKF_PV.GetPz(), gammaKF_DecayVtx.GetPz(), gammaKF_DecayVtx2.GetPz());

    if (std::fabs(v0eta) > max_eta_v0 || v0pt < min_pt_v0) {
      return;
    }

    if (isITSonlyTrack(ele) && isITSonlyTrack(pos) && v0pt > max_pt_v0_itsonly) {
      return;
    }

    KFParticle kfp_pos_DecayVtx = kfp_pos;  // Don't set Primary Vertex
    KFParticle kfp_ele_DecayVtx = kfp_ele;  // Don't set Primary Vertex
    kfp_pos_DecayVtx.TransportToPoint(xyz); // Don't set Primary Vertex
    kfp_ele_DecayVtx.TransportToPoint(xyz); // Don't set Primary Vertex

    V0PhotonCandidate v0photoncandidate(gammaKF_DecayVtx, kfp_pos_DecayVtx, kfp_ele_DecayVtx, collision, cospa_kf, d_bz);

    if (!ele.hasITS() && !pos.hasITS()) { // V0s with TPConly-TPConly
      if (max_r_itsmft_ss < rxy && rxy < maxX + margin_r_tpc) {
        if (v0photoncandidate.GetPCA() > max_dcav0dau_tpc_inner_fc) {
          return;
        }
      } else {
        if (v0photoncandidate.GetPCA() > max_dcav0dau_tpconly) {
          return;
        }
      }
    } else { // V0s with ITS hits
      if (rxy < max_r_req_its) {
        if (v0photoncandidate.GetPCA() > max_dcav0dau_itsibss) {
          return;
        }
      } else {
        if (v0photoncandidate.GetPCA() > max_dcav0dau_its) {
          return;
        }
      }
    }

    if (isITSonlyTrack(pos) && v0photoncandidate.GetPosPt() > maxpt_itsonly) {
      return;
    }

    if (isITSonlyTrack(ele) && v0photoncandidate.GetElePt() > maxpt_itsonly) {
      return;
    }

    if (v0photoncandidate.GetChi2NDF() > 6e+3) { // protection for uint16.
      return;
    }

    if (std::fabs(v0photoncandidate.GetDcaXYToPV()) > max_dcatopv_xy_v0 || std::fabs(v0photoncandidate.GetDcaZToPV()) > max_dcatopv_z_v0) {
      return;
    }

    if (!checkAP(v0photoncandidate.GetAlpha(), v0photoncandidate.GetQt(), max_alpha_ap, max_qt_ap)) { // store only photon conversions
      return;
    }
    pca_map[std::make_tuple(v0.globalIndex(), collision.globalIndex(), pos.globalIndex(), ele.globalIndex())] = v0photoncandidate.GetPCA();
    cospa_map[std::make_tuple(v0.globalIndex(), collision.globalIndex(), pos.globalIndex(), ele.globalIndex())] = v0photoncandidate.GetCosPA();

    if (applyPCMMl) {
      bool isSelectedML = false;
      std::vector<float> mlInputFeatures = emMlResponse.getInputFeatures(v0photoncandidate, pos, ele);
      if (use2DBinning) {
        if (std::string(centTypePCMMl) == "CentFT0C") {
          isSelectedML = emMlResponse.isSelectedMl(mlInputFeatures, v0photoncandidate.GetPt(), v0photoncandidate.GetCentFT0C(), outputML);
        } else if (std::string(centTypePCMMl) == "CentFT0A") {
          isSelectedML = emMlResponse.isSelectedMl(mlInputFeatures, v0photoncandidate.GetPt(), v0photoncandidate.GetCentFT0A(), outputML);
        } else if (std::string(centTypePCMMl) == "CentFT0M") {
          isSelectedML = emMlResponse.isSelectedMl(mlInputFeatures, v0photoncandidate.GetPt(), v0photoncandidate.GetCentFT0M(), outputML);
        } else {
          LOG(fatal) << "Unsupported centTypePCMMl: " << centTypePCMMl << " , please choose from CentFT0C, CentFT0A, CentFT0M.";
        }
      } else {
        isSelectedML = emMlResponse.isSelectedMl(mlInputFeatures, v0photoncandidate.GetPt(), outputML);
      }
      if (filltable) {
        registry.fill(HIST("V0/hBDTvalueBeforeCutVsPt"), v0photoncandidate.GetPt(), outputML[0]);
      }
      if (!isSelectedML) {
        return;
      }
      if (filltable) {
        registry.fill(HIST("V0/hBDTvalueAfterCutVsPt"), v0photoncandidate.GetPt(), outputML[0]);
      }
    }

    if (filltable) {
      registry.fill(HIST("V0/hAP"), v0photoncandidate.GetAlpha(), v0photoncandidate.GetQt());
      registry.fill(HIST("V0/hConversionPointXY"), gammaKF_DecayVtx.GetX(), gammaKF_DecayVtx.GetY());
      registry.fill(HIST("V0/hConversionPointRZ"), gammaKF_DecayVtx.GetZ(), rxy);
      registry.fill(HIST("V0/hPt"), v0photoncandidate.GetPt());
      registry.fill(HIST("V0/hEtaPhi"), v0phi, v0eta);
      registry.fill(HIST("V0/hCosPA"), v0photoncandidate.GetCosPA());
      registry.fill(HIST("V0/hCosPA_Rxy"), rxy, v0photoncandidate.GetCosPA());
      registry.fill(HIST("V0/hPCA"), v0photoncandidate.GetPCA());
      registry.fill(HIST("V0/hPCA_CosPA"), v0photoncandidate.GetCosPA(), v0photoncandidate.GetPCA());
      registry.fill(HIST("V0/hPCA_Rxy"), rxy, v0photoncandidate.GetPCA());
      registry.fill(HIST("V0/hDCAxyz"), v0photoncandidate.GetDcaXYToPV(), v0photoncandidate.GetDcaZToPV());
      registry.fill(HIST("V0/hPCA_diffX"), v0photoncandidate.GetPCA(), std::min(pTrack.getX(), nTrack.getX()) - rxy); // trackiu.x() - rxy should be positive
      registry.fill(HIST("V0/hPhiV"), v0photoncandidate.GetPhiV());

      float cospaXY_kf = cospaXY_KF(gammaKF_DecayVtx, KFPV);
      float cospaRZ_kf = cospaRZ_KF(gammaKF_DecayVtx, KFPV);
      // LOGF(info, "cospa_kf = %f, cospaXY_kf = %f, cospaRZ_kf = %f", cospa_kf, cospaXY_kf, cospaRZ_kf);
      registry.fill(HIST("V0/hCosPAXY_Rxy"), rxy, cospaXY_kf);
      registry.fill(HIST("V0/hCosPARZ_Rxy"), rxy, cospaRZ_kf);

      for (const auto& leg : {kfp_pos_DecayVtx, kfp_ele_DecayVtx}) {
        float legpt = RecoDecay::sqrtSumOfSquares(leg.GetPx(), leg.GetPy());
        float legeta = RecoDecay::eta(std::array{leg.GetPx(), leg.GetPy(), leg.GetPz()});
        float legphi = RecoDecay::constrainAngle(RecoDecay::phi(leg.GetPx(), leg.GetPy()));
        registry.fill(HIST("V0Leg/hPt"), legpt);
        registry.fill(HIST("V0Leg/hEtaPhi"), legphi, legeta);
      } // end of leg loop
      for (const auto& leg : {pos, ele}) {
        registry.fill(HIST("V0Leg/hdEdx_Pin"), leg.tpcInnerParam(), leg.tpcSignal());
        registry.fill(HIST("V0Leg/hTPCNsigmaEl"), leg.tpcInnerParam(), leg.tpcNSigmaEl());
      } // end of leg loop
      for (const auto& leg : {pTrack, nTrack}) {
        registry.fill(HIST("V0Leg/hXZ"), leg.getZ(), leg.getX());
        registry.fill(HIST("V0Leg/hRelDeltaPt"), leg.getPt(), leg.getPt() * std::sqrt(leg.getSigma1Pt2()));
      } // end of leg loop
      registry.fill(HIST("V0Leg/hDCAxyz"), posdcaXY, posdcaZ);
      registry.fill(HIST("V0Leg/hDCAxyz"), eledcaXY, eledcaZ);

      ROOT::Math::PxPyPzMVector vpos_sv(kfp_pos_DecayVtx.GetPx(), kfp_pos_DecayVtx.GetPy(), kfp_pos_DecayVtx.GetPz(), o2::constants::physics::MassElectron);
      ROOT::Math::PxPyPzMVector vele_sv(kfp_ele_DecayVtx.GetPx(), kfp_ele_DecayVtx.GetPy(), kfp_ele_DecayVtx.GetPz(), o2::constants::physics::MassElectron);
      ROOT::Math::PxPyPzMVector v0_sv = vpos_sv + vele_sv;
      registry.fill(HIST("V0/hMeeSV_Rxy"), rxy, v0_sv.M());

      v0photonskf(collision.globalIndex(), v0.globalIndex(), v0legs.lastIndex() + 1, v0legs.lastIndex() + 2,
                  gammaKF_DecayVtx.GetX(), gammaKF_DecayVtx.GetY(), gammaKF_DecayVtx.GetZ(),
                  gammaKF_PV.GetPx(), gammaKF_PV.GetPy(), gammaKF_PV.GetPz(),
                  v0_sv.M(), v0photoncandidate.GetDcaXYToPV(), v0photoncandidate.GetDcaZToPV(),
                  cospa_kf, cospaXY_kf, cospaRZ_kf,
                  v0photoncandidate.GetPCA(), v0photoncandidate.GetAlpha(), v0photoncandidate.GetQt(), v0photoncandidate.GetChi2NDF());
      v0photonsphiv(v0photoncandidate.GetPhiV());

      // v0photonskfcov(gammaKF_PV.GetCovariance(9), gammaKF_PV.GetCovariance(14), gammaKF_PV.GetCovariance(20), gammaKF_PV.GetCovariance(13), gammaKF_PV.GetCovariance(19), gammaKF_PV.GetCovariance(18));

      fillTrackTable<isMC>(pos, pTrack, kfp_pos_DecayVtx, posdcaXY, posdcaZ); // positive leg first
      fillTrackTable<isMC>(ele, nTrack, kfp_ele_DecayVtx, eledcaXY, eledcaZ); // negative leg second
    } // end of fill table
  }

  Preslice<aod::V0s> perCollision = o2::aod::v0::collisionId;
  std::map<std::tuple<int64_t, int64_t, int64_t, int64_t>, float> pca_map;      // (v0.globalIndex(), collision.globalIndex(), pos.globalIndex(), ele.globalIndex()) -> pca
  std::map<std::tuple<int64_t, int64_t, int64_t, int64_t>, float> cospa_map;    // (v0.globalIndex(), collision.globalIndex(), pos.globalIndex(), ele.globalIndex()) -> cospa
  std::vector<std::pair<int64_t, int64_t>> stored_v0Ids;                        // (pos.globalIndex(), ele.globalIndex())
  std::vector<std::tuple<int64_t, int64_t, int64_t, int64_t>> stored_fullv0Ids; // (v0.globalIndex(), collision.globalIndex(), pos.globalIndex(), ele.globalIndex())
  std::unordered_map<int64_t, int> nv0_map;                                     // map collisionId -> nv0

  template <bool isMC, bool isTriggerAnalysis, bool enableFilter, typename TCollisions, typename TV0s, typename TTracks, typename TBCs>
  void build(TCollisions const& collisions, TV0s const& v0s, TTracks const&, TBCs const&)
  {
    for (const auto& collision : collisions) {
      if constexpr (isMC) {
        if (!collision.has_mcCollision()) {
          continue;
        }
      }

      if (!collision.isSelected()) {
        continue;
      }

      if constexpr (isTriggerAnalysis) {
        if (collision.swtaliastmp_raw() == 0) {
          continue;
        }
      }

      nv0_map[collision.globalIndex()] = 0;

      const auto& bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      registry.fill(HIST("hCollisionCounter"), 1);

      updateCCDB(bc); // delay update until is needed

      const auto& v0s_per_coll = v0s.sliceBy(perCollision, collision.globalIndex());
      // LOGF(info, "n v0 = %d", v0s_per_coll.size());
      for (const auto& v0 : v0s_per_coll) {
        // LOGF(info, "collision.globalIndex() = %d, v0.globalIndex() = %d, v0.posTrackId() = %d, v0.negTrackId() = %d", collision.globalIndex(), v0.globalIndex(), v0.posTrackId() , v0.negTrackId());
        fillV0Table<isMC, TBCs, TCollisions, TTracks>(v0, false);
      } // end of v0 loop
    } // end of collision loop

    stored_v0Ids.reserve(pca_map.size());     // number of photon candidates per DF
    stored_fullv0Ids.reserve(pca_map.size()); // number of photon candidates per DF

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
        // auto v0 = v0s.rawIteratorAt(v0Id);
        // auto collision = collisions.rawIteratorAt(collisionId);
        // auto pos = tracks.rawIteratorAt(posId);
        // auto ele = tracks.rawIteratorAt(eleId);
        // LOGF(info, "!accept! | collision id = %d | v0id1 = %d , posid1 = %d , eleid1 = %d , pca1 = %f , cospa = %f", collisionId, v0Id, posId, eleId, v0pca, cospa);

        // fillV0Table<isMC, TCollisions, TTracks>(v0, true);
        stored_v0Ids.emplace_back(std::make_pair(posId, eleId));
        stored_fullv0Ids.emplace_back(std::make_tuple(v0Id, collisionId, posId, eleId));
        nv0_map[collisionId]++;
      }
    } // end of pca_map loop
    // LOGF(info, "pca_map.size() = %d", pca_map.size());

    for (const auto& fullv0Id : stored_fullv0Ids) {
      auto v0Id = std::get<0>(fullv0Id);
      // auto collisionId = std::get<1>(fullv0Id);
      // auto posId = std::get<2>(fullv0Id);
      // auto eleId = std::get<3>(fullv0Id);
      // LOGF(info, "!accept! | collision id = %d | v0id = %d , posid = %d , eleid = %d", collisionId, v0Id, posId, eleId);

      auto v0 = v0s.rawIteratorAt(v0Id);
      if constexpr (enableFilter) {
        auto collision_tmp = v0.template collision_as<TCollisions>(); // collision where this v0 belongs.
        if (!(collision_tmp.neeuls() >= 1 || collision_tmp.neeuls() + nv0_map[collision_tmp.globalIndex()] >= 2)) {
          continue;
        }
        // LOGF(info, "collision_tmp.globalIndex() = %d, collision_tmp.neeuls() = %d, nv0_map = %d", collision_tmp.globalIndex(), collision_tmp.neeuls(), nv0_map[collision_tmp.globalIndex()]);
      }

      fillV0Table<isMC, TBCs, TCollisions, TTracks>(v0, true);
    } // end of fullv0Id loop

    for (const auto& collision : collisions) {
      if constexpr (isMC) {
        if (!collision.has_mcCollision()) {
          continue;
        }
      }
      // events_ngpcm(nv0_map[collision.globalIndex()]);
    } // end of collision loop

    pca_map.clear();
    cospa_map.clear();
    nv0_map.clear();
    stored_v0Ids.clear();
    stored_v0Ids.shrink_to_fit();
    stored_fullv0Ids.clear();
    stored_fullv0Ids.shrink_to_fit();
  } // end of build

  //! type of V0. 0: built solely for cascades (does not pass standard V0 cuts), 1: standard 2, 3: photon-like with TPC-only use. Regular analysis should always use type 1 or 3.
  Filter v0Filter = o2::aod::v0::v0Type > (uint8_t)0;
  using filteredV0s = soa::Filtered<aod::V0s>;

  void processRec(MyCollisions const& collisions, filteredV0s const& v0s, MyTracksIU const& tracks, aod::BCsWithTimestamps const& bcs)
  {
    build<false, false, false>(collisions, v0s, tracks, bcs);
  }
  PROCESS_SWITCH(PhotonConversionBuilder, processRec, "process reconstructed info for data", true);

  void processRec_SWT(MyCollisionsWithSWT const& collisions, filteredV0s const& v0s, MyTracksIU const& tracks, aod::BCsWithTimestamps const& bcs)
  {
    build<false, true, false>(collisions, v0s, tracks, bcs);
  }
  PROCESS_SWITCH(PhotonConversionBuilder, processRec_SWT, "process reconstructed info for data", false);

  void processMC(MyCollisionsMC const& collisions, filteredV0s const& v0s, MyTracksIUMC const& tracks, aod::BCsWithTimestamps const& bcs)
  {
    build<true, false, false>(collisions, v0s, tracks, bcs);
  }
  PROCESS_SWITCH(PhotonConversionBuilder, processMC, "process reconstructed info for MC", false);

  void processRec_OnlyIfDielectron(soa::Join<MyCollisions, aod::EMEventsNee> const& collisions, filteredV0s const& v0s, MyTracksIU const& tracks, aod::BCsWithTimestamps const& bcs)
  {
    build<false, false, true>(collisions, v0s, tracks, bcs);
  }
  PROCESS_SWITCH(PhotonConversionBuilder, processRec_OnlyIfDielectron, "process reconstructed info for data", false);

  void processRec_SWT_OnlyIfDielectron(soa::Join<MyCollisionsWithSWT, aod::EMEventsNee> const& collisions, filteredV0s const& v0s, MyTracksIU const& tracks, aod::BCsWithTimestamps const& bcs)
  {
    build<false, true, true>(collisions, v0s, tracks, bcs);
  }
  PROCESS_SWITCH(PhotonConversionBuilder, processRec_SWT_OnlyIfDielectron, "process reconstructed info for data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<PhotonConversionBuilder>(cfgc, TaskName{"photon-conversion-builder"})};
}
