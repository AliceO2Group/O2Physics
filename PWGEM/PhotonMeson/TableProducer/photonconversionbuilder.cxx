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
#define HomogeneousField // o2-linter: disable=name/macro (name coming from KFParticle, not us, needed for KFParticle::SetField)
#endif

#include "PWGEM/Dilepton/Utils/PairUtilities.h"
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
#include <Framework/Concepts.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/HelixHelper.h>
#include <ReconstructionDataFormats/PID.h>

#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TPDGCode.h>

#include <KFPTrack.h>
#include <KFPVertex.h>
#include <KFParticle.h>

#include <GPUROOTCartesianFwd.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <set>
#include <string>
#include <tuple>
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
// using MyCollisionsWithSWT = soa::Join<MyCollisions, aod::EMSWTriggerBitsTMP>;
using MyCollisionsMC = soa::Join<MyCollisions, aod::McCollisionLabels>;

using MyTracksIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCFullEl, aod::pidTPCFullPi>;
using MyTracksIUMC = soa::Join<MyTracksIU, aod::McTrackLabels, aod::mcTPCTuneOnData>;

enum MatCorrType {
  None = 0,
  TGeo = 1,
  LUT = 2
};

enum TrackPropMode {
  kProper = 0,
  kFast = 1,
  kBoth = 2
};

enum V0DeduplicationMode {
  Pairwise = 0,       // old algorithm, V0 with best PCA + cosPA at the same time in a group wins
  GreedyMatching = 1, // new algorithm, best score wins (score combines PCA and cosPA -- see detailed implementation in PCMUtilities.cxx)
  GroupMatching = 2,  // new algorithm, best combination of leg-disjoint V0s wins, not the single best V0
  KeepAll = 3         // no deduplication, keep all V0s
};

// Struct needed for deduplication mode 1. Used to keep track  of v0 properties and a score based on pca and cosPA
struct V0CandidateHelper {
  int v0ID = -1;
  int colID = -1;
  int posID = -1;
  int eleID = -1;
  float cosPA = -1.f;
  float pca = -1.f;
  float score = -1.f;
  float mee = 0.f; // e+e- mass at the secondary vertex (GeV/c^2)

  V0CandidateHelper() = default;

  V0CandidateHelper(int v0, int col, int pos, int ele, float cpa, float p, float s, float m = 0.f)
    : v0ID(v0), colID(col), posID(pos), eleID(ele), cosPA(cpa), pca(p), score(s), mee(m)
  {
  }

  bool operator==(const V0CandidateHelper& other) const { return score == other.score; }
  bool operator<(const V0CandidateHelper& other) const { return score < other.score; }
  bool operator>(const V0CandidateHelper& other) const { return score > other.score; }
};

// (v0.globalIndex(), collision.globalIndex(), pos.globalIndex(), ele.globalIndex())
using CandKey = std::tuple<int64_t, int64_t, int64_t, int64_t>;

// Everything the truth diagnosis need,
struct DedupDiag {
  std::map<CandKey, std::vector<CandKey>> blockersByKey; // rejected candidate -> ALL candidates that took one of its legs
  std::vector<V0CandidateHelper> snapshot;               // ALL candidates, before deduplication
  std::vector<CandKey> stored;                           // the survivors
  // group matching only: rejected candidate -> (how many candidates the best solution CONTAINING
  // it would have cost, by how much its best solution was worse in total score).
  // deficit == 0 means an equally large alternative existed and the score alone decided.
  std::map<CandKey, std::pair<int, float>> lossMargin;
};

struct PhotonConversionBuilder {
  Produces<aod::V0PhotonsKF> v0photonskf;
  Produces<aod::V0Legs> v0legs;
  Produces<aod::V0LegsXYZ> v0legsXYZ;
  Produces<aod::V0LegsDeDxMC> v0legsDeDxMC;
  Produces<aod::V0PhotonsPhiVPsi> v0photonsphivpsi;
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
  Configurable<int> modeTrackPropagation{"modeTrackPropagation", 0, "0: use real track propagation, including material, 1: use fast approximation using only geometry, 2: Use real track propagation and make comparison to fast propagation (only for debugging and testing)"};
  Configurable<int> deduplicationMode{"deduplicationMode", 0, "0: Pairwise deduplication, 1: Based on Greedy matching (best score wins), 2: Based on Group matching (crossed pairs are both kept), 3: Keep all V0s, our default in the config is mode 0, however if a wrong configuration is used, the default will be mode 1 (Greedy matching) to avoid crashes"};
  Configurable<float> deduplicationScoreWeight{"deduplicationScoreWeight", 0.5f, "0.:only pca goes into the score, 1: only cosPA goes itno score, any number in between is a mix of pca and cosPA"};
  Configurable<bool> cfgDedupTruthMaps{"cfgDedupTruthMaps", true, "fill the fine (dEta, dPhi, q) truth maps in MC"};
  Configurable<int> dedupMaxGroupSize{"dedupMaxGroupSize", 12, "group matching (mode 2): conflict groups up to this size are solved exactly, larger ones greedily (capped at 16)"};

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
  Configurable<int> centTypePCMMl{"centTypePCMMl", 2, "Centrality type for 2D ML application: FT0M:0, FT0A:1, FT0C:2"};
  Configurable<std::vector<int>> cutDirPCMMl{"cutDirPCMMl", std::vector<int>{o2::analysis::em_cuts_ml::vecCutDir}, "Whether to reject score values greater or smaller than the threshold"};
  Configurable<std::vector<std::string>> namesInputFeatures{"namesInputFeatures", std::vector<std::string>{"feature1", "feature2"}, "Names of ML model input features"};
  Configurable<std::vector<std::string>> modelPathsCCDB{"modelPathsCCDB", std::vector<std::string>{"path_ccdb/BDT_PCM/"}, "Paths of models on CCDB"};
  Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"ModelHandler_onnx_PCM.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
  Configurable<std::vector<std::string>> labelsBinsPCMMl{"labelsBinsPCMMl", std::vector<std::string>{"bin 0", "bin 1"}, "Labels for bins"};
  Configurable<std::vector<std::string>> labelsCutScoresPCMMl{"labelsCutScoresPCMMl", std::vector<std::string>{o2::analysis::em_cuts_ml::labelsCutScore}, "Labels for cut scores"};
  Configurable<std::vector<double>> binsPtPCMMl{"binsPtPCMMl", std::vector<double>{0.0, +1e+10}, "pT bin limits for ML application"};
  Configurable<std::vector<double>> binsCentPCMMl{"binsCentPCMMl", std::vector<double>{0.0, 100.0}, "Centrality bin limits for ML application"};
  Configurable<std::vector<double>> cutsPCMMlFlat{"cutsPCMMlFlat", {0.5}, "Flattened ML cuts: [bin0_score0, bin0_score1, ..., binN_scoreM]"};

  Configurable<float> propV0LegsRadius{"propV0LegsRadius", 60.f, "Radius to which the V0 legs are propagated to calculate psipair and phiV"};

  o2::analysis::EmMlResponsePCM<float> emMlResponse;
  std::vector<float> outputML;
  V0PhotonCandidate v0photoncandidate;
  o2::ccdb::CcdbApi ccdbApi;

  int mRunNumber{};
  float d_bz{};
  float maxSnp{};  // max sine phi for propagation
  float maxStep{}; // max step size (cm) for propagation
  Service<o2::ccdb::BasicCCDBManager> ccdb{};
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
      {"V0/hPhiVPsiPair", "phiV vs. psi pair;#psi_{pair} (rad.);#phi_{V} (rad.)", {HistType::kTH2F, {{500, -o2::constants::math::PI, o2::constants::math::PI}, {500, 0.0f, o2::constants::math::TwoPI}}}},
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

    static constexpr std::array<const char*, 4> DedupNames = {"pairwise", "greedy matching", "group matching", "keep all"};
    if (deduplicationMode < 0 || deduplicationMode >= static_cast<int>(DedupNames.size())) {
      LOG(fatal) << "unknown deduplicationMode " << deduplicationMode.value;
    }
    LOGF(info, "photon-conversion-builder: deduplicationMode = %d (%s), score weight = %.2f",
         deduplicationMode.value, DedupNames[deduplicationMode.value], deduplicationScoreWeight.value);

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    switch (useMatCorrType) {
      case MatCorrType::TGeo:
        LOGF(info, "TGeo correction requested, loading geometry");
        if (!o2::base::GeometryManager::isGeometryLoaded()) {
          ccdb->get<TGeoManager>(geoPath);
        }
        matCorr = o2::base::Propagator::MatCorrType::USEMatCorrTGeo;
        break;
      case MatCorrType::LUT:
        LOGF(info, "LUT correction requested, loading LUT");
        lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(lutPath));
        matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
        break;
      default:
        LOGF(info, "no correction requested, loading LUT by default!");
        lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(lutPath));
        matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
        break;
    }

    if (applyPCMMl) {
      if (use2DBinning) {
        int binsNPt = static_cast<int>(binsPtPCMMl->size()) - 1;
        int binsNCent = static_cast<int>(binsCentPCMMl->size()) - 1;
        int binsN = binsNPt * binsNCent;
        if (binsN * static_cast<int>(cutDirPCMMl->size()) != static_cast<int>(cutsPCMMlFlat->size())) {
          LOG(fatal) << "Mismatch in number of bins and cuts provided for 2D ML application: binsN * cutDirPCMMl: " << binsN * static_cast<int>(cutDirPCMMl->size()) << " bins vs. cutsPCMMlFlat: " << cutsPCMMlFlat->size() << " cuts";
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
      if (nClassesPCMMl == o2::analysis::NClassesML::kTwo) {
        registry.add("V0/hBDTBackgroundScoreBeforeCutVsPt", "BDT background score before cut vs pT; pT (GeV/c); BDT background score", {HistType::kTH2F, {{1000, 0.0f, 20.0f}, {1000, 0.0f, 1.0f}}});
        registry.add("V0/hBDTBackgroundScoreAfterCutVsPt", "BDT background score after cut vs pT; pT (GeV/c); BDT background score", {HistType::kTH2F, {{1000, 0.0f, 20.0f}, {1000, 0.0f, 1.0f}}});
        registry.add("V0/hBDTSignalScoreBeforeCutVsPt", "BDT signal score before cut vs pT; pT (GeV/c); BDT signal score", {HistType::kTH2F, {{1000, 0.0f, 20.0f}, {1000, 0.0f, 1.0f}}});
        registry.add("V0/hBDTSignalScoreAfterCutVsPt", "BDT signal score after cut vs pT; pT (GeV/c); BDT signal score", {HistType::kTH2F, {{1000, 0.0f, 20.0f}, {1000, 0.0f, 1.0f}}});
      } else if (nClassesPCMMl == o2::analysis::NClassesML::kThree) {
        registry.add("V0/hBDTBackgroundScoreBeforeCutVsPt", "BDT background score before cut vs pT; pT (GeV/c); BDT background score", {HistType::kTH2F, {{1000, 0.0f, 20.0f}, {1000, 0.0f, 1.0f}}});
        registry.add("V0/hBDTBackgroundScoreAfterCutVsPt", "BDT background score after cut vs pT; pT (GeV/c); BDT background score", {HistType::kTH2F, {{1000, 0.0f, 20.0f}, {1000, 0.0f, 1.0f}}});
        registry.add("V0/hBDTPrimaryPhotonScoreBeforeCutVsPt", "BDT primary photon score before cut vs pT; pT (GeV/c); BDT primary photon score", {HistType::kTH2F, {{1000, 0.0f, 20.0f}, {1000, 0.0f, 1.0f}}});
        registry.add("V0/hBDTPrimaryPhotonScoreAfterCutVsPt", "BDT primary photon score after cut vs pT; pT (GeV/c); BDT primary photon score", {HistType::kTH2F, {{1000, 0.0f, 20.0f}, {1000, 0.0f, 1.0f}}});
        registry.add("V0/hBDTSecondaryPhotonScoreBeforeCutVsPt", "BDT secondary photon score before cut vs pT; pT (GeV/c); BDT secondary photon score", {HistType::kTH2F, {{1000, 0.0f, 20.0f}, {1000, 0.0f, 1.0f}}});
        registry.add("V0/hBDTSecondaryPhotonScoreAfterCutVsPt", "BDT secondary photon score after cut vs pT; pT (GeV/c); BDT secondary photon score", {HistType::kTH2F, {{1000, 0.0f, 20.0f}, {1000, 0.0f, 1.0f}}});
      } else {
        registry.add("V0/hBDTScoreBeforeCutVsPt", "BDT score before cut vs pT; pT (GeV/c); BDT score", {HistType::kTH2F, {{1000, 0.0f, 20.0f}, {1000, 0.0f, 1.0f}}});
        registry.add("V0/hBDTScoreAfterCutVsPt", "BDT score after cut vs pT; pT (GeV/c); BDT score", {HistType::kTH2F, {{1000, 0.0f, 20.0f}, {1000, 0.0f, 1.0f}}});
      }
    }

    // Compare proper propagation and fast geometrical propagation
    if (modeTrackPropagation == TrackPropMode::kBoth) {
      registry.add("V0Leg/hDCAxyPropagationCompare", "Comparison of DCA_{xy} propagation;DCA_{xy} (cm) proper propagation; DCA_{xy} (cm) geom. propagation", {HistType::kTH2F, {{200, -10., 10.}, {200, -10., 10.}}});
      registry.add("V0Leg/hDCAzPropagationCompare", "Comparison of DCA_{z} propagation;DCA_{z} (cm) proper propagation; DCA_{z} (cm) geom. propagation", {HistType::kTH2F, {{200, -10., 10.}, {200, -10., 10.}}});
      registry.add("V0/hPhivPropagationCompare", "Comparison of #phi_{v};#phi_{v} proper propagation; #phi_{v} geom. propagation", {HistType::kTH2F, {{100, 0., 1.6}, {100, 0., 1.6}}});
      registry.add("V0/hPsiPairPropagationCompare", "Comparison of #Psi_{pair};#Psi_{pair} proper propagation; #Psi_{pair} geom. propagation", {HistType::kTH2F, {{100, 0., 1.6}, {100, 0., 1.6}}});
    }

    // Make sure the deduplicationScoreWeight is between 0 and 1
    if (deduplicationScoreWeight > 1.f) { // o2-linter: disable=magic-number (score has to be below unity)
      LOG(warning) << "deduplicationScoreWeight is larger than unity which is not allowed";
    }
    if (deduplicationScoreWeight < 0.f) { // o2-linter: disable=magic-number (score has to be above zero)
      LOG(warning) << "deduplicationScoreWeight is smaller than zero which is not allowed";
    }

    if (doprocessMC) {
      const AxisSpec axQ{60, 0.f, 0.3f, "q_{inv}^{true} (GeV/c)"};
      const AxisSpec axDEta{80, -1.6f, 1.6f, "#Delta#eta_{#gamma#gamma}^{true}"};
      const AxisSpec axClass{3, -0.5f, 2.5f, "0 = true, 1 = cross-leg fake, 2 = other fake"};
      // One folder per pair fate, with the SAME histogram names inside. The post-processing then
      // only varies the folder, and everything belonging to one fate sits together.
      // NB: the "lost..." folders are diagnostic selections, NOT an exclusive decomposition -
      // lostByPartnerFake is a subset of lostByCrossFake, and a pair can appear in several of them.
      for (const auto& nm : {"before", "bothSurvive", "bothSurviveSameColl", "oneLost", "bothLost",
                             "lostByCrossFake", "lostByPartnerFake", "lostByOtherFake", "lostByTrue"}) {
        const std::string dir = std::string("MCDedup/Pairs/") + nm + "/";
        registry.add((dir + "hQ").c_str(), "truth photon pairs with >=1 true V0 candidate before deduplication;q_{inv}^{true} (GeV/c);pairs", kTH1F, {axQ}, true);
        registry.add((dir + "hDEta").c_str(), "truth photon pairs with >=1 true V0 candidate before deduplication;#Delta#eta^{true};pairs", kTH1F, {axDEta}, true);
      }
      registry.add("MCDedup/Candidates/hClassStage", "candidate composition;class;0 = before, 1 = stored", kTH2F, {axClass, {2, -0.5f, 1.5f}}, true);

      const AxisSpec axBlocker{4, -0.5f, 3.5f, "0 = true, 1 = cross-leg fake, 2 = other fake, 3 = no blocker (cut, not dedup)"};

      registry.add("MCDedup/Photons/hFate", "fate of the true photons;0 = alive, 1 = lost (blocked), 2 = lost (no blocker);photons", kTH1F, {{3, -0.5f, 2.5f}}, true);
      registry.add("MCDedup/Photons/hBlockerClass", "blockers of a killed true photon (ENTRIES = blockers, not photons);class of the blocker;blocker entries", kTH1F, {axBlocker}, true);
      // only meaningful for mode 0/1, where the smaller score always wins. Mode 2 compares SETS and
      // may deliberately keep the worse-scoring candidates, so it is not filled there.
      registry.add("MCDedup/Photons/hKillMargin", "mode 0/1 only;S_{victim} - S_{winner};kills", kTH1F, {{200, 0.f, 1.f}}, true);
      registry.add("MCDedup/Photons/hVictimVsWinnerPCA", "PCA of the two;PCA_{victim} (cm);PCA_{winner} (cm)", kTH2F, {{60, 0.f, 3.f}, {60, 0.f, 3.f}}, true);
      // mode 2 only: why exactly was this true photon not part of the winning set?
      // x = 0 -> an equally large solution containing it existed, only the total score decided.
      //          That is a pure ranking problem and the part that a better discriminant can win back.
      // x > 0 -> keeping it would have cost x candidates. Structural, no arbitration can repair it.
      registry.add("MCDedup/Photons/hLostMargin", "why the true photon lost (mode 2);cardinality deficit;#Sigma S(best set with it) - #Sigma S(winner)",
                   kTH2F, {{5, -0.5f, 4.5f}, {100, 0.f, 2.f}}, true);
      if (cfgDedupTruthMaps) {
        const AxisSpec axDEtaFine{640, -1.6f, 1.6f, "#Delta#eta^{true}"};
        const AxisSpec axDPhi{144, -o2::constants::math::PI, o2::constants::math::PI, "#Delta#varphi^{true} (rad)"};
        for (const auto& nm : {"before", "bothSurvive", "bothLost", "lostByPartnerFake"}) {
          // the fine map lands in the same folder as hQ and hDEta of that fate
          registry.add((std::string("MCDedup/Pairs/") + nm + "/hMap").c_str(), "truth photon pairs", kTHnSparseF, {axDEtaFine, axDPhi, axQ}, true);
        }
      }
    }
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    // In case override, don't proceed, please - no CCDB access required
    if (d_bz_input > -990) { // o2-linter: disable=magic-number (override value)
      d_bz = d_bz_input;
      o2::parameters::GRPMagField grpmag;
      if (std::fabs(d_bz) > 1e-5) {                   // o2-linter: disable=magic-number (override value)
        grpmag.setL3Current(30000.f / (d_bz / 5.0f)); // o2-linter: disable=magic-number (override value)
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

    if (useMatCorrType == 2) { // o2-linter: disable=magic-number (material budget correction)
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

  float cospaXY_KF(const KFParticle& kfp, const KFParticle& PV)
  {
    float lx = kfp.GetX() - PV.GetX(); // flight length X
    float ly = kfp.GetY() - PV.GetY(); // flight length Y

    float px = kfp.GetPx();
    float py = kfp.GetPy();
    float cospaXY = RecoDecay::dotProd(std::array{lx, ly}, std::array{px, py}) / (RecoDecay::sqrtSumOfSquares(lx, ly) * RecoDecay::sqrtSumOfSquares(px, py));
    if (cospaXY < -1.f) {
      return -1.f;
    }
    if (cospaXY > 1.f) {
      return 1.f;
    }
    return cospaXY;
  }

  float cospaRZ_KF(const KFParticle& kfp, const KFParticle& PV)
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
    }
    if (cospaRZ > 1.f) {
      return 1.f;
    }
    return cospaRZ;
  }

  template <bool isMC, typename TTrack, typename TShiftedTrack>
  void fillTrackTable(TTrack const& track, TShiftedTrack const& shiftedtrack, const float dcaXY, const float dcaZ)
  {
    v0legs(track.collisionId(), track.globalIndex(), track.sign(),
           shiftedtrack.GetPx(), shiftedtrack.GetPy(), shiftedtrack.GetPz(), dcaXY, dcaZ,
           track.tpcNClsFindable(), track.tpcNClsFindableMinusFound(), track.tpcNClsFindableMinusCrossedRows(), track.tpcNClsShared(),
           track.tpcChi2NCl(), track.tpcInnerParam(), track.tpcSignal(),
           track.tpcNSigmaEl(), track.tpcNSigmaPi(),
           track.itsClusterSizes(), track.itsChi2NCl(), track.detectorMap());

    v0legsXYZ(shiftedtrack.GetX(), shiftedtrack.GetY(), shiftedtrack.GetZ());

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
    std::array<float, 2> dcaInfo{};

    auto pTrack = getTrackParCov(pos);
    if (moveTPCTracks && isTPConlyTrack(pos) && !mVDriftMgr.moveTPCTrack<TBCs, TCollisions>(collision, pos, pTrack)) {
      LOGP(error, "failed correction for positive tpc track");
      return;
    }
    auto pTrackC = pTrack;
    o2::math_utils::Point3D<float> vtxPrim{
      collision.posX(),
      collision.posY(),
      collision.posZ()};
    if (modeTrackPropagation != TrackPropMode::kProper) {
      dcaInfo = CalculateDCAFast(pTrackC, vtxPrim, d_bz);
    }
    pTrackC.setPID(o2::track::PID::Electron);

    std::array<float, 2> dcaInfoFast = dcaInfo;
    if (modeTrackPropagation != TrackPropMode::kFast) {
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, pTrackC, 2.f, matCorr, &dcaInfo);
      if (filltable && modeTrackPropagation == TrackPropMode::kBoth) {
        registry.fill(HIST("V0Leg/hDCAxyPropagationCompare"), dcaInfo[0], dcaInfoFast[0]);
        registry.fill(HIST("V0Leg/hDCAzPropagationCompare"), dcaInfo[1], dcaInfoFast[1]);
      }
    }
    auto posdcaXY = dcaInfo[0];
    auto posdcaZ = dcaInfo[1];

    auto nTrack = getTrackParCov(ele);
    if (moveTPCTracks && isTPConlyTrack(ele) && !mVDriftMgr.moveTPCTrack<TBCs, TCollisions>(collision, ele, nTrack)) {
      LOGP(error, "failed correction for negative tpc track");
      return;
    }
    auto nTrackC = nTrack;
    if (modeTrackPropagation != TrackPropMode::kProper) {
      dcaInfo = CalculateDCAFast(nTrackC, vtxPrim, d_bz);
    }
    nTrackC.setPID(o2::track::PID::Electron);

    dcaInfoFast = dcaInfo;
    if (modeTrackPropagation != TrackPropMode::kFast) {
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, nTrackC, 2.f, matCorr, &dcaInfo);
      if (filltable && modeTrackPropagation == TrackPropMode::kBoth) {
        registry.fill(HIST("V0Leg/hDCAxyPropagationCompare"), dcaInfo[0], dcaInfoFast[0]);
        registry.fill(HIST("V0Leg/hDCAzPropagationCompare"), dcaInfo[1], dcaInfoFast[1]);
      }
    }
    auto eledcaXY = dcaInfo[0];
    auto eledcaZ = dcaInfo[1];

    if (std::fabs(posdcaXY) < dcapostopv || std::fabs(eledcaXY) < dcanegtopv) {
      return;
    }

    std::array<float, 3> xyz = {0.f, 0.f, 0.f};
    Vtx_recalculationParCov(o2::base::Propagator::Instance(), pTrack, nTrack, xyz, matCorr);
    float rxy_tmp = RecoDecay::sqrtSumOfSquares(xyz[0], xyz[1]);
    if (rxy_tmp > maxX + margin_r_tpc) {
      return;
    }
    if (rxy_tmp < std::fabs(xyz[2]) * std::tan(2 * std::atan(std::exp(-max_eta_v0))) - margin_z) {
      return; // RZ line cut
    }

    float phiv = 999.f;
    float psipair = 999.f;
    float phivFast = 999.f;
    float psipairFast = 999.f;
    float baseR = std::hypot(xyz[0], xyz[1]);
    // This method uses the track Helix instead of the full propagation.
    // Hence, it is only an approximation but much faster
    if (modeTrackPropagation != TrackPropMode::kProper) {

      o2::track::TrackAuxPar helixPosEle(nTrack, d_bz);
      o2::track::TrackAuxPar helixPosPos(pTrack, d_bz);

      float diffX = helixPosEle.xC - helixPosPos.xC;
      float diffY = helixPosEle.yC - helixPosPos.yC;
      auto phiHelix = RecoDecay::constrainAngle<float, float>(std::atan2(diffY, diffX) - o2::constants::math::PI / 2.);

      // Electron
      float arcLenghtEle = helixPosEle.rC * 0.9 > propV0LegsRadius ? std::asin(propV0LegsRadius / helixPosEle.rC) * helixPosEle.rC : o2::constants::math::PI / 2.2 * helixPosEle.rC; // This assumes that the photon momentum vector is a tangent of the circle // o2-linter: disable=magic-number (geometrical assumption for propagation)
      auto propTrackEle = getPropMomentumFromTrackHelix(arcLenghtEle, ele, helixPosEle, d_bz / 10., phiHelix - ele.phi());
      // Positron
      float arcLenghtPos = helixPosPos.rC * 0.9 > propV0LegsRadius ? std::asin(propV0LegsRadius / helixPosPos.rC) * helixPosPos.rC : o2::constants::math::PI / 2.2 * helixPosPos.rC; // This assumes that the photon momentum vector is a tangent of the circle // o2-linter: disable=magic-number (geometrical assumption for propagation)
      auto propTrackPos = getPropMomentumFromTrackHelix(arcLenghtPos, pos, helixPosPos, d_bz / 10., phiHelix - pos.phi());

      phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(propTrackPos[0], propTrackPos[1], propTrackPos[2], propTrackEle[0], propTrackEle[1], propTrackEle[2], pos.sign(), ele.sign(), d_bz);
      psipair = o2::aod::pwgem::dilepton::utils::pairutil::getPsiPair(propTrackPos[0], propTrackPos[1], propTrackPos[2], propTrackEle[0], propTrackEle[1], propTrackEle[2]);

      // Store values for later comparison
      phivFast = phiv;
      psipairFast = psipair;
    }
    // This uses the full propagation including material effects
    if (modeTrackPropagation != TrackPropMode::kFast) {
      std::array<float, 3> offsetsR = {propV0LegsRadius, 30.f, 10.f};
      bool pPropagatedSuccess = false;
      bool nPropagatedSuccess = false;
      auto pTrackProp = pTrack;
      auto nTrackProp = nTrack;
      for (const auto& offsetR : offsetsR) {
        pTrackProp = pTrack;
        pTrackProp.setPID(o2::track::PID::Electron);
        nTrackProp = nTrack;
        nTrackProp.setPID(o2::track::PID::Electron);

        o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, pTrackProp, 2.f, matCorr, &dcaInfo);
        o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, nTrackProp, 2.f, matCorr, &dcaInfo);

        pPropagatedSuccess = o2::base::Propagator::Instance()->propagateToR(pTrackProp, baseR + offsetR);
        nPropagatedSuccess = o2::base::Propagator::Instance()->propagateToR(nTrackProp, baseR + offsetR);

        if (pPropagatedSuccess && nPropagatedSuccess) {
          KFPTrack kfp_track_posProp = createKFPTrackFromTrackParCov(pTrackProp, pos.sign(), pos.tpcNClsFound(), pos.tpcChi2NCl());
          KFPTrack kfp_track_eleProp = createKFPTrackFromTrackParCov(nTrackProp, ele.sign(), ele.tpcNClsFound(), ele.tpcChi2NCl());
          phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(kfp_track_posProp.GetPx(), kfp_track_posProp.GetPy(), kfp_track_posProp.GetPz(), kfp_track_eleProp.GetPx(), kfp_track_eleProp.GetPy(), kfp_track_eleProp.GetPz(), pos.sign(), ele.sign(), d_bz);
          psipair = o2::aod::pwgem::dilepton::utils::pairutil::getPsiPair(kfp_track_posProp.GetPx(), kfp_track_posProp.GetPy(), kfp_track_posProp.GetPz(), kfp_track_eleProp.GetPx(), kfp_track_eleProp.GetPy(), kfp_track_eleProp.GetPz());
          break;
        }
        LOG(debug) << "Propagation to offset" << offsetR << " cm failed for " << (pPropagatedSuccess ? "negative" : "positive") << " track. Trying smaller offset.";
      }
      if (filltable && modeTrackPropagation == TrackPropMode::kBoth) {
        registry.fill(HIST("V0/hPhivPropagationCompare"), phiv, phivFast);
        registry.fill(HIST("V0/hPsiPairPropagationCompare"), psipair, psipairFast);
      }
    }

    if (phiv == 999.f || psipair == 999.f) { // o2-linter: disable=magic-number (nonsensical default value)
      LOG(debug) << "Propagation failed for all radii (" << propV0LegsRadius << ", 30, 10 cm). Using default values for phiv and psipair (999.f).";
    }

    KFPTrack kfp_track_pos = createKFPTrackFromTrackParCov(pTrack, pos.sign(), pos.tpcNClsFound(), pos.tpcChi2NCl());
    KFPTrack kfp_track_ele = createKFPTrackFromTrackParCov(nTrack, ele.sign(), ele.tpcNClsFound(), ele.tpcChi2NCl());
    KFParticle kfp_pos(kfp_track_pos, kPositron);
    KFParticle kfp_ele(kfp_track_ele, kElectron);
    std::array<const KFParticle*, 2> GammaDaughters{&kfp_pos, &kfp_ele};

    KFParticle gammaKF;
    gammaKF.SetConstructMethod(2);
    gammaKF.Construct(GammaDaughters.data(), 2);
    if (kfMassConstrain > -0.1) { // o2-linter: disable=magic-number (nonsensical default value)
      gammaKF.SetNonlinearMassConstraint(kfMassConstrain);
    }
    KFPVertex kfpVertex = createKFPVertexFromCollision(collision);
    KFParticle KFPV(kfpVertex);

    // Transport the gamma to the recalculated decay vertex
    KFParticle gammaKF_DecayVtx = gammaKF; // with respect to (0,0,0)
    gammaKF_DecayVtx.TransportToPoint(xyz.data());

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

    KFParticle kfp_pos_DecayVtx = kfp_pos;         // Don't set Primary Vertex
    KFParticle kfp_ele_DecayVtx = kfp_ele;         // Don't set Primary Vertex
    kfp_pos_DecayVtx.TransportToPoint(xyz.data()); // Don't set Primary Vertex
    kfp_ele_DecayVtx.TransportToPoint(xyz.data()); // Don't set Primary Vertex

    float cospaXYKF = cospaXY_KF(gammaKF_DecayVtx, KFPV);
    float cospaRZKF = cospaRZ_KF(gammaKF_DecayVtx, KFPV);
    auto centType = static_cast<CentType>(centTypePCMMl.value);
    v0photoncandidate.setPhotonCandidate(gammaKF_DecayVtx, gammaKF_PV, pos, kfp_pos_DecayVtx, ele, kfp_ele_DecayVtx, collision, cospaXYKF, cospaRZKF, cospaXYKF, psipair, phiv, centType, posdcaXY, eledcaXY, posdcaZ, eledcaZ);

    if (!ele.hasITS() && !pos.hasITS()) { // V0s with TPConly-TPConly
      if (max_r_itsmft_ss < rxy && rxy < maxX + margin_r_tpc) {
        if (v0photoncandidate.getPCA() > max_dcav0dau_tpc_inner_fc) {
          return;
        }
      } else {
        if (v0photoncandidate.getPCA() > max_dcav0dau_tpconly) {
          return;
        }
      }
    } else { // V0s with ITS hits
      if (rxy < max_r_req_its) {
        if (v0photoncandidate.getPCA() > max_dcav0dau_itsibss) {
          return;
        }
      } else {
        if (v0photoncandidate.getPCA() > max_dcav0dau_its) {
          return;
        }
      }
    }

    if (isITSonlyTrack(pos) && v0photoncandidate.getPosPt() > maxpt_itsonly) {
      return;
    }

    if (isITSonlyTrack(ele) && v0photoncandidate.getElePt() > maxpt_itsonly) {
      return;
    }

    if (v0photoncandidate.getChi2NDF() > 6e+3) { // protection for uint16. // o2-linter: disable=magic-number (protection for uint16)
      return;
    }

    if (std::fabs(v0photoncandidate.getDcaXYToPV()) > max_dcatopv_xy_v0 || std::fabs(v0photoncandidate.getDcaZToPV()) > max_dcatopv_z_v0) {
      return;
    }

    if (!checkAP(v0photoncandidate.getAlpha(), v0photoncandidate.getQt(), max_alpha_ap, max_qt_ap)) { // store only photon conversions
      return;
    }
    pca_map[std::make_tuple(v0.globalIndex(), collision.globalIndex(), pos.globalIndex(), ele.globalIndex())] = v0photoncandidate.getPCA();
    cospa_map[std::make_tuple(v0.globalIndex(), collision.globalIndex(), pos.globalIndex(), ele.globalIndex())] = v0photoncandidate.getCosPA();

    ROOT::Math::PxPyPzMVector vposDedup(kfp_pos_DecayVtx.GetPx(), kfp_pos_DecayVtx.GetPy(), kfp_pos_DecayVtx.GetPz(), o2::constants::physics::MassElectron);
    ROOT::Math::PxPyPzMVector veleDedup(kfp_ele_DecayVtx.GetPx(), kfp_ele_DecayVtx.GetPy(), kfp_ele_DecayVtx.GetPz(), o2::constants::physics::MassElectron);
    const auto meeSV = static_cast<float>((vposDedup + veleDedup).M());

    const auto score = getScoreV0(v0photoncandidate.getCosPA(), v0photoncandidate.getPCA(), deduplicationScoreWeight);
    V0CandidateHelper v0Helper(v0.globalIndex(), collision.globalIndex(), pos.globalIndex(), ele.globalIndex(),
                               v0photoncandidate.getCosPA(), v0photoncandidate.getPCA(), score, meeSV);
    vecV0Dedup.emplace_back(v0Helper);

    if (applyPCMMl) {
      bool isSelectedML = false;
      std::vector<float> mlInputFeatures = emMlResponse.getInputFeatures(v0photoncandidate, pos, ele);
      if (use2DBinning) {
        isSelectedML = emMlResponse.isSelectedMl(mlInputFeatures, v0photoncandidate.getPt(), v0photoncandidate.getCent(), outputML);
      } else {
        isSelectedML = emMlResponse.isSelectedMl(mlInputFeatures, v0photoncandidate.getPt(), outputML);
      }
      if (filltable) {
        if (nClassesPCMMl == o2::analysis::NClassesML::kTwo) {
          registry.fill(HIST("V0/hBDTBackgroundScoreBeforeCutVsPt"), v0photoncandidate.getPt(), outputML[0]);
          registry.fill(HIST("V0/hBDTSignalScoreBeforeCutVsPt"), v0photoncandidate.getPt(), outputML[1]);
        } else if (nClassesPCMMl == o2::analysis::NClassesML::kThree) {
          registry.fill(HIST("V0/hBDTPrimaryPhotonScoreBeforeCutVsPt"), v0photoncandidate.getPt(), outputML[0]);
          registry.fill(HIST("V0/hBDTSecondaryPhotonScoreBeforeCutVsPt"), v0photoncandidate.getPt(), outputML[1]);
          registry.fill(HIST("V0/hBDTBackgroundScoreBeforeCutVsPt"), v0photoncandidate.getPt(), outputML[2]);
        } else {
          registry.fill(HIST("V0/hBDTScoreBeforeCutVsPt"), v0photoncandidate.getPt(), outputML[0]);
        }
      }
      if (!isSelectedML) {
        return;
      }
      if (filltable) {
        if (nClassesPCMMl == o2::analysis::NClassesML::kTwo) {
          registry.fill(HIST("V0/hBDTBackgroundScoreAfterCutVsPt"), v0photoncandidate.getPt(), outputML[0]);
          registry.fill(HIST("V0/hBDTSignalScoreAfterCutVsPt"), v0photoncandidate.getPt(), outputML[1]);
        } else if (nClassesPCMMl == o2::analysis::NClassesML::kThree) {
          registry.fill(HIST("V0/hBDTPrimaryPhotonScoreAfterCutVsPt"), v0photoncandidate.getPt(), outputML[0]);
          registry.fill(HIST("V0/hBDTSecondaryPhotonScoreAfterCutVsPt"), v0photoncandidate.getPt(), outputML[1]);
          registry.fill(HIST("V0/hBDTBackgroundScoreAfterCutVsPt"), v0photoncandidate.getPt(), outputML[2]);
        } else {
          registry.fill(HIST("V0/hBDTScoreAfterCutVsPt"), v0photoncandidate.getPt(), outputML[0]);
        }
      }
    }

    if (filltable) {
      registry.fill(HIST("V0/hAP"), v0photoncandidate.getAlpha(), v0photoncandidate.getQt());
      registry.fill(HIST("V0/hConversionPointXY"), gammaKF_DecayVtx.GetX(), gammaKF_DecayVtx.GetY());
      registry.fill(HIST("V0/hConversionPointRZ"), gammaKF_DecayVtx.GetZ(), rxy);
      registry.fill(HIST("V0/hPt"), v0photoncandidate.getPt());
      registry.fill(HIST("V0/hEtaPhi"), v0phi, v0eta);
      registry.fill(HIST("V0/hCosPA"), v0photoncandidate.getCosPA());
      registry.fill(HIST("V0/hCosPA_Rxy"), rxy, v0photoncandidate.getCosPA());
      registry.fill(HIST("V0/hPCA"), v0photoncandidate.getPCA());
      registry.fill(HIST("V0/hPCA_CosPA"), v0photoncandidate.getCosPA(), v0photoncandidate.getPCA());
      registry.fill(HIST("V0/hPCA_Rxy"), rxy, v0photoncandidate.getPCA());
      registry.fill(HIST("V0/hDCAxyz"), v0photoncandidate.getDcaXYToPV(), v0photoncandidate.getDcaZToPV());
      registry.fill(HIST("V0/hPCA_diffX"), v0photoncandidate.getPCA(), std::min(pTrack.getX(), nTrack.getX()) - rxy); // trackiu.x() - rxy should be positive
      registry.fill(HIST("V0/hPhiVPsiPair"), v0photoncandidate.getPsiPair(), v0photoncandidate.getPhiV());

      // LOGF(info, "cospa_kf = %f, cospaXY_kf = %f, cospaRZ_kf = %f", cospa_kf, cospaXY_kf, cospaRZ_kf);
      registry.fill(HIST("V0/hCosPAXY_Rxy"), rxy, cospaXYKF);
      registry.fill(HIST("V0/hCosPARZ_Rxy"), rxy, cospaRZKF);

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
                  v0_sv.M(), v0photoncandidate.getDcaXYToPV(), v0photoncandidate.getDcaZToPV(),
                  cospa_kf, cospaXYKF, cospaRZKF,
                  v0photoncandidate.getPCA(), v0photoncandidate.getAlpha(), v0photoncandidate.getQt(), v0photoncandidate.getChi2NDF());
      v0photonsphivpsi(v0photoncandidate.getPhiV(), v0photoncandidate.getPsiPair());

      // v0photonskfcov(gammaKF_PV.GetCovariance(9), gammaKF_PV.GetCovariance(14), gammaKF_PV.GetCovariance(20), gammaKF_PV.GetCovariance(13), gammaKF_PV.GetCovariance(19), gammaKF_PV.GetCovariance(18));

      fillTrackTable<isMC>(pos, kfp_pos_DecayVtx, posdcaXY, posdcaZ); // positive leg first
      fillTrackTable<isMC>(ele, kfp_ele_DecayVtx, eledcaXY, eledcaZ); // negative leg second
    } // end of fill table
  }

  Preslice<aod::V0s> perCollision = o2::aod::v0::collisionId;
  std::map<std::tuple<int64_t, int64_t, int64_t, int64_t>, float> pca_map;      // (v0.globalIndex(), collision.globalIndex(), pos.globalIndex(), ele.globalIndex()) -> pca
  std::map<std::tuple<int64_t, int64_t, int64_t, int64_t>, float> cospa_map;    // (v0.globalIndex(), collision.globalIndex(), pos.globalIndex(), ele.globalIndex()) -> cospa
  std::vector<V0CandidateHelper> vecV0Dedup;                                    // vector with all V0Candidates that is used to sort them by score (see struct V0CandidateHelper for more details of content)
  std::vector<std::pair<int64_t, int64_t>> stored_v0Ids;                        // (pos.globalIndex(), ele.globalIndex())
  std::vector<std::tuple<int64_t, int64_t, int64_t, int64_t>> stored_fullv0Ids; // (v0.globalIndex(), collision.globalIndex(), pos.globalIndex(), ele.globalIndex())
  std::unordered_map<int64_t, int> nv0_map;                                     // map collisionId -> nv0

  template <bool isMC, bool isTriggerAnalysis, bool enableFilter, typename TCollisions, typename TV0s, typename TTracks, typename TBCs>
  void build(TCollisions const& collisions, TV0s const& v0s, TTracks const& tracks, TBCs const&, DedupDiag* diag = nullptr)
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

      // if constexpr (isTriggerAnalysis) {
      //   if (collision.triggerMask_raw() == 0) {
      //     continue;
      //   }
      // }

      nv0_map[collision.globalIndex()] = 0;

      const auto& bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      registry.fill(HIST("hCollisionCounter"), 1);

      updateCCDB(bc); // delay update until is needed

      vecV0Dedup.reserve(30000); // rough estimate for number of V0s per DF
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
    if (deduplicationMode == V0DeduplicationMode::Pairwise) {
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

          if (collisionId != collisionId_tmp && eleId == eleId_tmp && posId == posId_tmp && cospa < cospa_tmp) {
            if (diag != nullptr) {
              diag->blockersByKey[key].push_back(key_tmp);
            }
            is_most_aligned_v0 = false;
            break;
          }
          if ((eleId == eleId_tmp || posId == posId_tmp) && v0pca > v0pca_tmp) {
            if (diag != nullptr) {
              diag->blockersByKey[key].push_back(key_tmp);
            }
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

    } else if (deduplicationMode == V0DeduplicationMode::KeepAll) {
      // Keep-all: Every candidate that survived the quality cuts is stored, INCLUDING the
      // collision duplicates.
      stored_v0Ids.clear();
      stored_fullv0Ids.clear();
      nv0_map.clear();
      for (const auto& v0Cand : vecV0Dedup) {
        stored_v0Ids.emplace_back(v0Cand.posID, v0Cand.eleID);
        stored_fullv0Ids.emplace_back(v0Cand.v0ID, v0Cand.colID, v0Cand.posID, v0Cand.eleID);
        nv0_map[v0Cand.colID]++;
      }
      std::sort(stored_fullv0Ids.begin(), stored_fullv0Ids.end(),
                [](const auto& a, const auto& b) {
                  return std::get<1>(a) < std::get<1>(b);
                });
    } else if (deduplicationMode == V0DeduplicationMode::GroupMatching) {
      const int maxEnum = std::min(dedupMaxGroupSize.value, 16); // o2-linter: disable=magic-number (subset check of photons)
      const int nCand = static_cast<int>(vecV0Dedup.size());
      stored_v0Ids.clear();
      stored_fullv0Ids.clear();
      nv0_map.clear();
      int nGreedyFallback = 0;

      // ---- conflict groups: candidates connected through a shared leg -----
      std::vector<int> parent(nCand);
      for (int i = 0; i < nCand; ++i) {
        parent[i] = i;
      }
      auto findRoot = [&parent](int x) {
        while (parent[x] != x) {
          parent[x] = parent[parent[x]]; // path halving
          x = parent[x];
        }
        return x;
      };
      std::unordered_map<int, int> firstCandOfTrack;
      for (int i = 0; i < nCand; ++i) {
        for (const int& trackId : {vecV0Dedup[i].posID, vecV0Dedup[i].eleID}) {
          auto [it, inserted] = firstCandOfTrack.try_emplace(trackId, i);
          if (!inserted) {
            parent[findRoot(i)] = findRoot(it->second);
          }
        }
      }
      std::map<int, std::vector<int>> groups;
      for (int i = 0; i < nCand; ++i) {
        groups[findRoot(i)].push_back(i);
      }

      // ---- decide each group on its own ----------------------------------
      for (auto& [root, members] : groups) { // o2-linter: disable=const-ref-in-for-loop (members is sorted in place below)
        const int nMembers = static_cast<int>(members.size());
        std::vector<char> selected(nMembers, 0);

        // Fix the order inside the group: best score first, ties by v0ID. This makes the
        // result independent of the order in which the candidates ended up in vecV0Dedup.
        std::sort(members.begin(), members.end(), [this](int a, int b) {
          if (vecV0Dedup[a].score != vecV0Dedup[b].score) {
            return vecV0Dedup[a].score < vecV0Dedup[b].score;
          }
          return vecV0Dedup[a].v0ID < vecV0Dedup[b].v0ID;
        });

        if (nMembers == 1) {
          selected[0] = 1; // no conflict -> always kept
        } else if (nMembers <= maxEnum) {
          // Exact. Instead of comparing track IDs for every subset, each leg of the group gets
          // a local bit number; "share a leg" is then a single bitwise AND.
          std::unordered_map<int, int> legBit;
          std::vector<uint64_t> legMask(nMembers, 0);
          for (int k = 0; k < nMembers; ++k) {
            for (const int& trackId : {vecV0Dedup[members[k]].posID, vecV0Dedup[members[k]].eleID}) {
              const int bit = legBit.try_emplace(trackId, static_cast<int>(legBit.size())).first->second;
              legMask[k] |= (1ull << bit); // at most 2 * 16 = 32 legs, fits into 64 bit
            }
          }
          uint32_t bestMask = 0;
          int bestCount = 0;
          float bestScore = 0.f;
          // per candidate: the best solution that CONTAINS it. Comparing that against the overall
          // best tells afterwards whether a rejected candidate lost on cardinality (structural,
          // no algorithm can help) or on total score in a tie (a ranking problem, i.e. fixable).
          std::vector<int> bestWithCount(nMembers, -1);
          std::vector<float> bestWithScore(nMembers, 0.f);
          for (uint32_t mask = 1; mask < (1u << nMembers); ++mask) {
            uint64_t usedLegsLocal = 0;
            float scoreSum = 0.f;
            int count = 0;
            bool disjoint = true;
            for (int k = 0; k < nMembers; ++k) {
              if ((mask & (1u << k)) == 0) {
                continue;
              }
              if ((usedLegsLocal & legMask[k]) != 0) {
                disjoint = false;
                break;
              }
              usedLegsLocal |= legMask[k];
              scoreSum += vecV0Dedup[members[k]].score;
              ++count;
            }
            if (!disjoint) {
              continue;
            }
            // lexicographic: cardinality first, then total score. Strict comparisons mean the
            // lowest mask wins a true tie, and since members is sorted by score those are the
            // better candidates.
            if (count > bestCount || (count == bestCount && scoreSum < bestScore)) {
              bestCount = count;
              bestScore = scoreSum;
              bestMask = mask;
            }
            if (diag != nullptr) {
              for (int k = 0; k < nMembers; ++k) {
                if ((mask & (1u << k)) == 0) {
                  continue;
                }
                if (count > bestWithCount[k] || (count == bestWithCount[k] && scoreSum < bestWithScore[k])) {
                  bestWithCount[k] = count;
                  bestWithScore[k] = scoreSum;
                }
              }
            }
          }
          for (int k = 0; k < nMembers; ++k) {
            selected[k] = ((bestMask & (1u << k)) != 0) ? 1 : 0;
            if (diag != nullptr && selected[k] == 0 && bestWithCount[k] >= 0) {
              const auto& cand = vecV0Dedup[members[k]];
              diag->lossMargin[std::make_tuple(static_cast<int64_t>(cand.v0ID), static_cast<int64_t>(cand.colID),
                                               static_cast<int64_t>(cand.posID), static_cast<int64_t>(cand.eleID))] =
                std::make_pair(bestCount - bestWithCount[k], bestWithScore[k] - bestScore);
            }
          }
        } else {
          // Too large to enumerate: greedy by score (members is already sorted that way). This
          // gives a maximal, not necessarily a maximum matching, so it is counted below.
          ++nGreedyFallback;
          std::vector<int> usedTracks;
          usedTracks.reserve(2 * nMembers);
          for (int k = 0; k < nMembers; ++k) {
            const auto& cand = vecV0Dedup[members[k]];
            if (std::find(usedTracks.begin(), usedTracks.end(), cand.posID) != usedTracks.end() ||
                std::find(usedTracks.begin(), usedTracks.end(), cand.eleID) != usedTracks.end()) {
              continue;
            }
            usedTracks.push_back(cand.posID);
            usedTracks.push_back(cand.eleID);
            selected[k] = 1;
          }
        }

        // ---- store the winners, book the blockers of the losers ----------
        for (int k = 0; k < nMembers; ++k) {
          const auto& cand = vecV0Dedup[members[k]];
          if (selected[k] != 0) {
            stored_v0Ids.emplace_back(cand.posID, cand.eleID);
            stored_fullv0Ids.emplace_back(cand.v0ID, cand.colID, cand.posID, cand.eleID);
            nv0_map[cand.colID]++;
            continue;
          }
          if (diag == nullptr) {
            continue;
          }
          // Record ALL winners that took a leg of this candidate. A cross-leg fake is typically
          // blocked by TWO different true photons (one per leg); storing only the first would
          // make the later MC cause-of-loss attribution depend on the loop order.
          const CandKey victimKey = std::make_tuple(static_cast<int64_t>(cand.v0ID), static_cast<int64_t>(cand.colID),
                                                    static_cast<int64_t>(cand.posID), static_cast<int64_t>(cand.eleID));
          for (int j = 0; j < nMembers; ++j) {
            if (selected[j] == 0) {
              continue;
            }
            const auto& winner = vecV0Dedup[members[j]];
            if (winner.posID == cand.posID || winner.eleID == cand.eleID) {
              diag->blockersByKey[victimKey].emplace_back(static_cast<int64_t>(winner.v0ID), static_cast<int64_t>(winner.colID),
                                                          static_cast<int64_t>(winner.posID), static_cast<int64_t>(winner.eleID));
            }
          }
        }
      }

      if (nGreedyFallback > 0) {
        LOGF(warning, "group matching: %d conflict group(s) larger than dedupMaxGroupSize = %d were solved greedily", nGreedyFallback, maxEnum);
      }

      // the table expects candidates ordered by collision
      std::sort(stored_fullv0Ids.begin(), stored_fullv0Ids.end(),
                [](const auto& a, const auto& b) {
                  return std::get<1>(a) < std::get<1>(b);
                });

    } else {
      // Sort best candidates first, depending on score
      std::sort(vecV0Dedup.begin(), vecV0Dedup.end());
      // container to keep track of which electron and positron tracks have already been used
      std::vector<bool> usedLegs(tracks.size(), false);
      std::unordered_map<int, CandKey> ownerOfLeg; // diagnostics only
      stored_v0Ids.clear();
      stored_fullv0Ids.clear();
      nv0_map.clear();
      for (const auto& v0Cand : vecV0Dedup) {
        const CandKey candKey = std::make_tuple(static_cast<int64_t>(v0Cand.v0ID), static_cast<int64_t>(v0Cand.colID),
                                                static_cast<int64_t>(v0Cand.posID), static_cast<int64_t>(v0Cand.eleID));
        if (usedLegs[v0Cand.posID] || usedLegs[v0Cand.eleID]) {
          if (diag != nullptr) {
            // both legs, not just the first one found: a cross-leg fake is usually blocked by
            // two different winners, and that is exactly the case being measured here
            auto& blockers = diag->blockersByKey[candKey];
            for (const int& legId : {v0Cand.posID, v0Cand.eleID}) {
              const auto itOwner = ownerOfLeg.find(legId);
              if (itOwner == ownerOfLeg.end()) {
                continue;
              }
              if (std::find(blockers.begin(), blockers.end(), itOwner->second) == blockers.end()) {
                blockers.push_back(itOwner->second); // the same winner may hold both legs
              }
            }
          }
          continue;
        }
        usedLegs[v0Cand.posID] = true;
        usedLegs[v0Cand.eleID] = true;
        if (diag != nullptr) {
          ownerOfLeg[v0Cand.posID] = candKey;
          ownerOfLeg[v0Cand.eleID] = candKey;
        }
        stored_v0Ids.emplace_back(v0Cand.posID, v0Cand.eleID);
        stored_fullv0Ids.emplace_back(v0Cand.v0ID, v0Cand.colID, v0Cand.posID, v0Cand.eleID);
        nv0_map[v0Cand.colID]++;
      }

      // sort accepted candidates by collision ID to be compatible with table
      std::sort(stored_fullv0Ids.begin(), stored_fullv0Ids.end(),
                [](const auto& a, const auto& b) {
                  return std::get<1>(a) < std::get<1>(b);
                });

      // End of deduplication
    }

    if (diag != nullptr) {
      diag->snapshot = vecV0Dedup;
      diag->stored = stored_fullv0Ids;
    }

    for (const auto& fullv0Id : stored_fullv0Ids) {
      auto v0Id = std::get<0>(fullv0Id);
      // auto collisionId = std::get<1>(fullv0Id);
      // auto posId = std::get<2>(fullv0Id);
      // auto eleId = std::get<3>(fullv0Id);
      // LOGF(info, "!accept! | collision id = %d | v0id = %d , posid = %d , eleid = %d", collisionId, v0Id, posId, eleId);

      auto v0 = v0s.rawIteratorAt(v0Id);
      if constexpr (enableFilter) {
        auto collision_tmp = v0.template collision_as<TCollisions>();                                               // collision where this v0 belongs.
        if (!(collision_tmp.neeuls() >= 1 || collision_tmp.neeuls() + nv0_map[collision_tmp.globalIndex()] >= 2)) { // o2-linter: disable=magic-number (nonsensical default value)
          continue;
        }
        // LOGF(info, "collision_tmp.globalIndex() = %d, collision_tmp.neeuls() = %d, nv0_map = %d", collision_tmp.globalIndex(), collision_tmp.neeuls(), nv0_map[collision_tmp.globalIndex()]);
      }

      fillV0Table<isMC, TBCs, TCollisions, TTracks>(v0, true);
    } // end of fullv0Id loop

    // for (const auto& collision : collisions) {
    //   if constexpr (isMC) {
    //     if (!collision.has_mcCollision()) {
    //       continue;
    //     }
    //   }
    //   // events_ngpcm(nv0_map[collision.globalIndex()]);
    // } // end of collision loop

    pca_map.clear();
    cospa_map.clear();
    nv0_map.clear();
    stored_v0Ids.clear();
    stored_v0Ids.shrink_to_fit();
    stored_fullv0Ids.clear();
    stored_fullv0Ids.shrink_to_fit();
    vecV0Dedup.clear();
  } // end of build

  Filter v0Filter = o2::aod::v0::v0Type > (uint8_t)0;
  using FilteredV0s = soa::Filtered<aod::V0s>;

  enum PairFate {
    kBefore = 0,
    kBothSurvive,
    kOneLost,
    kBothLost,
    // NB: the four "lostBy..." entries are diagnostic selections, NOT an exclusive decomposition.
    // kLostByPartnerFake is a SUBSET of kLostByCrossFake, and cross + other + true do not add up
    // to oneLost + bothLost (a pair can lose two photons to two different blockers).
    kLostByCrossFake,
    kLostByPartnerFake,
    kLostByOtherFake,
    kLostByTrue
  };

  void fillTruthPair(int fate, float q, float dEta)
  {
    switch (fate) {
      case kBefore:
        registry.fill(HIST("MCDedup/Pairs/before/hQ"), q);
        registry.fill(HIST("MCDedup/Pairs/before/hDEta"), dEta);
        break;
      case kBothSurvive:
        registry.fill(HIST("MCDedup/Pairs/bothSurvive/hQ"), q);
        registry.fill(HIST("MCDedup/Pairs/bothSurvive/hDEta"), dEta);
        break;
      case kOneLost:
        registry.fill(HIST("MCDedup/Pairs/oneLost/hQ"), q);
        registry.fill(HIST("MCDedup/Pairs/oneLost/hDEta"), dEta);
        break;
      case kBothLost:
        registry.fill(HIST("MCDedup/Pairs/bothLost/hQ"), q);
        registry.fill(HIST("MCDedup/Pairs/bothLost/hDEta"), dEta);
        break;
      case kLostByCrossFake:
        registry.fill(HIST("MCDedup/Pairs/lostByCrossFake/hQ"), q);
        registry.fill(HIST("MCDedup/Pairs/lostByCrossFake/hDEta"), dEta);
        break;
      case kLostByPartnerFake:
        registry.fill(HIST("MCDedup/Pairs/lostByPartnerFake/hQ"), q);
        registry.fill(HIST("MCDedup/Pairs/lostByPartnerFake/hDEta"), dEta);
        break;
      case kLostByOtherFake:
        registry.fill(HIST("MCDedup/Pairs/lostByOtherFake/hQ"), q);
        registry.fill(HIST("MCDedup/Pairs/lostByOtherFake/hDEta"), dEta);
        break;
      case kLostByTrue:
        registry.fill(HIST("MCDedup/Pairs/lostByTrue/hQ"), q);
        registry.fill(HIST("MCDedup/Pairs/lostByTrue/hDEta"), dEta);
        break;
      default:
        break;
    }
  }

  void fillTruthPairMap(int fate, float dEta, float dPhi, float q)
  {
    if (!cfgDedupTruthMaps) {
      return;
    }
    switch (fate) {
      case kBefore:
        registry.fill(HIST("MCDedup/Pairs/before/hMap"), dEta, dPhi, q);
        break;
      case kBothSurvive:
        registry.fill(HIST("MCDedup/Pairs/bothSurvive/hMap"), dEta, dPhi, q);
        break;
      case kBothLost:
        registry.fill(HIST("MCDedup/Pairs/bothLost/hMap"), dEta, dPhi, q);
        break;
      case kLostByPartnerFake:
        registry.fill(HIST("MCDedup/Pairs/lostByPartnerFake/hMap"), dEta, dPhi, q);
        break;
      default:
        break;
    }
  }

  template <o2::soa::is_table TTracks, o2::soa::is_table TMCParticles>
  void fillDedupTruthDiagnostics(TTracks const& tracks, TMCParticles const& mcparticles, DedupDiag const& diag)
  {
    const int nCand = static_cast<int>(diag.snapshot.size());
    if (nCand == 0) {
      return;
    }
    const std::set<CandKey> storedSet(diag.stored.begin(), diag.stored.end());

    std::vector<CandKey> key(nCand);
    std::vector<uint8_t> cls(nCand, static_cast<uint8_t>(kV0OtherFake));
    std::vector<int64_t> motherPos(nCand, -1), motherEle(nCand, -1);
    std::map<CandKey, int> indexOfKey;
    for (int i = 0; i < nCand; ++i) {
      const auto& c = diag.snapshot[i];
      key[i] = std::make_tuple(static_cast<int64_t>(c.v0ID), static_cast<int64_t>(c.colID),
                               static_cast<int64_t>(c.posID), static_cast<int64_t>(c.eleID));
      indexOfKey[key[i]] = i;
      cls[i] = static_cast<uint8_t>(classifyV0Truth(tracks.rawIteratorAt(c.posID), tracks.rawIteratorAt(c.eleID),
                                                    mcparticles, motherPos[i], motherEle[i]));
      registry.fill(HIST("MCDedup/Candidates/hClassStage"), static_cast<float>(cls[i]), 0.f);
      if (storedSet.contains(key[i]) > 0) {
        registry.fill(HIST("MCDedup/Candidates/hClassStage"), static_cast<float>(cls[i]), 1.f);
      }
    }

    // ---- per true photon: which photon passed and which one was removed
    // IMPORTANT: a photon is lost only when ALL of its candidates died, so all of them are asked
    // for blockers, not just the best one - different candidates have different opponents.
    struct PhotonFate {
      bool alive = false;
      int best = -1;                // best-score candidate, only used for the margin plots
      std::vector<int> cands;       // all true candidates of this photon
      std::set<int64_t> aliveColls; // collisions with a surviving candidate
    };
    std::map<int64_t, PhotonFate> fate;
    for (int i = 0; i < nCand; ++i) {
      if (cls[i] != kV0True) {
        continue;
      }
      auto& f = fate[motherPos[i]];
      f.cands.push_back(i);
      if (storedSet.contains(key[i]) > 0) {
        f.alive = true;
        f.aliveColls.insert(std::get<1>(key[i]));
      }
      if (f.best < 0 || diag.snapshot[i].score < diag.snapshot[f.best].score) {
        f.best = i;
      }
    }

    // motherId -> indices of ALL candidates that took a leg from this photon. Empty means the
    // photon disappeared without anybody claiming a leg: then the mode cut, it did not dedup.
    std::map<int64_t, std::vector<int>> blockersOfPhoton;
    for (const auto& [mid, f] : fate) {
      if (f.alive) {
        registry.fill(HIST("MCDedup/Photons/hFate"), 0.f);
        continue;
      }
      auto& blockers = blockersOfPhoton[mid];
      for (const int& ci : f.cands) {
        const auto itBlocker = diag.blockersByKey.find(key[ci]);
        if (itBlocker == diag.blockersByKey.end()) {
          continue;
        }
        for (const auto& bKey : itBlocker->second) {
          const auto itIndex = indexOfKey.find(bKey);
          if (itIndex == indexOfKey.end()) {
            continue;
          }
          if (std::find(blockers.begin(), blockers.end(), itIndex->second) == blockers.end()) {
            blockers.push_back(itIndex->second);
          }
        }
      }
      // Why was it not in the winning set? Take the most favourable of its candidates: the one
      // that could have been kept at the smallest cost in cardinality.
      int bestDeficit = -1;
      float bestExcess = 0.f;
      for (const int& ci : f.cands) {
        const auto itMargin = diag.lossMargin.find(key[ci]);
        if (itMargin == diag.lossMargin.end()) {
          continue;
        }
        if (bestDeficit < 0 || itMargin->second.first < bestDeficit ||
            (itMargin->second.first == bestDeficit && itMargin->second.second < bestExcess)) {
          bestDeficit = itMargin->second.first;
          bestExcess = itMargin->second.second;
        }
      }
      if (bestDeficit >= 0) {
        registry.fill(HIST("MCDedup/Photons/hLostMargin"), static_cast<float>(bestDeficit), bestExcess);
      }

      // ONE entry per photon, so the photon-level survival rate is comparable between modes.
      // hVictimBlockerClass below counts one entry per BLOCKER and a victim can have two.
      registry.fill(HIST("MCDedup/Photons/hFate"), blockers.empty() ? 2.f : 1.f);

      if (blockers.empty()) {
        registry.fill(HIST("MCDedup/Photons/hBlockerClass"), 3.f); // o2-linter: disable=magic-number (bin 3 = no blocker)
        continue;
      }
      for (const int& w : blockers) {
        registry.fill(HIST("MCDedup/Photons/hBlockerClass"), static_cast<float>(cls[w]));
      }
      // margin and PCA comparison relate the best candidate to its blocker - that is only the
      // quantity that decided in mode 0/1. Mode 2 compares sets, so it is skipped there.
      if (f.best >= 0 && deduplicationMode.value < static_cast<int>(V0DeduplicationMode::GroupMatching)) {
        const int w = blockers.front();
        registry.fill(HIST("MCDedup/Photons/hKillMargin"), diag.snapshot[f.best].score - diag.snapshot[w].score);
        registry.fill(HIST("MCDedup/Photons/hVictimVsWinnerPCA"), diag.snapshot[f.best].pca, diag.snapshot[w].pca);
      }
    }

    // ---- truth pairs, grouped by MC collision ----------------------------
    std::map<int, std::vector<int64_t>> mothersByMcColl;
    for (const auto& [mid, f] : fate) {
      mothersByMcColl[mcparticles.iteratorAt(mid).mcCollisionId()].push_back(mid);
    }
    constexpr float MaxQTruth = 0.3f;
    for (const auto& [mcCol, mids] : mothersByMcColl) {
      for (size_t a = 0; a < mids.size(); ++a) {
        for (size_t b = a + 1; b < mids.size(); ++b) {
          const auto g1 = mcparticles.iteratorAt(mids[a]);
          const auto g2 = mcparticles.iteratorAt(mids[b]);
          const float e1 = std::hypot(g1.px(), g1.py(), g1.pz());
          const float e2 = std::hypot(g2.px(), g2.py(), g2.pz());
          const float qTrue = std::sqrt(std::max(0.f, 2.f * (e1 * e2 - g1.px() * g2.px() - g1.py() * g2.py() - g1.pz() * g2.pz())));
          if (qTrue > MaxQTruth) {
            continue;
          }
          const float dEta = g1.eta() - g2.eta();
          const float dPhi = RecoDecay::constrainAngle(static_cast<float>(g1.phi() - g2.phi()), -o2::constants::math::PI);
          fillTruthPair(kBefore, qTrue, dEta);
          fillTruthPairMap(kBefore, dEta, dPhi, qTrue);

          const auto& f1 = fate.at(mids[a]);
          const auto& f2 = fate.at(mids[b]);
          const int nLost = (f1.alive ? 0 : 1) + (f2.alive ? 0 : 1);
          if (nLost == 0) {
            fillTruthPair(kBothSurvive, qTrue, dEta);
            fillTruthPairMap(kBothSurvive, dEta, dPhi, qTrue);
            bool sameColl = false;
            for (const auto& c1 : f1.aliveColls) {
              sameColl = sameColl || (f2.aliveColls.count(c1) > 0);
            }
            if (sameColl) {
              registry.fill(HIST("MCDedup/Pairs/bothSurviveSameColl/hQ"), qTrue);
            }
            continue;
          }
          if (nLost == 1) {
            fillTruthPair(kOneLost, qTrue, dEta);
          } else {
            fillTruthPair(kBothLost, qTrue, dEta);
            fillTruthPairMap(kBothLost, dEta, dPhi, qTrue);
          }

          // Who killed them? These four flags are DIAGNOSES, not a decomposition: partnerFake is
          // a subset of crossFake, and one pair can end up in several categories at once.
          bool byCross = false, byPartner = false, byOther = false, byTrue = false;
          auto attribute = [&](int64_t victim, int64_t partner) {
            if (fate.at(victim).alive) {
              return;
            }
            const auto itB = blockersOfPhoton.find(victim);
            if (itB == blockersOfPhoton.end()) {
              return;
            }
            for (const int& w : itB->second) {
              if (cls[w] == kV0CrossLegFake) {
                byCross = true;
                // partner fake = a cross-leg fake using a leg of the OTHER photon of this very
                // pair -> the double-kill mechanism
                if (motherPos[w] == partner || motherEle[w] == partner) {
                  byPartner = true;
                }
              } else if (cls[w] == kV0OtherFake) {
                byOther = true;
              } else if (cls[w] == kV0True) {
                byTrue = true; // real photon against real photon: not repairable by any algorithm
              }
            }
          };
          attribute(mids[a], mids[b]);
          attribute(mids[b], mids[a]);
          if (byCross) {
            fillTruthPair(kLostByCrossFake, qTrue, dEta);
          }
          if (byPartner) {
            fillTruthPair(kLostByPartnerFake, qTrue, dEta);
            fillTruthPairMap(kLostByPartnerFake, dEta, dPhi, qTrue);
          }
          if (byOther) {
            fillTruthPair(kLostByOtherFake, qTrue, dEta);
          }
          if (byTrue) {
            fillTruthPair(kLostByTrue, qTrue, dEta);
          }
        }
      }
    }
  }

  void processRec(MyCollisions const& collisions, FilteredV0s const& v0s, MyTracksIU const& tracks, aod::BCsWithTimestamps const& bcs)
  {
    build<false, false, false>(collisions, v0s, tracks, bcs);
  }
  PROCESS_SWITCH(PhotonConversionBuilder, processRec, "process reconstructed info for data", true);

  // void processRec_SWT(MyCollisionsWithSWT const& collisions, FilteredV0s const& v0s, MyTracksIU const& tracks, aod::BCsWithTimestamps const& bcs)
  // {
  //   build<false, true, false>(collisions, v0s, tracks, bcs);
  // }
  // PROCESS_SWITCH(PhotonConversionBuilder, processRec_SWT, "process reconstructed info for data", false);

  void processMC(MyCollisionsMC const& collisions, FilteredV0s const& v0s, MyTracksIUMC const& tracks,
                 aod::BCsWithTimestamps const& bcs, aod::McParticles const& mcparticles)
  {
    DedupDiag diag;
    build<true, false, false>(collisions, v0s, tracks, bcs, &diag);
    fillDedupTruthDiagnostics(tracks, mcparticles, diag);
  }
  PROCESS_SWITCH(PhotonConversionBuilder, processMC, "process reconstructed info for MC", false);

  void processRec_OnlyIfDielectron(soa::Join<MyCollisions, aod::EMEventsNee> const& collisions, FilteredV0s const& v0s, MyTracksIU const& tracks, aod::BCsWithTimestamps const& bcs)
  {
    build<false, false, true>(collisions, v0s, tracks, bcs);
  }
  PROCESS_SWITCH(PhotonConversionBuilder, processRec_OnlyIfDielectron, "process reconstructed info for data", false);

  // void processRec_SWT_OnlyIfDielectron(soa::Join<MyCollisionsWithSWT, aod::EMEventsNee> const& collisions, FilteredV0s const& v0s, MyTracksIU const& tracks, aod::BCsWithTimestamps const& bcs)
  // {
  //   build<false, true, true>(collisions, v0s, tracks, bcs);
  // }
  // PROCESS_SWITCH(PhotonConversionBuilder, processRec_SWT_OnlyIfDielectron, "process reconstructed info for data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& context)
{
  return WorkflowSpec{
    adaptAnalysisTask<PhotonConversionBuilder>(context, TaskName{"photon-conversion-builder"})};
}
