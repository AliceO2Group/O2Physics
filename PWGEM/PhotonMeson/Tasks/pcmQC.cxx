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

/// \file pcmQC.cxx
/// \brief Quality, purity, material and reconstruction monitoring of the selected V0 photon candidates
///
///   - processQC:      kinematics, conversion-point and leg-track QA of the
///                     reconstructed photons. Material-budget configurables
///                     add conversion-point maps, restricted to selectable
///                     detector regions for material studies.
///   - processQCML:    the same QC with the ML-based photon selection (BDT
///                     models from CCDB) instead of the classical V0 cuts
///   - processPCMQCMC: MC-truth of the QC observables for the purity
///                     assessment, plus resolution assessment.
///                     Configurables activate collision-association studies
///                     and the same regional material-budget mode as in
///                     processQC, there with per-region resolution on top.
///   - processPCMQCMCML: the same MC QC with the ML-based photon selection
///   - processGen:     generator-level distributions as the denominator
///                     reference (requires the binned generated-pT derived
///                     data as input).
///   - processRecoQA:  runs on AO2Ds, before any skimming: track-level QA
///                     plus a truth-based check of the SVertexer candidate
///                     building and deduplication, separating genuine
///                     reconstruction quality from reconstruction artifacts.
///
/// \author Daiki Sekihata <daiki.sekihata@cern.ch> and author Stefanie Mrozinski <stefanie.mrozinski@cern.ch>

#include "PWGEM/Dilepton/Utils/MCUtilities.h"
#include "PWGEM/PhotonMeson/Core/EMPhotonEventCut.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCandidate.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/DataModel/EventTables.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/MCUtilities.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/MathConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DataFormatsParameters/GRPObject.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <MathUtils/Utils.h>

#include <TH1.h>
#include <TH2.h>
#include <TPDGCode.h>

#include <array>
#include <cmath>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::photon;
using namespace o2::aod::pwgem::photonmeson::utils::mcutil;
using namespace o2::aod::pwgem::dilepton::utils::mcutil;

using MyCollisions = soa::Join<aod::PMEvents, aod::EMEventsAlias, aod::EMEventsMult_000, aod::EMEventsCent_000>;
using MyCollision = MyCollisions::iterator;

using MyV0Photons = soa::Join<aod::V0PhotonsKF, aod::V0KFEMEventIds>;
using MyV0Photon = MyV0Photons::iterator;

using MyV0PhotonsML = soa::Join<MyV0Photons, aod::V0PhotonsPhiVPsi>;
using MyV0PhotonML = MyV0PhotonsML::iterator;

// MC Joins
using MyCollisionsMC = soa::Join<aod::PMEvents, aod::EMEventsAlias, aod::EMEventsMult_000, aod::EMEventsCent_000, aod::EMMCEventLabels>;
using MyCollisionMC = MyCollisionsMC::iterator;

using MyMCCollisions = soa::Join<aod::EMMCEvents, aod::BinnedGenPts>;
using MyMCCollision = MyMCCollisions::iterator;

using MyMCV0Legs = soa::Join<aod::V0Legs, aod::V0LegMCLabels>;
using MyMCV0Leg = MyMCV0Legs::iterator;

// AO2D level (processRecoQA): tracks before any skimming, with their MC labels
using MyTracksMC = soa::Join<aod::TracksIU, aod::TracksExtra, aod::McTrackLabels>;

struct PCMQC {
  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};
  Configurable<float> d_bz_input{"d_bz_input", -999, "bz field in kG, -999 is automatic"};

  EMPhotonEventCut fEMEventCut;
  struct : ConfigurableGroup {
    std::string prefix = "eventcut_group";
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", 10.f, "max. Zvtx"};
    Configurable<bool> cfgRequireSel8{"cfgRequireSel8", true, "require sel8 in event cut"};
    Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "require FT0AND in event cut"};
    Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", true, "require No time frame border in event cut"};
    Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", true, "require no ITS readout frame border in event cut"};
    Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "require no same bunch pileup in event cut"};
    Configurable<bool> cfgRequireVertexITSTPC{"cfgRequireVertexITSTPC", false, "require Vertex ITSTPC in event cut"}; // ITS-TPC matched track contributes PV.
    Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", false, "require good Zvtx between FT0 vs. PV in event cut"};
    Configurable<int> cfgTrackOccupancyMin{"cfgTrackOccupancyMin", -2, "min. occupancy"};
    Configurable<int> cfgTrackOccupancyMax{"cfgTrackOccupancyMax", 1000000000, "max. occupancy"};
    Configurable<float> cfgFT0COccupancyMin{"cfgFT0COccupancyMin", -2, "min. FT0C occupancy"};
    Configurable<float> cfgFT0COccupancyMax{"cfgFT0COccupancyMax", 1000000000, "max. FT0C occupancy"};
    Configurable<bool> cfgRequireNoCollInTimeRangeStandard{"cfgRequireNoCollInTimeRangeStandard", false, "require no collision in time range standard"};
    Configurable<bool> cfgRequireNoCollInTimeRangeStrict{"cfgRequireNoCollInTimeRangeStrict", false, "require no collision in time range strict"};
    Configurable<bool> cfgRequireNoCollInITSROFStandard{"cfgRequireNoCollInITSROFStandard", false, "require no collision in time range standard"};
    Configurable<bool> cfgRequireNoCollInITSROFStrict{"cfgRequireNoCollInITSROFStrict", false, "require no collision in time range strict"};
    Configurable<bool> cfgRequireNoHighMultCollInPrevRof{"cfgRequireNoHighMultCollInPrevRof", false, "require no HM collision in previous ITS ROF"};
  } eventcuts;

  V0PhotonCut fV0PhotonCut;
  // v0 photon cuts
  struct : ConfigurableGroup {
    std::string prefix = "v0cut_group";
    Configurable<bool> cfgRequireV0WithITSTPC{"cfgRequireV0WithITSTPC", false, "flag to select V0s with ITS-TPC matched tracks"};
    Configurable<bool> cfgRequireV0WithITSonly{"cfgRequireV0WithITSonly", false, "flag to select V0s with ITSonly tracks"};
    Configurable<bool> cfgRequireV0WithTPConly{"cfgRequireV0WithTPConly", false, "flag to select V0s with TPConly tracks"};
    Configurable<float> cfgMinPtV0{"cfgMinPtV0", 0.1, "min pT for v0 photons at PV"};
    Configurable<float> cfgMaxPtV0{"cfgMaxPtV0", 1e+10, "max pT for v0 photons at PV"};
    Configurable<float> cfgMinEtaV0{"cfgMinEtaV0", -0.8, "min eta for v0 photons at PV"};
    Configurable<float> cfgMaxEtaV0{"cfgMaxEtaV0", +0.8, "max eta for v0 photons at PV"};
    Configurable<float> cfgMinV0Radius{"cfgMinV0Radius", 4.0, "min v0 radius"};
    Configurable<float> cfgMaxV0Radius{"cfgMaxV0Radius", 90.0, "max v0 radius"};
    Configurable<float> cfgMidLV0Radius{"cfgMidLV0Radius", -1.0, "middle low v0 radius for rejection if >0"};
    Configurable<float> cfgMidHV0Radius{"cfgMidHV0Radius", -1.0, "middle high v0 radius for rejection if >0"};
    Configurable<float> cfgMaxAlphaAP{"cfgMaxAlphaAP", 0.95, "max alpha for AP cut"};
    Configurable<float> cfgMaxQtAP{"cfgMaxQtAP", 0.01, "max qT for AP cut"};
    Configurable<float> cfgMinV0CosPA{"cfgMinV0CosPA", 0.999, "min V0 CosPA"};
    Configurable<float> cfgMaxPCA{"cfgMaxPCA", 1.5, "max distance btween 2 legs"};
    Configurable<float> cfgMaxChi2KF{"cfgMaxChi2KF", 1e+10, "max chi2/ndf with KF"};
    Configurable<bool> cfgRejectV0OnITSib{"cfgRejectV0OnITSib", true, "flag to reject V0s on ITSib"};
  } v0cuts;

  // single track cuts
  struct : ConfigurableGroup {
    std::string prefix = "trackcut_group";
    Configurable<int> cfgMinNClustersTPC{"cfgMinNClustersTPC", 0, "min ncluster tpc"};
    Configurable<int> cfgMinNCrossedRows{"cfgMinNCrossedRows", 40, "min crossed rows"};
    Configurable<float> cfgMaxFracSharedClustersTPC{"cfgMaxFracSharedClustersTPC", 999.f, "max fraction of shared clusters in TPC"};
    Configurable<float> cfgMaxChi2TPC{"cfgMaxChi2TPC", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfgMaxChi2ITS{"cfgMaxChi2ITS", 36.0, "max chi2/NclsITS"};
    Configurable<float> cfgMinTPCNsigmaEl{"cfgMinTPCNsigmaEl", -3.0, "min. TPC n sigma for electron"};
    Configurable<float> cfgMaxTPCNsigmaEl{"cfgMaxTPCNsigmaEl", +3.0, "max. TPC n sigma for electron"};
    Configurable<bool> cfgDisableITSonlyTracks{"cfgDisableITSonlyTracks", false, "disable ITSonly tracks in V0 legs"};
    Configurable<bool> cfgDisableTPConlyTracks{"cfgDisableTPConlyTracks", false, "disable TPConly tracks in V0 legs"};
    Configurable<bool> cfgDoDEdxPostCalibration{"cfgDoDEdxPostCalibration", false, "flag to enable dEdx post calibration"};
  } trackcuts;

  // pT-dependent loss QA
  struct : ConfigurableGroup {
    std::string prefix = "qaSettings_group";
    Configurable<bool> cfgDoPtDependentLossQA{"cfgDoPtDependentLossQA", false, "fill the cut variables vs. pT before AND after the V0 selection - the after/before ratio shows which candidates the cuts remove"};
  } qaSettingsGroup;

  // PCM ML inference
  struct : ConfigurableGroup {
    std::string prefix = "mlcut_group";
    Configurable<bool> cfgApplyPCMMl{"cfgApplyPCMMl", false, "Flag to apply ML selections"};
    Configurable<bool> cfgUse2DBinning{"cfgUse2DBinning", false, "Flag to enable/disable 2D binning for ML application"};
    Configurable<bool> cfgLoadModelsFromCCDB{"cfgLoadModelsFromCCDB", true, "Flag to enable or disable the loading of models from CCDB"};
    Configurable<int> cfgTimestampCCDB{"cfgTimestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
    Configurable<int> cfgNClassesPCMMl{"cfgNClassesPCMMl", static_cast<int>(o2::analysis::em_cuts_ml::NCutScores), "Number of classes in ML model"};
    Configurable<std::vector<int>> cfgCutDirPCMMl{"cfgCutDirPCMMl", std::vector<int>{o2::analysis::em_cuts_ml::vecCutDir}, "Whether to reject score values greater or smaller than the threshold"};
    Configurable<std::vector<std::string>> cfgNamesInputFeatures{"cfgNamesInputFeatures", std::vector<std::string>{"feature1", "feature2"}, "Names of ML model input features"};
    Configurable<std::vector<std::string>> cfgModelPathsCCDB{"cfgModelPathsCCDB", std::vector<std::string>{"path_ccdb/BDT_PCM/"}, "Paths of models on CCDB"};
    Configurable<std::vector<std::string>> cfgOnnxFileNames{"cfgOnnxFileNames", std::vector<std::string>{"ModelHandler_onnx_PCM.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
    Configurable<std::vector<std::string>> cfgLabelsBinsPCMMl{"cfgLabelsBinsPCMMl", std::vector<std::string>{"bin 0", "bin 1"}, "Labels for bins"};
    Configurable<std::vector<std::string>> cfgLabelsCutScoresPCMMl{"cfgLabelsCutScoresPCMMl", std::vector<std::string>{o2::analysis::em_cuts_ml::labelsCutScore}, "Labels for cut scores"};
    Configurable<std::vector<double>> cfgBinsPtPCMMl{"cfgBinsPtPCMMl", std::vector<double>{0.0, +1e+10}, "pT bin limits for ML application"};
    Configurable<std::vector<double>> cfgBinsCentPCMMl{"cfgBinsCentPCMMl", std::vector<double>{o2::analysis::em_cuts_ml::vecBinsCent}, "Centrality bin limits for ML application"};
    Configurable<std::vector<double>> cfgCutsPCMMlFlat{"cfgCutsPCMMlFlat", {0.5}, "Flattened ML cuts: [bin0_score0, bin0_score1, ..., binN_scoreM]"};
  } mlcuts;

  struct : ConfigurableGroup {
    std::string prefix = "mcAnalysisModeSettings_group";
    Configurable<bool> cfgDoDetailedResolution{"cfgDoDetailedResolution", false, "purpose of the configurable is to choose if we have more THnSparses included in order to see if some regions have better or worse resolution"};
    Configurable<bool> cfgDoCollisionAssociationQA{"cfgDoCollisionAssociationQA", false, "include the QA Plots which give an idea how photon distribution changes with different subclasses or what the effect of wrong collision association is"};
    Configurable<bool> cfgRequireTrueAssociation{"cfgRequireTrueAssociation", false, "flag to require true mc collision association for the truth-classified QC observables"};
  } mcAnalysisModeSettings;

  struct : ConfigurableGroup {
    std::string prefix = "materialBudgetSettings_group";
    Configurable<bool> cfgDoMaterialDistribution{"cfgDoMaterialDistribution", false, "add THnSparses of the conversion-point distribution for the material budget study"};
    Configurable<bool> cfgDoWiresDetail{"cfgDoWiresDetail", false, "add detailed conversion-point histograms around the wire region"};
    Configurable<bool> cfgDoMFTDetail{"cfgDoMFTDetail", false, "add detailed conversion-point histograms around the MFT region"};
    Configurable<bool> cfgDoITSDetail{"cfgDoITSDetail", false, "add detailed conversion-point histograms around the ITS layers"};
    Configurable<bool> cfgDoTPCInnerBarrelDetail{"cfgDoTPCInnerBarrelDetail", false, "add detailed conversion-point histograms around the TPC inner barrel"};
  } materialBudgetSettingsGroup;

  struct : ConfigurableGroup {
    std::string prefix = "genSettings_group";
    Configurable<float> cfgMaxRGen{"cfgMaxRGen", 90.f, "maximum conversion radius for generated photons (fiducial cut, matches the reco acceptance)"};
    Configurable<float> cfgMarginZMC{"cfgMarginZMC", 7.0, "margin for z cut in cm for MC"};
  } genSettingsGroup;

  struct : ConfigurableGroup {
    std::string prefix = "recoQASettings_group";
    Configurable<bool> cfgDoDetailedTrackQA{"cfgDoDetailedTrackQA", false, "add track-quality histograms of surviving vs. lost conversion legs"};
    Configurable<float> cfgRMinGen{"cfgRMinGen", 1.f, "min true conversion radius (cm)"};
    Configurable<float> cfgRMaxGen{"cfgRMaxGen", 90.f, "max true conversion radius (cm)"};
    Configurable<float> cfgLegEtaMax{"cfgLegEtaMax", 0.9f, "max |eta| of true MC legs"};
    Configurable<float> cfgPhotonEtaMax{"cfgPhotonEtaMax", 0.9f, "max |eta| of true MC photon"};
    Configurable<float> cfgPhotonPtMin{"cfgPhotonPtMin", 0.1f, "min pT of true MC photon"};
  } recoQASettingsGroup;

  o2::ccdb::CcdbApi ccdbApi;
  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb{};
  int mRunNumber = 0;
  float d_bz = 0;
  static constexpr std::array<std::string_view, 2> event_types = {"before/", "after/"};
  static constexpr std::array<std::string_view, 5> mcphoton_types = {"primary/", "fromWD/", "fromHS/", "fromPi0Dalitz/", "fromEtaDalitz/"};
  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  void init(InitContext& context)
  {

    (void)context;
    addhistograms();
    DefineEMEventCut();
    DefinePCMCut();

    mRunNumber = 0;
    d_bz = 0;

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
  }

  template <typename TCollision>
  void initCCDB(TCollision const& collision)
  {
    if (mRunNumber == collision.runNumber()) {
      return;
    }

    // In case override, don't proceed, please - no CCDB access required
    if (d_bz_input > -990) { // o2-linter: disable=magic-number (dummy value to indicate override)
      d_bz = d_bz_input;
      o2::parameters::GRPMagField grpmag;
      if (std::fabs(d_bz) > 1e-5) {                   // o2-linter: disable=magic-number (dummy value to indicate override)
        grpmag.setL3Current(30000.f / (d_bz / 5.0f)); // o2-linter: disable=magic-number (dummy value to indicate override)
      }
      mRunNumber = collision.runNumber();
      return;
    }

    auto run3grp_timestamp = collision.timestamp();
    o2::parameters::GRPObject* grpo = nullptr;
    o2::parameters::GRPMagField* grpmag = nullptr;
    if (!skipGRPOquery) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    }
    if (grpo) {
      // Fetch magnetic field from ccdb for current collision
      d_bz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      // Fetch magnetic field from ccdb for current collision
      d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    }
    fV0PhotonCut.SetD_Bz(d_bz);
    mRunNumber = collision.runNumber();
  }

  void addhistograms()
  {
    // event info
    auto hCollisionCounter = fRegistry.add<TH1>("Event/before/hCollisionCounter", "collision counter;;Number of events", kTH1F, {{10, 0.5, 10.5}}, false);
    hCollisionCounter->GetXaxis()->SetBinLabel(1, "all");
    hCollisionCounter->GetXaxis()->SetBinLabel(2, "No TF border");
    hCollisionCounter->GetXaxis()->SetBinLabel(3, "No ITS ROF border");
    hCollisionCounter->GetXaxis()->SetBinLabel(4, "No Same Bunch Pileup");
    hCollisionCounter->GetXaxis()->SetBinLabel(5, "Is Vertex ITSTPC");
    hCollisionCounter->GetXaxis()->SetBinLabel(6, "Is Good Zvtx FT0vsPV");
    hCollisionCounter->GetXaxis()->SetBinLabel(7, "FT0AND");
    hCollisionCounter->GetXaxis()->SetBinLabel(8, "sel8");
    hCollisionCounter->GetXaxis()->SetBinLabel(9, "|Z_{vtx}| < 10 cm");
    hCollisionCounter->GetXaxis()->SetBinLabel(10, "accepted");

    fRegistry.add("Event/before/hZvtx", "vertex z; Z_{vtx} (cm)", kTH1F, {{100, -50, +50}}, false);
    fRegistry.add("Event/before/hMultNTracksPV", "hMultNTracksPV; N_{track} to PV", kTH1F, {{6001, -0.5, 6000.5}}, false);
    fRegistry.add("Event/before/hMultNTracksPVeta1", "hMultNTracksPVeta1; N_{track} to PV", kTH1F, {{6001, -0.5, 6000.5}}, false);
    fRegistry.add("Event/before/hMultFT0", "hMultFT0;mult. FT0A;mult. FT0C", kTH2F, {{300, 0, 6000}, {300, 0, 6000}}, false);
    fRegistry.add("Event/before/hCentFT0A", "hCentFT0A;centrality FT0A (%)", kTH1F, {{110, 0, 110}}, false);
    fRegistry.add("Event/before/hCentFT0C", "hCentFT0C;centrality FT0C (%)", kTH1F, {{110, 0, 110}}, false);
    fRegistry.add("Event/before/hCentFT0M", "hCentFT0M;centrality FT0M (%)", kTH1F, {{110, 0, 110}}, false);
    fRegistry.add("Event/before/hCentFT0MvsMultNTracksPV", "hCentFT0MvsMultNTracksPV;centrality FT0M (%);N_{track} to PV", kTH2F, {{110, 0, 110}, {600, 0, 6000}}, false);
    fRegistry.add("Event/before/hMultFT0MvsMultNTracksPV", "hMultFT0MvsMultNTracksPV;mult. FT0M;N_{track} to PV", kTH2F, {{600, 0, 6000}, {600, 0, 6000}}, false);
    fRegistry.addClone("Event/before/", "Event/after/");

    // v0 info
    fRegistry.add("V0/hPt", "pT;p_{T,#gamma} (GeV/c)", kTH1F, {{2000, 0.0f, 20}}, false);
    fRegistry.add("V0/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{90, 0, o2::constants::math::TwoPI}, {200, -1.0f, 1.0f}}, false);
    fRegistry.add("V0/hXY", "conversion point in XY;V_{x} (cm);V_{y} (cm)", kTH2F, {{400, -100.0f, 100.0f}, {400, -100.0f, 100.0f}}, false);
    fRegistry.add("V0/hRZ", "conversion point in RZ;Z (cm);R_{xy} (cm)", kTH2F, {{200, -100, 100}, {200, 0.0f, 100.0f}}, false);
    fRegistry.add("V0/hCosPA", "V0CosPA;cosine pointing angle in 3D", kTH1F, {{100, 0.99f, 1.0f}}, false);
    fRegistry.add("V0/hCosPAXY", "V0CosPA;cosine pointing angle in XY", kTH1F, {{100, 0.99f, 1.0f}}, false);
    fRegistry.add("V0/hCosPARZ", "V0CosPA;cosine pointing angle in RZ", kTH1F, {{100, 0.99f, 1.0f}}, false);
    fRegistry.add("V0/hPCA", "distance between 2 legs;PCA (cm)", kTH1F, {{500, 0.0f, 5.0f}}, false);
    fRegistry.add("V0/hDCAxyz", "DCA to PV;DCA_{xy} (cm);DCA_{z} (cm)", kTH2F, {{200, -5.f, +5.f}, {200, -5.f, +5.f}}, false);
    fRegistry.add("V0/hDCAz_Pt", "DCA_{z} to PV vs. p_{T};DCA_{z} (cm);p_{T} (GeV/c)", kTH2F, {{200, -5.f, +5.f}, {2000, 0.0f, 20}}, false);
    fRegistry.add("V0/hAPplot", "AP plot;#alpha;q_{T} (GeV/c)", kTH2F, {{200, -1.0f, +1.0f}, {250, 0.0f, 0.25f}}, false);
    fRegistry.add("V0/hRxyVsPt", "conversion radius vs. pT;p_{T,#gamma} (GeV/c);R_{xy} (cm)", kTH2F, {{100, 0, 10}, {200, 0, 100}}, false);
    fRegistry.add("V0/hEtaVsPt", "#eta vs. pT;p_{T,#gamma} (GeV/c);#eta", kTH2F, {{100, 0, 10}, {200, -1.0f, 1.0f}}, false);
    fRegistry.add("V0/hPhiVsPt", "#varphi vs. pT;p_{T,#gamma} (GeV/c);#varphi (rad.)", kTH2F, {{100, 0, 10}, {90, 0, o2::constants::math::TwoPI}}, false);
    fRegistry.add("V0/hsAlphaQtPt", "Armenteros vs. pT;#alpha;q_{T} (GeV/c);p_{T,#gamma} (GeV/c)", kTHnSparseF, {{100, -1.0f, +1.0f}, {125, 0.0f, 0.25f}, {100, 0, 10}}, false);
    fRegistry.add("V0/hPsiPairVsPt", "#psi_{pair} vs. pT;p_{T,#gamma} (GeV/c);#psi_{pair} (rad.)", kTH2F, {{100, 0, 10}, {200, -0.5f, +0.5f}}, false);
    fRegistry.add("V0/hMassGamma", "hMassGamma;R_{xy} (cm);m_{ee} (GeV/c^{2})", kTH2F, {{200, 0.0f, 100.0f}, {100, 0.0f, 0.1f}}, false);
    fRegistry.add("V0/hKFChi2vsM", "KF chi2 vs. m_{ee};m_{ee} (GeV/c^{2});KF chi2/NDF", kTH2F, {{100, 0.0f, 0.1f}, {100, 0.f, 100.0f}}, false);
    fRegistry.add("V0/hKFChi2vsR", "KF chi2 vs. conversion point in XY;R_{xy} (cm);KF chi2/NDF", kTH2F, {{200, 0.0f, 100.0f}, {100, 0.f, 100.0f}}, false);
    fRegistry.add("V0/hKFChi2vsX", "KF chi2 vs. conversion point in X;X (cm);KF chi2/NDF", kTH2F, {{200, -100.0f, 100.0f}, {100, 0.f, 100.0f}}, false);
    fRegistry.add("V0/hKFChi2vsY", "KF chi2 vs. conversion point in Y;Y (cm);KF chi2/NDF", kTH2F, {{200, -100.0f, 100.0f}, {100, 0.f, 100.0f}}, false);
    fRegistry.add("V0/hKFChi2vsZ", "KF chi2 vs. conversion point in Z;Z (cm);KF chi2/NDF", kTH2F, {{200, -100.0f, 100.0f}, {100, 0.f, 100.0f}}, false);
    fRegistry.add("V0/hsConvPoint", "photon conversion point;r_{xy} (cm);#varphi (rad.);#eta;", kTHnSparseF, {{100, 0.0f, 100}, {90, 0, o2::constants::math::TwoPI}, {80, -2, +2}}, false);
    fRegistry.add("V0/hNgamma", "Number of #gamma candidates per collision", kTH1F, {{101, -0.5f, 100.5f}});

    if (mlcuts.cfgApplyPCMMl) {
      if (mlcuts.cfgNClassesPCMMl == 2) { // o2-linter: disable=magic-number (BDT class)
        fRegistry.add("V0/hBDTBackgroundScoreVsPt", "BDT background score vs pT; pT (GeV/c); BDT background score", {HistType::kTH2F, {{1000, 0.0f, 20.0f}, {1000, 0.0f, 1.0f}}});
        fRegistry.add("V0/hBDTSignalScoreVsPt", "BDT signal score vs pT; pT (GeV/c); BDT signal score", {HistType::kTH2F, {{1000, 0.0f, 20.0f}, {1000, 0.0f, 1.0f}}});
        fRegistry.add("V0/hPhiVPsi", "#varphi vs. #psi angle;#psi (rad.); #varphi (rad.)", kTH2F, {{200, -o2::constants::math::PI, o2::constants::math::PI}, {200, 0, o2::constants::math::TwoPI}}, false);
      } else if (mlcuts.cfgNClassesPCMMl == 3) { // o2-linter: disable=magic-number (BDT class)
        fRegistry.add("V0/hBDTBackgroundScoreVsPt", "BDT background score vs pT; pT (GeV/c); BDT background score", {HistType::kTH2F, {{1000, 0.0f, 20.0f}, {1000, 0.0f, 1.0f}}});
        fRegistry.add("V0/hBDTPrimaryPhotonScoreVsPt", "BDT primary photon score vs pT; pT (GeV/c); BDT primary photon score", {HistType::kTH2F, {{1000, 0.0f, 20.0f}, {1000, 0.0f, 1.0f}}});
        fRegistry.add("V0/hBDTSecondaryPhotonScoreVsPt", "BDT secondary photon score vs pT; pT (GeV/c); BDT secondary photon score", {HistType::kTH2F, {{1000, 0.0f, 20.0f}, {1000, 0.0f, 1.0f}}});
        fRegistry.add("V0/hPhiVPsi", "#varphi vs. #psi angle;#psi (rad.); #varphi (rad.)", kTH2F, {{200, -o2::constants::math::PI, o2::constants::math::PI}, {200, 0, o2::constants::math::TwoPI}}, false);
      } else {
        fRegistry.add("V0/hBDTScoreVsPt", "BDT score vs pT; pT (GeV/c); BDT score", {HistType::kTH2F, {{1000, 0.0f, 20.0f}, {1000, 0.0f, 1.0f}}});
        fRegistry.add("V0/hPhiVPsi", "#varphi vs. #psi angle;#psi (rad.); #varphi (rad.)", kTH2F, {{200, -o2::constants::math::PI, o2::constants::math::PI}, {200, 0, o2::constants::math::TwoPI}}, false);
      }
    }

    // v0leg info
    fRegistry.add("V0Leg/hPt", "pT;p_{T,e} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);
    fRegistry.add("V0Leg/hQoverPt", "q/pT;q/p_{T} (GeV/c)^{-1}", kTH1F, {{1000, -50, 50}}, false);
    fRegistry.add("V0Leg/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{90, 0, o2::constants::math::TwoPI}, {200, -1.0f, 1.0f}}, false);
    fRegistry.add("V0Leg/hDCAxyz", "DCA xy vs. z;DCA_{xy} (cm);DCA_{z} (cm)", kTH2F, {{200, -50.0f, 50.0f}, {200, -50.0f, 50.0f}}, false);
    fRegistry.add("V0Leg/hNclsTPC", "number of TPC clusters", kTH1F, {{161, -0.5, 160.5}}, false);
    fRegistry.add("V0Leg/hNcrTPC", "number of TPC crossed rows", kTH1F, {{161, -0.5, 160.5}}, false);
    fRegistry.add("V0Leg/hChi2TPC", "chi2/number of TPC clusters", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("V0Leg/hTPCdEdx", "TPC dE/dx;p_{in} (GeV/c);TPC dE/dx (a.u.)", kTH2F, {{1000, 0, 10}, {200, 0, 200}}, false);
    fRegistry.add("V0Leg/hTPCNsigmaEl", "TPC n sigma el;p_{in} (GeV/c);n #sigma_{e}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("V0Leg/hTPCNsigmaPi", "TPC n sigma pi;p_{in} (GeV/c);n #sigma_{#pi}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("V0Leg/hTPCNsigmaElVsEta", "TPC n sigma el vs. eta;#eta;n #sigma_{e}^{TPC}", kTH2F, {{40, -1.0f, 1.0f}, {100, -5, +5}}, false); // eta-dependence = dE/dx calibration check (TOF/ITS n#sigma not stored in the V0Leg skim)
    fRegistry.add("V0Leg/hTPCNcr2Nf", "TPC Ncr/Nfindable", kTH1F, {{200, 0, 2}}, false);
    fRegistry.add("V0Leg/hTPCNcls2Nf", "TPC Ncls/Nfindable", kTH1F, {{200, 0, 2}}, false);
    fRegistry.add("V0Leg/hTPCNclsShared", "TPC Ncls shared/Ncls;p_{T} (GeV/c);N_{cls}^{shared}/N_{cls} in TPC", kTH2F, {{1000, 0, 10}, {100, 0, 1}}, false);
    fRegistry.add("V0Leg/hNclsITS", "number of ITS clusters", kTH1F, {{8, -0.5, 7.5}}, false);
    fRegistry.add("V0Leg/hChi2ITS", "chi2/number of ITS clusters", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("V0Leg/hITSClusterMap", "ITS cluster map", kTH1F, {{128, -0.5, 127.5}}, false);
    fRegistry.add("V0Leg/hMeanClusterSizeITS", "mean cluster size ITS;<cluster size> on ITS #times cos(#lambda)", kTH2F, {{1000, 0, 10}, {160, 0, 16}}, false);
    if (trackcuts.cfgDoDEdxPostCalibration) {
      fRegistry.add("V0Leg/hPvsConvPointvsTPCNsigmaElvsEta_Pos", "momentum of pos leg vs. conversion point of V0 vs. TPC n sigma pos vs. eta of pos leg; p (GeV/c); r_{xy} (cm); n #sigma_{e}^{TPC}; #eta", kTHnSparseF, {{200, 0, 20}, {100, 0, 100}, {500, -5, 5}, {200, -1, +1}}, false);
      fRegistry.add("V0Leg/hPvsConvPointvsTPCNsigmaElvsEta_Ele", "momentum of neg leg vs. conversion point of V0 vs. TPC n sigma el vs. eta of neg leg; p (GeV/c); r_{xy} (cm); n #sigma_{e}^{TPC}; #eta", kTHnSparseF, {{200, 0, 20}, {100, 0, 100}, {500, -5, 5}, {200, -1, +1}}, false);
    }

    if (materialBudgetSettingsGroup.cfgDoMaterialDistribution) {
      fRegistry.add("MaterialBudget/hs", "conversion point;z (cm);r_{xy} (cm);#eta;#varphi (rad.);p_{T,#gamma} (GeV/c)", kTHnSparseF, {{200, -100, 100}, {100, 0, 100}, {80, -2, +2}, {90, 0, o2::constants::math::TwoPI}, {100, 0, 10}}, false);
    }
    if (materialBudgetSettingsGroup.cfgDoWiresDetail) {
      fRegistry.add("MaterialBudget/Wires/hsLeft", "conversion point near left wire;x (cm);y (cm);z (cm);#varphi (rad.);r_{xy} (cm);p_{T,#gamma} (GeV/c)", kTHnSparseF, {{200, -20, 20}, {200, -20, 20}, {80, -20, 20}, {40, 3.15, 3.4}, {80, 0, 20}, {100, 0, 10}}, false);
      fRegistry.add("MaterialBudget/Wires/hsRight", "conversion point near right wire;x (cm);y (cm);z (cm);#varphi (rad.);r_{xy} (cm);p_{T,#gamma} (GeV/c)", kTHnSparseF, {{200, -20, 20}, {200, -20, 20}, {80, -20, 20}, {40, 6.00, 6.15}, {80, 0, 20}, {100, 0, 10}}, false);
    }
    if (materialBudgetSettingsGroup.cfgDoITSDetail) {
      fRegistry.add("MaterialBudget/ITS/hs", "conversion point around the ITS layers;x (cm);y (cm);z (cm);#varphi (rad.);r_{xy} (cm);p_{T,#gamma} (GeV/c)", kTHnSparseF, {{160, -40, 40}, {160, -40, 40}, {80, -40, 40}, {90, 0, o2::constants::math::TwoPI}, {150, 0, 60}, {100, 0, 10}}, false);
    }
    if (materialBudgetSettingsGroup.cfgDoMFTDetail) {
      fRegistry.add("MaterialBudget/MFT/hs", "conversion point around the MFT region;x (cm);y (cm);z (cm);#varphi (rad.);r_{xy} (cm);p_{T,#gamma} (GeV/c)", kTHnSparseF, {{160, -80, 80}, {160, -80, 80}, {40, -40, 40}, {90, 0, o2::constants::math::TwoPI}, {40, 40, 60}, {100, 0, 10}}, false);
    }
    if (materialBudgetSettingsGroup.cfgDoTPCInnerBarrelDetail) {
      fRegistry.add("MaterialBudget/TPC/hs", "conversion point around the TPC inner barrel;x (cm);y (cm);z (cm);#varphi (rad.);r_{xy} (cm);p_{T,#gamma} (GeV/c)", kTHnSparseF, {{160, -90, 90}, {160, -90, 90}, {80, -40, 40}, {90, 0, o2::constants::math::TwoPI}, {80, 60, 80}, {100, 0, 10}}, false);
    }

    if (doprocessGen) {
      std::vector<double> ptbins;
      ptbins.reserve(72);
      for (int i = 0; i < 2; i++) {                // o2-linter: disable=magic-number (binning)
        ptbins.emplace_back(0.05 * (i - 0) + 0.0); // from 0 to 0.05 GeV/c, every 0.05 GeV/c
      }
      for (int i = 2; i < 51; i++) {              // o2-linter: disable=magic-number (binning)
        ptbins.emplace_back(0.1 * (i - 2) + 0.1); // from 0.1 to 4.9 GeV/c, every 0.1 GeV/c
      }
      for (int i = 51; i < 61; i++) {              // o2-linter: disable=magic-number (binning)
        ptbins.emplace_back(0.5 * (i - 51) + 5.0); // from 5 to 9.5 GeV/c, every 0.5 GeV/c
      }
      for (int i = 61; i < 72; i++) {               // o2-linter: disable=magic-number (binning)
        ptbins.emplace_back(1.0 * (i - 61) + 10.0); // from 10 to 20 GeV/c, every 1 GeV/c
      }
      const AxisSpec axisPtGen{ptbins, "p_{T,#gamma} (GeV/c)"};
      const AxisSpec axisRapidityGen{{0.0, +0.8, +0.9}, "rapidity |y_{#gamma}|"};

      fRegistry.add("Generated/hPt", "pT;p_{T} (GeV/c)", kTH1D, {axisPtGen}, true);
      fRegistry.add("Generated/hPtY", "Generated info", kTH2D, {axisPtGen, axisRapidityGen}, true);
      fRegistry.add("Generated/hPt_ConversionPhoton", "converted photon pT;p_{T} (GeV/c)", kTH1D, {axisPtGen}, true);
      fRegistry.add("Generated/hY_ConversionPhoton", "converted photon y;rapidity y", kTH1F, {{40, -2.0f, 2.0f}}, true);
      fRegistry.add("Generated/hPhi_ConversionPhoton", "converted photon #varphi;#varphi (rad.)", kTH1F, {{180, 0, o2::constants::math::TwoPI}}, true);
      fRegistry.add("Generated/hXY", "conversion point in XY MC;V_{x} (cm);V_{y} (cm)", kTH2F, {{800, -100.0f, 100.0f}, {800, -100.0f, 100.0f}}, true);
      fRegistry.add("Generated/hRZ", "conversion point in RZ MC;V_{z} (cm);R_{xy} (cm)", kTH2F, {{400, -100.0f, 100.0f}, {400, 0.f, 100.0f}}, true);
      fRegistry.add("Generated/hRPhi", "conversion point of #varphi vs. R_{xy} MC;#varphi (rad.);R_{xy} (cm);N_{e}", kTH2F, {{360, 0.0f, o2::constants::math::TwoPI}, {400, 0, 100}}, true);
      fRegistry.add("Generated/hsConvPoint", "photon conversion point;r_{xy} (cm);#varphi (rad.);#eta;", kTHnSparseF, {{100, 0.0f, 100}, {90, 0, o2::constants::math::TwoPI}, {80, -2, +2}}, true);
      fRegistry.add("Generated/hR_ConversionPhoton_wideR", "converted photon R_{xy} up to the EMCal, before fiducial cuts;R_{xy} (cm);N_{#gamma}", kTH1F, {{250, 0.0f, 500.0f}}, true);
      fRegistry.add("Generated/hRZ_wideR", "conversion point in RZ up to the EMCal, before fiducial cuts;V_{z} (cm);R_{xy} (cm)", kTH2F, {{200, -400.0f, 400.0f}, {250, 0.0f, 500.0f}}, true);
      fRegistry.add("Generated/hsRZPhi_wideR", "conversion point up to the EMCal, before fiducial cuts;V_{z} (cm);#varphi (rad.);R_{xy} (cm)", kTHnSparseF, {{200, -400.0f, 400.0f}, {360, 0.0f, o2::constants::math::TwoPI}, {250, 0.0f, 500.0f}}, true);
      if (mcAnalysisModeSettings.cfgDoDetailedResolution) {
        fRegistry.add("Generated/hPtEtaPhi", "Photon pt vs eta, and phi;p_{T, gen} (GeV/c);#eta;#phi", kTHnSparseF, {{200, 0., 20.}, {18, -0.9, 0.9}, {36, 0, o2::constants::math::TwoPI}}, false);
      }
    }

    if (doprocessPCMQCMC || doprocessPCMQCMCML) {
      fRegistry.add("V0/primary/hPt", "pT;p_{T,#gamma} (GeV/c)", kTH1F, {{2000, 0.0f, 20}}, false);
      fRegistry.add("V0/primary/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{90, 0, o2::constants::math::TwoPI}, {200, -1.0f, 1.0f}}, false);
      fRegistry.add("V0/primary/hXY", "conversion point in XY;V_{x} (cm);V_{y} (cm)", kTH2F, {{400, -100.0f, 100.0f}, {400, -100.0f, 100.0f}}, false);
      fRegistry.add("V0/primary/hRZ", "conversion point in RZ;Z (cm);R_{xy} (cm)", kTH2F, {{200, -100, 100}, {200, 0.0f, 100.0f}}, false);
      fRegistry.add("V0/primary/hCosPA", "V0CosPA;cosine pointing angle in 3D", kTH1F, {{100, 0.99f, 1.0f}}, false);
      fRegistry.add("V0/primary/hCosPAXY", "V0CosPA;cosine pointing angle in XY", kTH1F, {{100, 0.99f, 1.0f}}, false);
      fRegistry.add("V0/primary/hCosPARZ", "V0CosPA;cosine pointing angle in RZ", kTH1F, {{100, 0.99f, 1.0f}}, false);
      fRegistry.add("V0/primary/hPCA", "distance between 2 legs;PCA (cm)", kTH1F, {{500, 0.0f, 5.0f}}, false);
      fRegistry.add("V0/primary/hDCAxyz", "DCA to PV;DCA_{xy} (cm);DCA_{z} (cm)", kTH2F, {{200, -5.f, +5.f}, {200, -5.f, +5.f}}, false);
      fRegistry.add("V0/primary/hDCAz_Pt", "DCA_{z} to PV vs. p_{T};DCA_{z} (cm);p_{T} (GeV/c)", kTH2F, {{200, -5.f, +5.f}, {2000, 0.0f, 20}}, false);
      fRegistry.add("V0/primary/hAPplot", "AP plot;#alpha;q_{T} (GeV/c)", kTH2F, {{200, -1.0f, +1.0f}, {250, 0.0f, 0.25f}}, false);
      fRegistry.add("V0/primary/hRxyVsPt", "conversion radius vs. pT;p_{T,#gamma} (GeV/c);R_{xy} (cm)", kTH2F, {{100, 0, 10}, {200, 0, 100}}, false);
      fRegistry.add("V0/primary/hEtaVsPt", "#eta vs. pT;p_{T,#gamma} (GeV/c);#eta", kTH2F, {{100, 0, 10}, {200, -1.0f, 1.0f}}, false);
      fRegistry.add("V0/primary/hPhiVsPt", "#varphi vs. pT;p_{T,#gamma} (GeV/c);#varphi (rad.)", kTH2F, {{100, 0, 10}, {90, 0, o2::constants::math::TwoPI}}, false);
      fRegistry.add("V0/primary/hsAlphaQtPt", "Armenteros vs. pT;#alpha;q_{T} (GeV/c);p_{T,#gamma} (GeV/c)", kTHnSparseF, {{100, -1.0f, +1.0f}, {125, 0.0f, 0.25f}, {100, 0, 10}}, false);
      fRegistry.add("V0/primary/hPsiPairVsPt", "#psi_{pair} vs. pT;p_{T,#gamma} (GeV/c);#psi_{pair} (rad.)", kTH2F, {{100, 0, 10}, {200, -0.5f, +0.5f}}, false);
      fRegistry.add("V0/primary/hMassGamma", "hMassGamma;R_{xy} (cm);m_{ee} (GeV/c^{2})", kTH2F, {{200, 0.0f, 100.0f}, {100, 0.0f, 0.1f}}, false);
      fRegistry.add("V0/primary/hKFChi2vsM", "KF chi2 vs. m_{ee};m_{ee} (GeV/c^{2});KF chi2/NDF", kTH2F, {{100, 0.0f, 0.1f}, {100, 0.f, 100.0f}}, false);
      fRegistry.add("V0/primary/hKFChi2vsR", "KF chi2 vs. conversion point in XY;R_{xy} (cm);KF chi2/NDF", kTH2F, {{200, 0.0f, 100.0f}, {100, 0.f, 100.0f}}, false);
      fRegistry.add("V0/primary/hKFChi2vsX", "KF chi2 vs. conversion point in X;X (cm);KF chi2/NDF", kTH2F, {{200, -100.0f, 100.0f}, {100, 0.f, 100.0f}}, false);
      fRegistry.add("V0/primary/hKFChi2vsY", "KF chi2 vs. conversion point in Y;Y (cm);KF chi2/NDF", kTH2F, {{200, -100.0f, 100.0f}, {100, 0.f, 100.0f}}, false);
      fRegistry.add("V0/primary/hKFChi2vsZ", "KF chi2 vs. conversion point in Z;Z (cm);KF chi2/NDF", kTH2F, {{200, -100.0f, 100.0f}, {100, 0.f, 100.0f}}, false);
      fRegistry.add("V0/primary/hNgamma", "Number of true #gamma per collision;N_{#gamma} per event;Number of events", kTH1F, {{101, -0.5f, 100.5f}});
      fRegistry.add("V0/primary/hConvPoint_diffX", "conversion point diff X MC;X_{MC} (cm);X_{rec} - X_{MC} (cm)", kTH2F, {{200, -100, +100}, {100, -50.0f, 50.0f}}, true);
      fRegistry.add("V0/primary/hConvPoint_diffY", "conversion point diff Y MC;Y_{MC} (cm);Y_{rec} - Y_{MC} (cm)", kTH2F, {{200, -100, +100}, {100, -50.0f, 50.0f}}, true);
      fRegistry.add("V0/primary/hConvPoint_diffZ", "conversion point diff Z MC;Z_{MC} (cm);Z_{rec} - Z_{MC} (cm)", kTH2F, {{200, -100, +100}, {100, -50.0f, 50.0f}}, true);
      fRegistry.add("V0/primary/hPtGen_DeltaPtOverPtGen", "photon p_{T} resolution;p_{T}^{gen} (GeV/c);(p_{T}^{rec} - p_{T}^{gen})/p_{T}^{gen}", kTH2F, {{1000, 0, 10}, {200, -1.0f, 1.0f}}, true);
      fRegistry.add("V0/primary/hPtGen_DeltaEta", "photon #eta resolution;p_{T}^{gen} (GeV/c);#eta^{rec} - #eta^{gen}", kTH2F, {{200, 0, 20}, {199, -0.099, 0.099}}, true);
      fRegistry.add("V0/primary/hPtGen_DeltaPhi", "photon #varphi resolution;p_{T}^{gen} (GeV/c);#varphi^{rec} - #varphi^{gen} (rad.)", kTH2F, {{200, 0, 20}, {199, -0.099, 0.099}}, true);
      fRegistry.add("V0/primary/hRxyGen_DeltaPtOverPtGen", "photon p_{T} resolution; R_{xy}^{gen} (cm);(p_{T}^{rec} - p_{T}^{gen})/p_{T}^{gen}", kTH2F, {{100, 0, 100}, {200, -1.0f, 1.0f}}, true);
      fRegistry.add("V0/primary/hRxyGen_DeltaEta", "photon #eta resolution;R_{xy}^{gen} (cm);#eta^{rec} - #eta^{gen}", kTH2F, {{100, 0, 100}, {100, -0.5f, 0.5f}}, true);
      fRegistry.add("V0/primary/hRxyGen_DeltaPhi", "photon #varphi resolution;R_{xy}^{gen} (cm);#varphi^{rec} - #varphi^{gen} (rad.)", kTH2F, {{100, 0, 100}, {100, -0.5f, 0.5f}}, true);
      fRegistry.add("V0/primary/hRxyGen_DeltaR", "photon R_{xy} resolution;R_{xy}^{gen} (cm);R_{xy}^{rec} - R_{xy}^{gen} (cm)", kTH2F, {{100, 0, 100}, {100, 0, 100}}, true);
      fRegistry.add("V0/primary/hXY_MC", "X vs. Y of true photon conversion point.;X (cm);Y (cm)", kTH2F, {{400, -100.0f, +100}, {400, -100, +100}}, true);
      fRegistry.add("V0/primary/hRZ_MC", "R vs. Z of true photon conversion point;Z (cm);R_{xy} (cm)", kTH2F, {{200, -100.0f, +100}, {200, 0, 100}}, true);
      fRegistry.add("V0/primary/hsConvPoint", "photon conversion point;r_{xy} (cm);#varphi (rad.);#eta;", kTHnSparseF, {{100, 0.0f, 100}, {90, 0, o2::constants::math::TwoPI}, {80, -2, +2}}, false);
      if (mcAnalysisModeSettings.cfgDoDetailedResolution) {
        fRegistry.add("V0/primary/hEtaPhiResol", "Photon eta-phi resolution;p_{T} (GeV/c);#eta-diff;#phi-diff", kTH3F, {{200, 0., 20.}, {99, -0.049, 0.049}, {99, -0.049, 0.049}}, false);
        fRegistry.add("V0/primary/hPtResolPtEtaPhi", "Photon resolution vs. pt, eta, and phi;p_{T_rec} - p_{T true} / p_{T true};p_{T} (GeV/c);#eta;#phi", kTHnSparseF, {{199, -0.995, 0.995}, {200, 0., 20.}, {18, -0.9, 0.9}, {36, 0, o2::constants::math::TwoPI}}, false);
        fRegistry.add("V0/primary/hMomResolPtEtaPhi", "Photon momentum resolution vs. p, eta, and phi;p_{rec} - p_{true} / p_{true};p_{T} (GeV/c);#eta;#phi", kTHnSparseF, {{199, -0.995, 0.995}, {200, 0., 20.}, {18, -0.9, 0.9}, {36, 0, o2::constants::math::TwoPI}}, false);
        fRegistry.add("V0/primary/hPtEtaPhi", "pt, eta, and phi;p_{T} (GeV/c);#eta;#phi", kTHnSparseF, {{200, 0., 20.}, {18, -0.9, 0.9}, {36, 0, o2::constants::math::TwoPI}}, false);
      }
      if (mlcuts.cfgApplyPCMMl) {
        if (mlcuts.cfgNClassesPCMMl == 2) { // o2-linter: disable=magic-number (BDT group)
          fRegistry.add("V0/primary/hBDTBackgroundScoreVsPt", "BDT background score vs pT; pT (GeV/c); BDT background score", {HistType::kTH2F, {{200, 0.0f, 20.0f}, {100, 0.0f, 1.0f}}});
          fRegistry.add("V0/primary/hBDTSignalScoreVsPt", "BDT signal score vs pT; pT (GeV/c); BDT signal score", {HistType::kTH2F, {{200, 0.0f, 20.0f}, {100, 0.0f, 1.0f}}});
        } else if (mlcuts.cfgNClassesPCMMl == 3) { // o2-linter: disable=magic-number (BDT group)
          fRegistry.add("V0/primary/hBDTBackgroundScoreVsPt", "BDT background score vs pT; pT (GeV/c); BDT background score", {HistType::kTH2F, {{200, 0.0f, 20.0f}, {100, 0.0f, 1.0f}}});
          fRegistry.add("V0/primary/hBDTPrimaryPhotonScoreVsPt", "BDT primary photon score vs pT; pT (GeV/c); BDT primary photon score", {HistType::kTH2F, {{200, 0.0f, 20.0f}, {100, 0.0f, 1.0f}}});
          fRegistry.add("V0/primary/hBDTSecondaryPhotonScoreVsPt", "BDT secondary photon score vs pT; pT (GeV/c); BDT secondary photon score", {HistType::kTH2F, {{200, 0.0f, 20.0f}, {100, 0.0f, 1.0f}}});
        } else {
          fRegistry.add("V0/primary/hBDTScoreVsPt", "BDT score vs pT; pT (GeV/c); BDT score", {HistType::kTH2F, {{200, 0.0f, 20.0f}, {100, 0.0f, 1.0f}}});
        }
        fRegistry.add("V0/primary/hPhiVPsi", "#varphi vs. #psi angle;#psi (rad.); #varphi (rad.)", kTH2F, {{200, -o2::constants::math::PI, o2::constants::math::PI}, {200, 0, o2::constants::math::TwoPI}}, false);
      }
      fRegistry.addClone("V0/primary/", "V0/fromWD/");                  // from weak decay
      fRegistry.addClone("V0/primary/", "V0/fromHS/");                  // from hadronic shower in detector materials
      fRegistry.addClone("V0/primary/", "V0/fromPi0Dalitz/");           // misidentified dielectron from pi0 dalitz decay
      fRegistry.addClone("V0/primary/", "V0/fromEtaDalitz/");           // misidentified dielectron from eta dalitz decay
      fRegistry.addClone("V0/primary/hPt", "V0/candidate/hPt");         // only for purity
      fRegistry.addClone("V0/primary/hEtaPhi", "V0/candidate/hEtaPhi"); // only for purity

      fRegistry.add("V0Leg/primary/hPt", "pT;p_{T,e} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);
      fRegistry.add("V0Leg/primary/hQoverPt", "q/pT;q/p_{T} (GeV/c)^{-1}", kTH1F, {{1000, -50, 50}}, false);
      fRegistry.add("V0Leg/primary/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{90, 0, o2::constants::math::TwoPI}, {200, -1.0f, 1.0f}}, false);
      fRegistry.add("V0Leg/primary/hDCAxyz", "DCA xy vs. z;DCA_{xy} (cm);DCA_{z} (cm)", kTH2F, {{200, -50.0f, 50.0f}, {200, -50.0f, 50.0f}}, false);
      fRegistry.add("V0Leg/primary/hNclsTPC", "number of TPC clusters", kTH1F, {{161, -0.5, 160.5}}, false);
      fRegistry.add("V0Leg/primary/hNcrTPC", "number of TPC crossed rows", kTH1F, {{161, -0.5, 160.5}}, false);
      fRegistry.add("V0Leg/primary/hChi2TPC", "chi2/number of TPC clusters", kTH1F, {{100, 0, 10}}, false);
      fRegistry.add("V0Leg/primary/hTPCdEdx", "TPC dE/dx;p_{in} (GeV/c);TPC dE/dx (a.u.)", kTH2F, {{1000, 0, 10}, {200, 0, 200}}, false);
      fRegistry.add("V0Leg/primary/hTPCNsigmaEl", "TPC n sigma el;p_{in} (GeV/c);n #sigma_{e}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("V0Leg/primary/hTPCNsigmaPi", "TPC n sigma pi;p_{in} (GeV/c);n #sigma_{#pi}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("V0Leg/primary/hTPCNsigmaElVsEta", "TPC n sigma el vs. eta;#eta;n #sigma_{e}^{TPC}", kTH2F, {{40, -1.0f, 1.0f}, {100, -5, +5}}, false);
      fRegistry.add("V0Leg/primary/hTPCNcr2Nf", "TPC Ncr/Nfindable", kTH1F, {{200, 0, 2}}, false);
      fRegistry.add("V0Leg/primary/hTPCNcls2Nf", "TPC Ncls/Nfindable", kTH1F, {{200, 0, 2}}, false);
      fRegistry.add("V0Leg/primary/hTPCNclsShared", "TPC Ncls shared/Ncls;p_{T} (GeV/c);N_{cls}^{shared}/N_{cls} in TPC", kTH2F, {{1000, 0, 10}, {100, 0, 1}}, false);
      fRegistry.add("V0Leg/primary/hNclsITS", "number of ITS clusters", kTH1F, {{8, -0.5, 7.5}}, false);
      fRegistry.add("V0Leg/primary/hChi2ITS", "chi2/number of ITS clusters", kTH1F, {{100, 0, 10}}, false);
      fRegistry.add("V0Leg/primary/hITSClusterMap", "ITS cluster map", kTH1F, {{128, -0.5, 127.5}}, false);
      fRegistry.add("V0Leg/primary/hMeanClusterSizeITS", "mean cluster size ITS;<cluster size> on ITS #times cos(#lambda)", kTH2F, {{1000, 0, 10}, {160, 0, 16}}, false);
      fRegistry.add("V0Leg/primary/hPtGen_DeltaPtOverPtGen", "electron p_{T} resolution;p_{T}^{gen} (GeV/c);(p_{T}^{rec} - p_{T}^{gen})/p_{T}^{gen}", kTH2F, {{1000, 0, 10}, {200, -1.0f, 1.0f}}, true);
      fRegistry.add("V0Leg/primary/hPtGen_DeltaEta", "electron #eta resolution;p_{T}^{gen} (GeV/c);#eta^{rec} - #eta^{gen}", kTH2F, {{1000, 0, 10}, {100, -0.5f, 0.5f}}, true);
      fRegistry.add("V0Leg/primary/hPtGen_DeltaPhi", "electron #varphi resolution;p_{T}^{gen} (GeV/c);#varphi^{rec} - #varphi^{gen} (rad.)", kTH2F, {{1000, 0, 10}, {100, -0.5f, 0.5f}}, true);
      fRegistry.add("V0Leg/primary/hRxyGen_DeltaPtOverPtGen", "electron p_{T} resolution; R_{xy}^{gen} (cm);(p_{T}^{rec} - p_{T}^{gen})/p_{T}^{gen}", kTH2F, {{100, 0, 100}, {200, -1.0f, 1.0f}}, true);
      fRegistry.add("V0Leg/primary/hRxyGen_DeltaEta", "electron #eta resolution;R_{xy}^{gen} (cm);#eta^{rec} - #eta^{gen}", kTH2F, {{100, 0, 100}, {100, -0.5f, 0.5f}}, true);
      fRegistry.add("V0Leg/primary/hRxyGen_DeltaPhi", "electron #varphi resolution;R_{xy}^{gen} (cm);#varphi^{rec} - #varphi^{gen} (rad.)", kTH2F, {{100, 0, 100}, {100, -0.5f, 0.5f}}, true);
      fRegistry.addClone("V0Leg/primary/", "V0Leg/fromWD/");
      fRegistry.addClone("V0Leg/primary/", "V0Leg/fromHS/");
      fRegistry.addClone("V0Leg/primary/", "V0Leg/fromPi0Dalitz/");
      fRegistry.addClone("V0Leg/primary/", "V0Leg/fromEtaDalitz/");
      fRegistry.addClone("V0Leg/primary/hPt", "V0Leg/candidate/hPt");         // only for purity
      fRegistry.addClone("V0Leg/primary/hEtaPhi", "V0Leg/candidate/hEtaPhi"); // only for purity

      // collision-association QA:
      if (mcAnalysisModeSettings.cfgDoCollisionAssociationQA) {
        fRegistry.add("CollisionAssociation/hDeltaCollId", "MC-event index difference;#Delta coll ID (rec - true);N_{#gamma}", kTH1F, {{21, -10.5, 10.5}}, true);
        fRegistry.add("CollisionAssociation/RightCollisions/hs", "correctly associated;R_{rec} (cm);p_{T} (GeV/c);#Delta p_{T}/p_{T,gen};#Delta z (cm);#Delta R (cm);#Delta coll ID", kTHnSparseF, {{100, 0, 100}, {200, 0, 10}, {200, -1, 1}, {200, -20, 20}, {200, -8, 8}, {21, -10.5, 10.5}}, true);
        fRegistry.add("CollisionAssociation/WrongCollisions/hs", "wrongly associated;R_{rec} (cm);p_{T} (GeV/c);#Delta p_{T}/p_{T,gen};#Delta z (cm);#Delta R (cm);#Delta coll ID", kTHnSparseF, {{100, 0, 100}, {200, 0, 10}, {200, -1, 1}, {200, -20, 20}, {200, -8, 8}, {21, -10.5, 10.5}}, true);
      }
    }

    // pT-dependent loss QA:
    if (qaSettingsGroup.cfgDoPtDependentLossQA) {
      fRegistry.add("V0/LossQA/before/hCosPAvsPt", "cosPA vs. pT;p_{T,#gamma} (GeV/c);cosine pointing angle", kTH2F, {{100, 0, 10}, {100, 0.99, 1.0}}, false);
      fRegistry.add("V0/LossQA/before/hPCAvsPt", "PCA vs. pT;p_{T,#gamma} (GeV/c);PCA (cm)", kTH2F, {{100, 0, 10}, {500, 0, 5}}, false);
      fRegistry.add("V0/LossQA/before/hQtAPvsPt", "AP q_{T} vs. pT;p_{T,#gamma} (GeV/c);q_{T} (GeV/c)", kTH2F, {{100, 0, 10}, {250, 0, 0.25}}, false);
      fRegistry.add("V0/LossQA/before/hRxyVsPt", "conversion radius vs. pT;p_{T,#gamma} (GeV/c);R_{xy} (cm)", kTH2F, {{100, 0, 10}, {200, 0, 100}}, false);
      fRegistry.add("V0/LossQA/before/hMeeVsPt", "m_{ee} vs. pT;p_{T,#gamma} (GeV/c);m_{ee} (GeV/c^{2})", kTH2F, {{100, 0, 10}, {100, 0, 0.1}}, false);
      fRegistry.addClone("V0/LossQA/before/", "V0/LossQA/after/");
    }

    if (doprocessRecoQA) {
      fRegistry.add("RecoQA/0_converted/hs", "true converted photons;p_{T,#gamma} (GeV/c);#eta;#varphi (rad.);R_{conv}^{true} (cm)", kTHnSparseF, {{100, 0, 10}, {36, -0.9, 0.9}, {36, 0, o2::constants::math::TwoPI}, {100, 0, 100}}, true);
      fRegistry.addClone("RecoQA/0_converted/hs", "RecoQA/1_bothLegsTracked/hs");
      fRegistry.addClone("RecoQA/0_converted/hs", "RecoQA/2_bothLegsInV0Legs/hs");
      fRegistry.addClone("RecoQA/0_converted/hs", "RecoQA/3_photonInV0PhotonsKF/hs");
      fRegistry.add("RecoQA/hStageVsPt", "reconstruction stage vs. pT;p_{T,#gamma} (GeV/c);stage (0=converted,1=tracked,2=V0Legs,3=V0PhotonsKF)", kTH2F, {{100, 0, 10}, {4, -0.5, 3.5}}, true);
      fRegistry.add("RecoQA/hMissingLeg_PtRconv", "legs tracked but lost before V0Legs;p_{T,leg}^{true} (GeV/c);R_{conv}^{true} (cm)", kTH2F, {{100, 0, 1}, {100, 0, 100}}, true);
      fRegistry.add("RecoQA/Duplicates/hNTracksPerTrueLeg", "reconstructed tracks per true conversion leg;N_{tracks} per true leg;N_{legs}", kTH1F, {{6, -0.5, 5.5}}, true);
      fRegistry.add("RecoQA/Duplicates/hNV0LegsPerTrueLeg", "V0Leg entries per true conversion leg;N_{V0Legs} per true leg;N_{legs}", kTH1F, {{6, -0.5, 5.5}}, true);
      fRegistry.add("RecoQA/Duplicates/hNV0CandsPerPhoton_Rconv", "correctly built V0 candidates per true photon;N_{cand} per true #gamma;R_{conv}^{true} (cm)", kTH2F, {{6, -0.5, 5.5}, {100, 0, 100}}, true);
      if (recoQASettingsGroup.cfgDoDetailedTrackQA) {
        fRegistry.add("RecoQA/LegQuality/survived/hNcrTPC", "TPC crossed rows;TPC crossed rows", kTH1F, {{161, -0.5, 160.5}}, true);
        fRegistry.add("RecoQA/LegQuality/survived/hChi2TPC", "chi2 per TPC cluster;#chi^{2}/N_{cls}^{TPC}", kTH1F, {{100, 0, 10}}, true);
        fRegistry.add("RecoQA/LegQuality/survived/hSharedFracTPC", "TPC shared-cluster fraction;N_{cls}^{shared}/N_{cls}", kTH1F, {{110, 0, 1.1}}, true);
        fRegistry.add("RecoQA/LegQuality/survived/hNclsTPC", "number of TPC clusters;N_{cls}^{TPC}", kTH1F, {{161, -0.5, 160.5}}, true);
        fRegistry.addClone("RecoQA/LegQuality/survived/", "RecoQA/LegQuality/lost/");
      }
    }
  }

  void DefineEMEventCut()
  {
    fEMEventCut = EMPhotonEventCut("fEMEventCut", "fEMEventCut");
    fEMEventCut.SetRequireSel8(eventcuts.cfgRequireSel8);
    fEMEventCut.SetRequireFT0AND(eventcuts.cfgRequireFT0AND);
    fEMEventCut.SetZvtxRange(-eventcuts.cfgZvtxMax, +eventcuts.cfgZvtxMax);
    fEMEventCut.SetRequireNoTFB(eventcuts.cfgRequireNoTFB);
    fEMEventCut.SetRequireNoITSROFB(eventcuts.cfgRequireNoITSROFB);
    fEMEventCut.SetRequireNoSameBunchPileup(eventcuts.cfgRequireNoSameBunchPileup);
    fEMEventCut.SetRequireVertexITSTPC(eventcuts.cfgRequireVertexITSTPC);
    fEMEventCut.SetRequireGoodZvtxFT0vsPV(eventcuts.cfgRequireGoodZvtxFT0vsPV);
    fEMEventCut.SetRequireNoCollInTimeRangeStandard(eventcuts.cfgRequireNoCollInTimeRangeStandard);
    fEMEventCut.SetRequireNoCollInTimeRangeStrict(eventcuts.cfgRequireNoCollInTimeRangeStrict);
    fEMEventCut.SetRequireNoCollInITSROFStandard(eventcuts.cfgRequireNoCollInITSROFStandard);
    fEMEventCut.SetRequireNoCollInITSROFStrict(eventcuts.cfgRequireNoCollInITSROFStrict);
    fEMEventCut.SetRequireNoHighMultCollInPrevRof(eventcuts.cfgRequireNoHighMultCollInPrevRof);
  }

  void DefinePCMCut()
  {
    fV0PhotonCut = V0PhotonCut("fV0PhotonCut", "fV0PhotonCut");

    // for v0
    fV0PhotonCut.SetV0PtRange(v0cuts.cfgMinPtV0, v0cuts.cfgMaxPtV0);
    fV0PhotonCut.SetV0EtaRange(v0cuts.cfgMinEtaV0, v0cuts.cfgMaxEtaV0);
    fV0PhotonCut.SetMinCosPA(v0cuts.cfgMinV0CosPA);
    fV0PhotonCut.SetMaxPCA(v0cuts.cfgMaxPCA);
    fV0PhotonCut.SetMaxChi2KF(v0cuts.cfgMaxChi2KF);
    fV0PhotonCut.SetRxyRange(v0cuts.cfgMinV0Radius, v0cuts.cfgMaxV0Radius, v0cuts.cfgMidLV0Radius, v0cuts.cfgMidHV0Radius);
    fV0PhotonCut.SetAPRange(v0cuts.cfgMaxAlphaAP, v0cuts.cfgMaxQtAP);
    fV0PhotonCut.RejectITSib(v0cuts.cfgRejectV0OnITSib);

    // for track
    fV0PhotonCut.SetMinNClustersTPC(trackcuts.cfgMinNClustersTPC);
    fV0PhotonCut.SetMinNCrossedRowsTPC(trackcuts.cfgMinNCrossedRows);
    fV0PhotonCut.SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    fV0PhotonCut.SetMaxFracSharedClustersTPC(trackcuts.cfgMaxFracSharedClustersTPC);
    fV0PhotonCut.SetChi2PerClusterTPC(0.0, trackcuts.cfgMaxChi2TPC);
    fV0PhotonCut.SetTPCNsigmaElRange(trackcuts.cfgMinTPCNsigmaEl, trackcuts.cfgMaxTPCNsigmaEl);
    fV0PhotonCut.SetChi2PerClusterITS(-1e+10, trackcuts.cfgMaxChi2ITS);
    fV0PhotonCut.SetNClustersITS(0, 7);
    fV0PhotonCut.SetMeanClusterSizeITSob(0.0, 16.0);
    fV0PhotonCut.SetDisableITSonly(trackcuts.cfgDisableITSonlyTracks);
    fV0PhotonCut.SetDisableTPConly(trackcuts.cfgDisableTPConlyTracks);
    fV0PhotonCut.SetRequireITSTPC(v0cuts.cfgRequireV0WithITSTPC);
    fV0PhotonCut.SetRequireITSonly(v0cuts.cfgRequireV0WithITSonly);
    fV0PhotonCut.SetRequireTPConly(v0cuts.cfgRequireV0WithTPConly);

    // for ML
    fV0PhotonCut.SetApplyMlCuts(mlcuts.cfgApplyPCMMl);
    fV0PhotonCut.SetUse2DBinning(mlcuts.cfgUse2DBinning);
    fV0PhotonCut.SetLoadMlModelsFromCCDB(mlcuts.cfgLoadModelsFromCCDB);
    fV0PhotonCut.SetNClassesMl(mlcuts.cfgNClassesPCMMl);
    fV0PhotonCut.SetMlTimestampCCDB(mlcuts.cfgTimestampCCDB);
    fV0PhotonCut.SetCcdbUrl(ccdburl);
    auto mCentralityTypeMlEnum = static_cast<CentType>(cfgCentEstimator.value);
    fV0PhotonCut.SetCentralityTypeMl(mCentralityTypeMlEnum);
    fV0PhotonCut.SetCutDirMl(mlcuts.cfgCutDirPCMMl);
    fV0PhotonCut.SetMlModelPathsCCDB(mlcuts.cfgModelPathsCCDB);
    fV0PhotonCut.SetMlOnnxFileNames(mlcuts.cfgOnnxFileNames);
    fV0PhotonCut.SetBinsPtMl(mlcuts.cfgBinsPtPCMMl);
    fV0PhotonCut.SetBinsCentMl(mlcuts.cfgBinsCentPCMMl);
    fV0PhotonCut.SetCutsMl(mlcuts.cfgCutsPCMMlFlat);
    fV0PhotonCut.SetNamesInputFeatures(mlcuts.cfgNamesInputFeatures);
    fV0PhotonCut.SetLabelsBinsMl(mlcuts.cfgLabelsBinsPCMMl);
    fV0PhotonCut.SetLabelsCutScoresMl(mlcuts.cfgLabelsCutScoresPCMMl);
    fV0PhotonCut.SetD_Bz(0.0f);

    if (mlcuts.cfgApplyPCMMl) {
      fV0PhotonCut.initV0MlModels(ccdbApi);
    }
  }

  template <const int ev_id, typename TCollision>
  void fillEventInfo(TCollision const& collision, const float /*weight*/ = 1.f)
  {
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 1.0);
    if (collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 2.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 3.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 4.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 5.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 6.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 7.0);
    }
    if (collision.sel8()) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 8.0);
    }
    if (std::fabs(collision.posZ()) < 10.0) { // o2-linter: disable=magic-number (vertex cut)
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 9.0);
    }
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hZvtx"), collision.posZ());

    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultNTracksPV"), collision.multNTracksPV());
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultNTracksPVeta1"), collision.multNTracksPVeta1());
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultFT0"), collision.multFT0A(), collision.multFT0C());
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCentFT0A"), collision.centFT0A());
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCentFT0C"), collision.centFT0C());
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCentFT0M"), collision.centFT0M());
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCentFT0MvsMultNTracksPV"), collision.centFT0M(), collision.multNTracksPV());
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultFT0MvsMultNTracksPV"), collision.multFT0A() + collision.multFT0C(), collision.multNTracksPV());
  }

  template <typename TV0>
  void fillV0Info(TV0 const& v0)
  {
    fRegistry.fill(HIST("V0/hPt"), v0.pt());
    fRegistry.fill(HIST("V0/hEtaPhi"), v0.phi(), v0.eta());
    fRegistry.fill(HIST("V0/hXY"), v0.vx(), v0.vy());
    fRegistry.fill(HIST("V0/hRZ"), v0.vz(), v0.v0radius());
    fRegistry.fill(HIST("V0/hCosPA"), v0.cospa());
    fRegistry.fill(HIST("V0/hCosPAXY"), v0.cospaXY());
    fRegistry.fill(HIST("V0/hCosPARZ"), v0.cospaRZ());
    fRegistry.fill(HIST("V0/hPCA"), v0.pca());
    fRegistry.fill(HIST("V0/hDCAxyz"), v0.dcaXYtopv(), v0.dcaZtopv());
    fRegistry.fill(HIST("V0/hDCAz_Pt"), v0.dcaZtopv(), v0.pt());
    fRegistry.fill(HIST("V0/hAPplot"), v0.alpha(), v0.qtarm());
    fRegistry.fill(HIST("V0/hRxyVsPt"), v0.pt(), v0.v0radius());
    fRegistry.fill(HIST("V0/hEtaVsPt"), v0.pt(), v0.eta());
    fRegistry.fill(HIST("V0/hPhiVsPt"), v0.pt(), v0.phi());
    fRegistry.fill(HIST("V0/hsAlphaQtPt"), v0.alpha(), v0.qtarm(), v0.pt());
    if constexpr (requires { v0.psipair(); }) {
      fRegistry.fill(HIST("V0/hPsiPairVsPt"), v0.pt(), v0.psipair());
    }
    fRegistry.fill(HIST("V0/hMassGamma"), v0.v0radius(), v0.mGamma());
    fRegistry.fill(HIST("V0/hKFChi2vsM"), v0.mGamma(), v0.chiSquareNDF());
    fRegistry.fill(HIST("V0/hKFChi2vsR"), v0.v0radius(), v0.chiSquareNDF());
    fRegistry.fill(HIST("V0/hKFChi2vsX"), v0.vx(), v0.chiSquareNDF());
    fRegistry.fill(HIST("V0/hKFChi2vsY"), v0.vy(), v0.chiSquareNDF());
    fRegistry.fill(HIST("V0/hKFChi2vsZ"), v0.vz(), v0.chiSquareNDF());

    float phi_cp = std::atan2(v0.vy(), v0.vx());
    o2::math_utils::bringTo02Pi(phi_cp);
    float eta_cp = std::atanh(v0.vz() / std::sqrt(std::pow(v0.vx(), 2) + std::pow(v0.vy(), 2) + std::pow(v0.vz(), 2)));
    fRegistry.fill(HIST("V0/hsConvPoint"), v0.v0radius(), phi_cp, eta_cp);

    // BDT response histogram can be filled here when apply BDT is true
    if (mlcuts.cfgApplyPCMMl) {
      const std::span<const float>& bdtValue = fV0PhotonCut.getBDTValue();
      float psipair = 999.f;
      float phiv = 999.f;
      if constexpr (requires { v0.psipair(); v0.phiv(); }) {
        psipair = v0.psipair();
        phiv = v0.phiv();
      }
      fRegistry.fill(HIST("V0/hPhiVPsi"), psipair, phiv);
      if (mlcuts.cfgNClassesPCMMl == 2 && bdtValue.size() == 2) { // o2-linter: disable=magic-number (BDT)
        fRegistry.fill(HIST("V0/hBDTBackgroundScoreVsPt"), v0.pt(), bdtValue[0]);
        fRegistry.fill(HIST("V0/hBDTSignalScoreVsPt"), v0.pt(), bdtValue[1]);
      } else if (mlcuts.cfgNClassesPCMMl == 3 && bdtValue.size() == 3) { // o2-linter: disable=magic-number (BDT)
        fRegistry.fill(HIST("V0/hBDTBackgroundScoreVsPt"), v0.pt(), bdtValue[0]);
        fRegistry.fill(HIST("V0/hBDTPrimaryPhotonScoreVsPt"), v0.pt(), bdtValue[1]);
        fRegistry.fill(HIST("V0/hBDTSecondaryPhotonScoreVsPt"), v0.pt(), bdtValue[2]);
      } else if (bdtValue.size() == 1) {
        fRegistry.fill(HIST("V0/hBDTScoreVsPt"), v0.pt(), bdtValue[0]);
      }
    }
  }

  template <typename TLeg>
  void fillV0LegInfo(TLeg const& leg)
  {
    fRegistry.fill(HIST("V0Leg/hPt"), leg.pt());
    fRegistry.fill(HIST("V0Leg/hQoverPt"), leg.sign() / leg.pt());
    fRegistry.fill(HIST("V0Leg/hEtaPhi"), leg.phi(), leg.eta());
    fRegistry.fill(HIST("V0Leg/hDCAxyz"), leg.dcaXY(), leg.dcaZ());
    fRegistry.fill(HIST("V0Leg/hNclsITS"), leg.itsNCls());
    fRegistry.fill(HIST("V0Leg/hNclsTPC"), leg.tpcNClsFound());
    fRegistry.fill(HIST("V0Leg/hNcrTPC"), leg.tpcNClsCrossedRows());
    fRegistry.fill(HIST("V0Leg/hTPCNcr2Nf"), leg.tpcCrossedRowsOverFindableCls());
    fRegistry.fill(HIST("V0Leg/hTPCNcls2Nf"), leg.tpcFoundOverFindableCls());
    fRegistry.fill(HIST("V0Leg/hTPCNclsShared"), leg.pt(), leg.tpcFractionSharedCls());
    fRegistry.fill(HIST("V0Leg/hChi2TPC"), leg.tpcChi2NCl());
    fRegistry.fill(HIST("V0Leg/hChi2ITS"), leg.itsChi2NCl());
    fRegistry.fill(HIST("V0Leg/hITSClusterMap"), leg.itsClusterMap());
    if (leg.hasITS()) {
      fRegistry.fill(HIST("V0Leg/hMeanClusterSizeITS"), leg.p(), leg.meanClusterSizeITS() * std::cos(std::atan(leg.tgl())));
    }
    fRegistry.fill(HIST("V0Leg/hTPCdEdx"), leg.tpcInnerParam(), leg.tpcSignal());
    fRegistry.fill(HIST("V0Leg/hTPCNsigmaEl"), leg.tpcInnerParam(), leg.tpcNSigmaEl());
    fRegistry.fill(HIST("V0Leg/hTPCNsigmaPi"), leg.tpcInnerParam(), leg.tpcNSigmaPi());
    fRegistry.fill(HIST("V0Leg/hTPCNsigmaElVsEta"), leg.eta(), leg.tpcNSigmaEl());
  }

  template <const int ev_id, typename TV0>
  void fillLossQAInfo(TV0 const& v0)
  {
    if (!qaSettingsGroup.cfgDoPtDependentLossQA) {
      return;
    }
    fRegistry.fill(HIST("V0/LossQA/") + HIST(event_types[ev_id]) + HIST("hCosPAvsPt"), v0.pt(), v0.cospa());
    fRegistry.fill(HIST("V0/LossQA/") + HIST(event_types[ev_id]) + HIST("hPCAvsPt"), v0.pt(), v0.pca());
    fRegistry.fill(HIST("V0/LossQA/") + HIST(event_types[ev_id]) + HIST("hQtAPvsPt"), v0.pt(), v0.qtarm());
    fRegistry.fill(HIST("V0/LossQA/") + HIST(event_types[ev_id]) + HIST("hRxyVsPt"), v0.pt(), v0.v0radius());
    fRegistry.fill(HIST("V0/LossQA/") + HIST(event_types[ev_id]) + HIST("hMeeVsPt"), v0.pt(), v0.mGamma());
  }

  template <typename TV0>
  void fillMaterialBudgetInfo(TV0 const& v0)
  {
    if (materialBudgetSettingsGroup.cfgDoMaterialDistribution) {
      fRegistry.fill(HIST("MaterialBudget/hs"), v0.vz(), v0.v0radius(), v0.eta(), v0.phi(), v0.pt());
    }
    if (materialBudgetSettingsGroup.cfgDoWiresDetail && v0.v0radius() < 20.f && std::fabs(v0.vz()) < 20.f) { // o2-linter: disable=magic-number (region window = booking range)
      if (3.15f < v0.phi() && v0.phi() < 3.4f) {                                                             // o2-linter: disable=magic-number (region window = booking range)
        fRegistry.fill(HIST("MaterialBudget/Wires/hsLeft"), v0.vx(), v0.vy(), v0.vz(), v0.phi(), v0.v0radius(), v0.pt());
      } else if (6.00f < v0.phi() && v0.phi() < 6.15f) { // o2-linter: disable=magic-number (region window = booking range)
        fRegistry.fill(HIST("MaterialBudget/Wires/hsRight"), v0.vx(), v0.vy(), v0.vz(), v0.phi(), v0.v0radius(), v0.pt());
      }
    }
    if (materialBudgetSettingsGroup.cfgDoITSDetail && v0.v0radius() < 60.f && std::fabs(v0.vz()) < 40.f) { // o2-linter: disable=magic-number (region window = booking range)
      fRegistry.fill(HIST("MaterialBudget/ITS/hs"), v0.vx(), v0.vy(), v0.vz(), v0.phi(), v0.v0radius(), v0.pt());
    }
    if (materialBudgetSettingsGroup.cfgDoMFTDetail && 40.f < v0.v0radius() && v0.v0radius() < 60.f && std::fabs(v0.vz()) < 40.f) { // o2-linter: disable=magic-number (region window = booking range)
      fRegistry.fill(HIST("MaterialBudget/MFT/hs"), v0.vx(), v0.vy(), v0.vz(), v0.phi(), v0.v0radius(), v0.pt());
    }
    if (materialBudgetSettingsGroup.cfgDoTPCInnerBarrelDetail && 60.f < v0.v0radius() && v0.v0radius() < 80.f && std::fabs(v0.vz()) < 40.f) { // o2-linter: disable=magic-number (region window = booking range)
      fRegistry.fill(HIST("MaterialBudget/TPC/hs"), v0.vx(), v0.vy(), v0.vz(), v0.phi(), v0.v0radius(), v0.pt());
    }
  }

  Preslice<MyV0Photons> perCollisionV0 = aod::v0photonkf::pmeventId;
  Preslice<MyV0PhotonsML> perCollisionV0ML = aod::v0photonkf::pmeventId;
  Filter collisionFilter_centrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  Filter collisionFilter_occupancy_track = eventcuts.cfgTrackOccupancyMin <= o2::aod::evsel::trackOccupancyInTimeRange && o2::aod::evsel::trackOccupancyInTimeRange < eventcuts.cfgTrackOccupancyMax;
  Filter collisionFilter_occupancy_ft0c = eventcuts.cfgFT0COccupancyMin <= o2::aod::evsel::ft0cOccupancyInTimeRange && o2::aod::evsel::ft0cOccupancyInTimeRange < eventcuts.cfgFT0COccupancyMax;
  using FilteredMyCollisions = soa::Filtered<MyCollisions>;
  using FilteredMyCollisionsMC = soa::Filtered<MyCollisionsMC>; // same filters, they act column-wise

  template <typename TV0Photon, typename TPerCollision>
  void process(FilteredMyCollisions const& collisions, TV0Photon const& v0photons, aod::V0Legs const&, TPerCollision const& perCollision)
  {
    for (const auto& collision : collisions) {
      initCCDB(collision);
      const std::array<float, 3> centralities = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      fillEventInfo<0>(collision);
      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      fillEventInfo<1>(collision);
      fRegistry.fill(HIST("Event/before/hCollisionCounter"), 10.0); // accepted
      fRegistry.fill(HIST("Event/after/hCollisionCounter"), 10.0);  // accepted

      fV0PhotonCut.SetCentrality(centralities[cfgCentEstimator]);
      int nv0 = 0;
      auto v0photons_coll = v0photons.sliceBy(perCollision, collision.globalIndex());
      for (const auto& v0 : v0photons_coll) {
        auto pos = v0.template posTrack_as<aod::V0Legs>();
        auto ele = v0.template negTrack_as<aod::V0Legs>();
        fillLossQAInfo<0>(v0); // all candidates offered to the selection
        if (!fV0PhotonCut.IsSelected<decltype(v0), aod::V0Legs>(v0)) {
          continue;
        }
        fillLossQAInfo<1>(v0);
        fillV0Info(v0);
        fillMaterialBudgetInfo(v0);
        for (const auto& leg : {pos, ele}) {
          fillV0LegInfo(leg);
        }
        if (trackcuts.cfgDoDEdxPostCalibration) {
          fRegistry.fill(HIST("V0Leg/hPvsConvPointvsTPCNsigmaElvsEta_Pos"), pos.p(), v0.v0radius(), pos.tpcNSigmaEl(), pos.eta());
          fRegistry.fill(HIST("V0Leg/hPvsConvPointvsTPCNsigmaElvsEta_Ele"), ele.p(), v0.v0radius(), ele.tpcNSigmaEl(), ele.eta());
        }
        nv0++;
      } // end of v0 loop
      fRegistry.fill(HIST("V0/hNgamma"), nv0);
    } // end of collision loop
  }
  void processQC(FilteredMyCollisions const& collisions, MyV0Photons const& v0photons, aod::V0Legs const& v0legs)
  {
    process(collisions, v0photons, v0legs, perCollisionV0);
  } // end of process

  void processQCML(FilteredMyCollisions const& collisions, MyV0PhotonsML const& v0photonsML, aod::V0Legs const& v0legs)
  {
    process(collisions, v0photonsML, v0legs, perCollisionV0ML);
  } // end of ML process

  template <int mctype, typename TV0, typename TMCV0, typename TMCLeg>
  void fillV0InfoMC(TV0 const& v0, TMCV0 const& mcphoton, TMCLeg const& mcleg)
  {
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hPt"), v0.pt());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hEtaPhi"), v0.phi(), v0.eta());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hXY"), v0.vx(), v0.vy());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hRZ"), v0.vz(), v0.v0radius());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hCosPA"), v0.cospa());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hCosPAXY"), v0.cospaXY());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hCosPARZ"), v0.cospaRZ());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hPCA"), v0.pca());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hDCAxyz"), v0.dcaXYtopv(), v0.dcaZtopv());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hDCAz_Pt"), v0.dcaZtopv(), v0.pt());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hAPplot"), v0.alpha(), v0.qtarm());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hRxyVsPt"), v0.pt(), v0.v0radius());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hEtaVsPt"), v0.pt(), v0.eta());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hPhiVsPt"), v0.pt(), v0.phi());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hsAlphaQtPt"), v0.alpha(), v0.qtarm(), v0.pt());
    if constexpr (requires { v0.psipair(); }) { // only the ML join carries psi_pair
      fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hPsiPairVsPt"), v0.pt(), v0.psipair());
    }
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hMassGamma"), v0.v0radius(), v0.mGamma());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hKFChi2vsM"), v0.mGamma(), v0.chiSquareNDF());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hKFChi2vsR"), v0.v0radius(), v0.chiSquareNDF());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hKFChi2vsX"), v0.vx(), v0.chiSquareNDF());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hKFChi2vsY"), v0.vy(), v0.chiSquareNDF());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hKFChi2vsZ"), v0.vz(), v0.chiSquareNDF());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hPtGen_DeltaPtOverPtGen"), mcphoton.pt(), (v0.pt() - mcphoton.pt()) / mcphoton.pt());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hPtGen_DeltaEta"), mcphoton.pt(), v0.eta() - mcphoton.eta());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hPtGen_DeltaPhi"), mcphoton.pt(), v0.phi() - mcphoton.phi());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hRxyGen_DeltaPtOverPtGen"), std::sqrt(std::pow(mcleg.vx(), 2) + std::pow(mcleg.vy(), 2)), (v0.pt() - mcphoton.pt()) / mcphoton.pt());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hRxyGen_DeltaEta"), std::sqrt(std::pow(mcleg.vx(), 2) + std::pow(mcleg.vy(), 2)), v0.eta() - mcphoton.eta());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hRxyGen_DeltaPhi"), std::sqrt(std::pow(mcleg.vx(), 2) + std::pow(mcleg.vy(), 2)), v0.phi() - mcphoton.phi());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hRxyGen_DeltaR"), std::sqrt(std::pow(mcleg.vx(), 2) + std::pow(mcleg.vy(), 2)), v0.v0radius() - std::sqrt(std::pow(mcleg.vx(), 2) + std::pow(mcleg.vy(), 2)));
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hConvPoint_diffX"), mcleg.vx(), v0.vx() - mcleg.vx());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hConvPoint_diffY"), mcleg.vy(), v0.vy() - mcleg.vy());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hConvPoint_diffZ"), mcleg.vz(), v0.vz() - mcleg.vz());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hXY_MC"), mcleg.vx(), mcleg.vy());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hRZ_MC"), mcleg.vz(), std::sqrt(std::pow(mcleg.vx(), 2) + std::pow(mcleg.vy(), 2)));

    if (mcAnalysisModeSettings.cfgDoDetailedResolution) {
      float resolPt = (v0.pt() - mcphoton.pt()) / mcphoton.pt();
      fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hPtResolPtEtaPhi"), resolPt, mcphoton.pt(), mcphoton.eta(), mcphoton.phi());
      float resolMom = (v0.p() - mcphoton.p()) / mcphoton.p();
      fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hMomResolPtEtaPhi"), resolMom, mcphoton.pt(), mcphoton.eta(), mcphoton.phi());
      fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hPtEtaPhi"), mcphoton.pt(), mcphoton.eta(), mcphoton.phi());
      fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hEtaPhiResol"), mcphoton.pt(), v0.eta() - mcphoton.eta(), v0.phi() - mcphoton.phi());
    }

    float phi_cp = std::atan2(v0.vy(), v0.vx());
    o2::math_utils::bringTo02Pi(phi_cp);
    float eta_cp = std::atanh(v0.vz() / std::sqrt(std::pow(v0.vx(), 2) + std::pow(v0.vy(), 2) + std::pow(v0.vz(), 2)));
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hsConvPoint"), v0.v0radius(), phi_cp, eta_cp);

    if (mlcuts.cfgApplyPCMMl) {
      const std::span<const float>& bdtValue = fV0PhotonCut.getBDTValue();
      float psipair = 999.f;
      float phiv = 999.f;
      if constexpr (requires { v0.psipair(); v0.phiv(); }) {
        psipair = v0.psipair();
        phiv = v0.phiv();
      }
      fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hPhiVPsi"), psipair, phiv);
      if (mlcuts.cfgNClassesPCMMl == 2 && bdtValue.size() == 2) { // o2-linter: disable=magic-number (BDT)
        fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hBDTBackgroundScoreVsPt"), v0.pt(), bdtValue[0]);
        fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hBDTSignalScoreVsPt"), v0.pt(), bdtValue[1]);
      } else if (mlcuts.cfgNClassesPCMMl == 3 && bdtValue.size() == 3) { // o2-linter: disable=magic-number (BDT)
        fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hBDTBackgroundScoreVsPt"), v0.pt(), bdtValue[0]);
        fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hBDTPrimaryPhotonScoreVsPt"), v0.pt(), bdtValue[1]);
        fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hBDTSecondaryPhotonScoreVsPt"), v0.pt(), bdtValue[2]);
      } else if (bdtValue.size() == 1) { // o2-linter: disable=magic-number (BDT)
        fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hBDTScoreVsPt"), v0.pt(), bdtValue[0]);
      }
    }
  }

  template <int mctype, typename TLeg>
  void fillV0LegInfoMC(TLeg const& leg)
  {
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hPt"), leg.pt());
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hQoverPt"), leg.sign() / leg.pt());
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hEtaPhi"), leg.phi(), leg.eta());
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hDCAxyz"), leg.dcaXY(), leg.dcaZ());
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hNclsITS"), leg.itsNCls());
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hNclsTPC"), leg.tpcNClsFound());
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hNcrTPC"), leg.tpcNClsCrossedRows());
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hTPCNcr2Nf"), leg.tpcCrossedRowsOverFindableCls());
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hTPCNcls2Nf"), leg.tpcFoundOverFindableCls());
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hTPCNclsShared"), leg.pt(), leg.tpcFractionSharedCls());
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hChi2TPC"), leg.tpcChi2NCl());
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hChi2ITS"), leg.itsChi2NCl());
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hITSClusterMap"), leg.itsClusterMap());
    if (leg.hasITS()) {
      fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hMeanClusterSizeITS"), leg.p(), leg.meanClusterSizeITS() * std::cos(std::atan(leg.tgl())));
    }
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hTPCdEdx"), leg.tpcInnerParam(), leg.tpcSignal());
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hTPCNsigmaEl"), leg.tpcInnerParam(), leg.tpcNSigmaEl());
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hTPCNsigmaPi"), leg.tpcInnerParam(), leg.tpcNSigmaPi());
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hTPCNsigmaElVsEta"), leg.eta(), leg.tpcNSigmaEl());
    auto mcleg = leg.template emmcparticle_as<aod::EMMCParticles>();
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hPtGen_DeltaPtOverPtGen"), mcleg.pt(), (leg.pt() - mcleg.pt()) / mcleg.pt());
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hPtGen_DeltaEta"), mcleg.pt(), leg.eta() - mcleg.eta());
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hPtGen_DeltaPhi"), mcleg.pt(), leg.phi() - mcleg.phi());
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hRxyGen_DeltaPtOverPtGen"), std::sqrt(std::pow(mcleg.vx(), 2) + std::pow(mcleg.vy(), 2)), (leg.pt() - mcleg.pt()) / mcleg.pt());
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hRxyGen_DeltaEta"), std::sqrt(std::pow(mcleg.vx(), 2) + std::pow(mcleg.vy(), 2)), leg.eta() - mcleg.eta());
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hRxyGen_DeltaPhi"), std::sqrt(std::pow(mcleg.vx(), 2) + std::pow(mcleg.vy(), 2)), leg.phi() - mcleg.phi());
  }

  template <typename TV0Photons, typename TPerCollision>
  void processMC(FilteredMyCollisionsMC const& collisions, TV0Photons const& v0photons, aod::EMMCParticles const& mcparticles, MyMCV0Legs const&, aod::EMMCEvents const&, TPerCollision const& percollision)
  {
    for (const auto& collision : collisions) {
      initCCDB(collision);
      const std::array<float, 3> centralities = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      fillEventInfo<0>(collision);
      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      fillEventInfo<1>(collision);
      fRegistry.fill(HIST("Event/before/hCollisionCounter"), 10.0); // accepted
      fRegistry.fill(HIST("Event/after/hCollisionCounter"), 10.0);  // accepted

      fV0PhotonCut.SetCentrality(centralities[cfgCentEstimator]); // set centrality for BDT response
      auto v0photons_coll = v0photons.sliceBy(percollision, collision.globalIndex());
      int ng_primary = 0, ng_wd = 0, ng_hs = 0, nee_pi0 = 0, nee_eta = 0;
      for (const auto& v0 : v0photons_coll) {
        auto pos = v0.template posTrack_as<MyMCV0Legs>();
        auto ele = v0.template negTrack_as<MyMCV0Legs>();
        auto posmc = pos.template emmcparticle_as<aod::EMMCParticles>();
        auto elemc = ele.template emmcparticle_as<aod::EMMCParticles>();

        fillLossQAInfo<0>(v0); // all candidates offered to the selection
        if (!fV0PhotonCut.IsSelected<decltype(v0), MyMCV0Legs>(v0)) {
          continue;
        }
        fillLossQAInfo<1>(v0);

        fillMaterialBudgetInfo(v0);

        fRegistry.fill(HIST("V0/candidate/hPt"), v0.pt());
        fRegistry.fill(HIST("V0/candidate/hEtaPhi"), v0.phi(), v0.eta());
        for (const auto& leg : {pos, ele}) {
          fRegistry.fill(HIST("V0Leg/candidate/hPt"), leg.pt());
          fRegistry.fill(HIST("V0Leg/candidate/hEtaPhi"), leg.phi(), leg.eta());
        }

        int photonid = FindCommonMotherFrom2Prongs(posmc, elemc, -11, 11, 22, mcparticles);
        int pi0id = FindCommonMotherFrom2Prongs(posmc, elemc, -11, 11, 111, mcparticles); // pi0 dalitz decay
        int etaid = FindCommonMotherFrom2Prongs(posmc, elemc, -11, 11, 221, mcparticles); // eta dalitz decay
        if (photonid < 0 && pi0id < 0 && etaid < 0) {
          continue;
        }

        if (photonid > 0) {
          auto mcphoton = mcparticles.iteratorAt(photonid);

          if (mcAnalysisModeSettings.cfgDoCollisionAssociationQA) {
            const int deltaCol = static_cast<int>(collision.emmceventId()) - static_cast<int>(mcphoton.emmceventId());
            const float rGen = std::sqrt(std::pow(elemc.vx(), 2) + std::pow(elemc.vy(), 2));
            const float deltaR = v0.v0radius() - rGen;
            const float deltaZ = v0.vz() - elemc.vz();
            const float relPtRes = mcphoton.pt() > 0.f ? (v0.pt() - mcphoton.pt()) / mcphoton.pt() : 0.f;
            fRegistry.fill(HIST("CollisionAssociation/hDeltaCollId"), deltaCol);
            if (deltaCol == 0) {
              fRegistry.fill(HIST("CollisionAssociation/RightCollisions/hs"), v0.v0radius(), v0.pt(), relPtRes, deltaZ, deltaR, deltaCol);
            } else {
              fRegistry.fill(HIST("CollisionAssociation/WrongCollisions/hs"), v0.v0radius(), v0.pt(), relPtRes, deltaZ, deltaR, deltaCol);
            }
          }

          if (mcAnalysisModeSettings.cfgRequireTrueAssociation && (mcphoton.emmceventId() != collision.emmceventId())) {
            continue;
          }

          if (mcphoton.isPhysicalPrimary() || mcphoton.producedByGenerator()) {
            fillV0InfoMC<0>(v0, mcphoton, elemc);
            for (const auto& leg : {pos, ele}) {
              fillV0LegInfoMC<0>(leg);
            }
            ng_primary++;
          } else if (IsFromWD(mcphoton.template emmcevent_as<aod::EMMCEvents>(), mcphoton, mcparticles) > 0) {
            fillV0InfoMC<1>(v0, mcphoton, elemc);
            for (const auto& leg : {pos, ele}) {
              fillV0LegInfoMC<1>(leg);
            }
            ng_wd++;
          } else {
            fillV0InfoMC<2>(v0, mcphoton, elemc);
            for (const auto& leg : {pos, ele}) {
              fillV0LegInfoMC<2>(leg);
            }
            ng_hs++;
          }
        } else if (pi0id > 0) {
          auto mcpi0 = mcparticles.iteratorAt(pi0id);
          if (mcAnalysisModeSettings.cfgRequireTrueAssociation && (mcpi0.emmceventId() != collision.emmceventId())) {
            continue;
          }
          if (mcpi0.isPhysicalPrimary() || mcpi0.producedByGenerator()) {
            fillV0InfoMC<3>(v0, mcpi0, elemc);
            for (const auto& leg : {pos, ele}) {
              fillV0LegInfoMC<3>(leg);
            }
            nee_pi0++;
          }
        } else if (etaid > 0) {
          auto mceta = mcparticles.iteratorAt(etaid);
          if (mcAnalysisModeSettings.cfgRequireTrueAssociation && (mceta.emmceventId() != collision.emmceventId())) {
            continue;
          }
          if (mceta.isPhysicalPrimary() || mceta.producedByGenerator()) {
            fillV0InfoMC<4>(v0, mceta, elemc);
            for (const auto& leg : {pos, ele}) {
              fillV0LegInfoMC<4>(leg);
            }
            nee_eta++;
          }
        }
      } // end of v0 loop
      fRegistry.fill(HIST("V0/primary/hNgamma"), ng_primary);
      fRegistry.fill(HIST("V0/fromWD/hNgamma"), ng_wd);
      fRegistry.fill(HIST("V0/fromHS/hNgamma"), ng_hs);
      fRegistry.fill(HIST("V0/fromPi0Dalitz/hNgamma"), nee_pi0);
      fRegistry.fill(HIST("V0/fromEtaDalitz/hNgamma"), nee_eta);
    } // end of collision loop
  } // end of MC process

  void processPCMQCMC(FilteredMyCollisionsMC const& collisions, MyV0Photons const& v0photons, aod::EMMCParticles const& mcparticles, MyMCV0Legs const& mcv0legs, aod::EMMCEvents const& mcevents)
  {
    processMC(collisions, v0photons, mcparticles, mcv0legs, mcevents, perCollisionV0);
  } // end of MC QC process

  void processPCMQCMCML(FilteredMyCollisionsMC const& collisions, MyV0PhotonsML const& v0photonsML, aod::EMMCParticles const& mcparticles, MyMCV0Legs const& mcv0legs, aod::EMMCEvents const& mcevents)
  {
    processMC(collisions, v0photonsML, mcparticles, mcv0legs, mcevents, perCollisionV0ML);
  } // end of MC QC process with ML cuts

  template <typename TBinnedData>
  void fillBinnedData(TBinnedData const& binned_data, const float weight = 1.f)
  {
    int xbin = 0, ybin = 0, zbin = 0;
    auto hPtY = fRegistry.get<TH2>(HIST("Generated/hPtY")); // 2D
    auto hPt = fRegistry.get<TH1>(HIST("Generated/hPt"));   // 1D

    for (int ibin = 0; ibin < hPtY->GetNcells(); ibin++) {
      int nentry = binned_data[ibin];
      hPtY->GetBinXYZ(ibin, xbin, ybin, zbin);
      float pt = hPtY->GetXaxis()->GetBinCenter(xbin);
      float y = hPtY->GetYaxis()->GetBinCenter(ybin);
      if (y > v0cuts.cfgMaxEtaV0) {
        continue;
      }

      for (int j = 0; j < nentry; j++) {
        hPtY->Fill(pt, y, weight);
        hPt->Fill(pt, weight);
      }
    }
  }

  PresliceUnsorted<aod::EMMCParticles> perMcCollision = aod::emmcparticle::emmceventId;
  void processGen(FilteredMyCollisionsMC const& collisions, MyMCCollisions const&, aod::EMMCParticles const& mcparticles)
  {

    for (const auto& collision : collisions) {
      const std::array<float, 3> centralities = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }
      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }

      auto mccollision = collision.template emmcevent_as<MyMCCollisions>();
      auto binned_data_gamma_gen = mccollision.generatedGamma();
      fillBinnedData(binned_data_gamma_gen, 1.f);

      auto mctracks_coll = mcparticles.sliceBy(perMcCollision, mccollision.globalIndex());
      for (const auto& mctrack : mctracks_coll) {
        if (std::fabs(mctrack.y()) > v0cuts.cfgMaxEtaV0) {
          continue;
        }
        if (std::abs(mctrack.pdgCode()) != PDG_t::kGamma || !(mctrack.isPhysicalPrimary() || mctrack.producedByGenerator())) {
          continue;
        }

        if (mcAnalysisModeSettings.cfgDoDetailedResolution) {
          fRegistry.fill(HIST("Generated/hPtEtaPhi"), mctrack.pt(), mctrack.eta(), mctrack.phi()); // for all generated photons before any kinematic cut
        }
        if (mctrack.daughtersIds().empty()) {
          continue;
        }
        auto daughter = mcparticles.iteratorAt(mctrack.daughtersIds()[0]); // choose ele or pos.
        float rxy_gen_e = std::sqrt(std::pow(daughter.vx(), 2) + std::pow(daughter.vy(), 2));
        float phi_cp = std::atan2(daughter.vy(), daughter.vx());
        o2::math_utils::bringTo02Pi(phi_cp);
        float eta_cp = std::atanh(daughter.vz() / std::sqrt(std::pow(daughter.vx(), 2) + std::pow(daughter.vy(), 2) + std::pow(daughter.vz(), 2)));

        fRegistry.fill(HIST("Generated/hR_ConversionPhoton_wideR"), rxy_gen_e);
        fRegistry.fill(HIST("Generated/hRZ_wideR"), daughter.vz(), rxy_gen_e);
        fRegistry.fill(HIST("Generated/hsRZPhi_wideR"), daughter.vz(), phi_cp, rxy_gen_e);
        if (rxy_gen_e < std::fabs(daughter.vz()) * std::tan(2 * std::atan(std::exp(-v0cuts.cfgMaxEtaV0))) - genSettingsGroup.cfgMarginZMC) {
          continue;
        }
        if (rxy_gen_e > genSettingsGroup.cfgMaxRGen) {
          continue;
        }

        fRegistry.fill(HIST("Generated/hPt_ConversionPhoton"), mctrack.pt());
        fRegistry.fill(HIST("Generated/hY_ConversionPhoton"), mctrack.y());
        fRegistry.fill(HIST("Generated/hPhi_ConversionPhoton"), mctrack.phi());
        fRegistry.fill(HIST("Generated/hsConvPoint"), rxy_gen_e, phi_cp, eta_cp);
        fRegistry.fill(HIST("Generated/hXY"), daughter.vx(), daughter.vy());
        fRegistry.fill(HIST("Generated/hRZ"), daughter.vz(), rxy_gen_e);
        fRegistry.fill(HIST("Generated/hRPhi"), phi_cp, rxy_gen_e);
      } // end of mctrack loop per collision
    } // end of collision loop
  }

  template <bool survived, typename TTrack>
  void fillLegQualityRecoQA(TTrack const& track)
  {
    if constexpr (survived) {
      fRegistry.fill(HIST("RecoQA/LegQuality/survived/hNcrTPC"), track.tpcNClsCrossedRows());
      fRegistry.fill(HIST("RecoQA/LegQuality/survived/hChi2TPC"), track.tpcChi2NCl());
      fRegistry.fill(HIST("RecoQA/LegQuality/survived/hSharedFracTPC"), track.tpcFractionSharedCls());
      fRegistry.fill(HIST("RecoQA/LegQuality/survived/hNclsTPC"), track.tpcNClsFound());
    } else {
      fRegistry.fill(HIST("RecoQA/LegQuality/lost/hNcrTPC"), track.tpcNClsCrossedRows());
      fRegistry.fill(HIST("RecoQA/LegQuality/lost/hChi2TPC"), track.tpcChi2NCl());
      fRegistry.fill(HIST("RecoQA/LegQuality/lost/hSharedFracTPC"), track.tpcFractionSharedCls());
      fRegistry.fill(HIST("RecoQA/LegQuality/lost/hNclsTPC"), track.tpcNClsFound());
    }
  }

  Preslice<aod::McParticles> perMcCollisionTruth = aod::mcparticle::mcCollisionId;
  void processRecoQA(aod::McCollisions const& mcCollisions, aod::McParticles const& mcparticles, MyTracksMC const& tracks, aod::V0Legs const& v0legs, aod::V0PhotonsKF const& v0photonsKF)
  {
    std::unordered_map<int, int> nTracksOfMcId;
    std::unordered_map<int, int> mcIdToBestTrack;
    std::unordered_map<int, float> mcIdBestCrossedRows;
    nTracksOfMcId.reserve(tracks.size());
    for (const auto& track : tracks) {
      if (!track.has_mcParticle()) {
        continue;
      }
      const int mcId = track.mcParticleId();
      nTracksOfMcId[mcId]++;
      const float ncr = track.tpcNClsCrossedRows();
      auto it = mcIdBestCrossedRows.find(mcId);
      if (it == mcIdBestCrossedRows.end() || ncr > it->second) {
        mcIdBestCrossedRows[mcId] = ncr;
        mcIdToBestTrack[mcId] = track.globalIndex();
      }
    }

    std::unordered_map<int, int> nV0LegsOfMcId;
    nV0LegsOfMcId.reserve(v0legs.size());
    for (const auto& leg : v0legs) {
      const int trackId = leg.trackId();
      if (trackId < 0 || trackId >= static_cast<int>(tracks.size())) {
        continue;
      }
      const auto track = tracks.iteratorAt(trackId);
      if (track.has_mcParticle()) {
        nV0LegsOfMcId[track.mcParticleId()]++;
      }
    }
    auto legMcMotherId = [&](auto const& leg) -> int {
      const int trackId = leg.trackId();
      if (trackId < 0 || trackId >= static_cast<int>(tracks.size())) {
        return -1;
      }
      const auto track = tracks.iteratorAt(trackId);
      if (!track.has_mcParticle()) {
        return -1;
      }
      const auto mc = track.mcParticle();
      if (!mc.has_mothers()) {
        return -1;
      }
      return mc.mothersIds()[0];
    };
    std::unordered_map<int, int> nV0CandsOfMcPhoton;
    for (const auto& v0 : v0photonsKF) {
      const int legPosIdx = v0.posTrackId();
      const int legNegIdx = v0.negTrackId();
      if (legPosIdx < 0 || legPosIdx >= static_cast<int>(v0legs.size()) || legNegIdx < 0 || legNegIdx >= static_cast<int>(v0legs.size())) {
        continue;
      }
      const int motherPos = legMcMotherId(v0legs.iteratorAt(legPosIdx));
      const int motherNeg = legMcMotherId(v0legs.iteratorAt(legNegIdx));
      if (motherPos >= 0 && motherPos == motherNeg) {
        nV0CandsOfMcPhoton[motherPos]++;
      }
    }
    for (const auto& mcCollision : mcCollisions) {
      const auto mcparticles_coll = mcparticles.sliceBy(perMcCollisionTruth, mcCollision.globalIndex());
      for (const auto& mctrack : mcparticles_coll) {
        if (mctrack.pdgCode() != PDG_t::kGamma || !mctrack.isPhysicalPrimary()) {
          continue;
        }
        if (std::fabs(mctrack.eta()) > recoQASettingsGroup.cfgPhotonEtaMax || mctrack.pt() < recoQASettingsGroup.cfgPhotonPtMin) {
          continue;
        }
        if (!mctrack.has_daughters()) {
          continue;
        }
        int idPos = -1, idNeg = -1;
        for (const auto& daughterId : mctrack.daughtersIds()) {
          if (daughterId < 0 || daughterId >= static_cast<int>(mcparticles.size())) {
            continue;
          }
          const auto daughter = mcparticles.iteratorAt(daughterId);
          if (daughter.pdgCode() == -PDG_t::kElectron) {
            idPos = daughterId;
          } else if (daughter.pdgCode() == PDG_t::kElectron) {
            idNeg = daughterId;
          }
        }
        if (idPos < 0 || idNeg < 0) {
          continue; // not an e+e- conversion
        }
        const auto daughterPos = mcparticles.iteratorAt(idPos);
        const auto daughterNeg = mcparticles.iteratorAt(idNeg);
        const float rConv = std::sqrt(std::pow(daughterPos.vx(), 2) + std::pow(daughterPos.vy(), 2));
        if (rConv < recoQASettingsGroup.cfgRMinGen || rConv > recoQASettingsGroup.cfgRMaxGen) {
          continue;
        }
        if (std::fabs(daughterPos.eta()) > recoQASettingsGroup.cfgLegEtaMax || std::fabs(daughterNeg.eta()) > recoQASettingsGroup.cfgLegEtaMax) {
          continue;
        }

        fRegistry.fill(HIST("RecoQA/0_converted/hs"), mctrack.pt(), mctrack.eta(), mctrack.phi(), rConv);
        fRegistry.fill(HIST("RecoQA/hStageVsPt"), mctrack.pt(), 0.f);

        auto countIn = [](std::unordered_map<int, int> const& m, int key) -> int {
          const auto it = m.find(key);
          return (it == m.end()) ? 0 : it->second;
        };
        const int nTrkPos = countIn(nTracksOfMcId, idPos);
        const int nTrkNeg = countIn(nTracksOfMcId, idNeg);
        const int nLegPos = countIn(nV0LegsOfMcId, idPos);
        const int nLegNeg = countIn(nV0LegsOfMcId, idNeg);
        const bool posTracked = nTrkPos > 0;
        const bool negTracked = nTrkNeg > 0;
        const bool posInV0Legs = nLegPos > 0;
        const bool negInV0Legs = nLegNeg > 0;

        fRegistry.fill(HIST("RecoQA/Duplicates/hNTracksPerTrueLeg"), nTrkPos);
        fRegistry.fill(HIST("RecoQA/Duplicates/hNTracksPerTrueLeg"), nTrkNeg);
        fRegistry.fill(HIST("RecoQA/Duplicates/hNV0LegsPerTrueLeg"), nLegPos);
        fRegistry.fill(HIST("RecoQA/Duplicates/hNV0LegsPerTrueLeg"), nLegNeg);
        fRegistry.fill(HIST("RecoQA/Duplicates/hNV0CandsPerPhoton_Rconv"), countIn(nV0CandsOfMcPhoton, static_cast<int>(mctrack.globalIndex())), rConv);

        if (posTracked && negTracked) {
          fRegistry.fill(HIST("RecoQA/1_bothLegsTracked/hs"), mctrack.pt(), mctrack.eta(), mctrack.phi(), rConv);
          fRegistry.fill(HIST("RecoQA/hStageVsPt"), mctrack.pt(), 1.f);
        }
        if (posInV0Legs && negInV0Legs) {
          fRegistry.fill(HIST("RecoQA/2_bothLegsInV0Legs/hs"), mctrack.pt(), mctrack.eta(), mctrack.phi(), rConv);
          fRegistry.fill(HIST("RecoQA/hStageVsPt"), mctrack.pt(), 2.f);
        }
        if (countIn(nV0CandsOfMcPhoton, static_cast<int>(mctrack.globalIndex())) > 0) {
          fRegistry.fill(HIST("RecoQA/3_photonInV0PhotonsKF/hs"), mctrack.pt(), mctrack.eta(), mctrack.phi(), rConv);
          fRegistry.fill(HIST("RecoQA/hStageVsPt"), mctrack.pt(), 3.f);
        }
        if (posTracked && !posInV0Legs) {
          fRegistry.fill(HIST("RecoQA/hMissingLeg_PtRconv"), daughterPos.pt(), rConv);
        }
        if (negTracked && !negInV0Legs) {
          fRegistry.fill(HIST("RecoQA/hMissingLeg_PtRconv"), daughterNeg.pt(), rConv);
        }

        if (recoQASettingsGroup.cfgDoDetailedTrackQA) {
          for (const auto& legMcId : {idPos, idNeg}) {
            if (countIn(nTracksOfMcId, legMcId) == 0) {
              continue;
            }
            const auto itTrack = mcIdToBestTrack.find(legMcId);
            if (itTrack == mcIdToBestTrack.end() || itTrack->second < 0 || itTrack->second >= static_cast<int>(tracks.size())) {
              continue;
            }
            const auto track = tracks.iteratorAt(itTrack->second);
            if (countIn(nV0LegsOfMcId, legMcId) > 0) {
              fillLegQualityRecoQA<true>(track);
            } else {
              fillLegQualityRecoQA<false>(track);
            }
          }
        }
      } // end of mc particle loop per true event
    } // end of true event loop
  } // end of reco QA process

  PROCESS_SWITCH(PCMQC, processQC, "run PCM QC", true);
  PROCESS_SWITCH(PCMQC, processQCML, "run PCM QC with ML", false);
  PROCESS_SWITCH(PCMQC, processPCMQCMC, "run PCM QC in MC", false);
  PROCESS_SWITCH(PCMQC, processPCMQCMCML, "run PCM QC in MC with ML cuts", false);
  PROCESS_SWITCH(PCMQC, processGen, "run generated information", false);
  PROCESS_SWITCH(PCMQC, processRecoQA, "run single-photon reconstruction QA on AO2Ds", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& context)
{
  return WorkflowSpec{
    adaptAnalysisTask<PCMQC>(context, TaskName{"pcm-qc"})};
}
