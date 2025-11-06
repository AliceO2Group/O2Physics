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
// This code will create data table for inputs to machine learning for electrons.
//    Please write to: daiki.sekihata@cern.ch

#include "PWGEM/Dilepton/DataModel/lmeeMLTables.h"
#include "PWGEM/Dilepton/Utils/PairUtilities.h"
#include "PWGEM/PhotonMeson/Utils/PCMUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "EventFiltering/Zorro.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DCAFitter/DCAFitterN.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "Math/Vector4D.h"

#include <array>
#include <random>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels>;
using MyCollision = MyCollisions::iterator;

using MyTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU,
                           aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                           aod::pidTOFFullEl, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;
using MyTrack = MyTracks::iterator;

struct TreeCreatorElectronMLDDA {
  SliceCache cache;
  Produces<o2::aod::EMTracksForMLPID> emprimarytracks; // flat table containing collision + track information
  Produces<o2::aod::EMPIDsEl> empidel;
  Produces<o2::aod::EMPIDsPi> empidpi;
  Produces<o2::aod::EMPIDsKa> empidka;
  Produces<o2::aod::EMPIDsPr> empidpr;

  // Basic checks
  HistogramRegistry registry{
    "registry",
    {
      {"Event/hEventCounter", "hEventCounter", {HistType::kTH1F, {{5, 0.5f, 5.5f}}}},
      {"Event/hNumContrib", "Number of contributors to PV;N_{contrib}^{PV};Entries", {HistType::kTH1F, {{65001, -0.5f, 65000.5f}}}},
      {"V0/hAP", "Armenteros Podolanski", {HistType::kTH2F, {{200, -1.f, +1.f}, {250, 0, 0.25}}}},
      {"V0/hXY_Gamma", "photon conversion point in XY;X (cm);Y (cm)", {HistType::kTH2F, {{400, -100, +100}, {400, -100, +100}}}},
      {"V0/hMassGamma_Rxy", "V0 mass gamma", {HistType::kTH2F, {{200, 0, 100}, {100, 0, 0.1}}}},
      {"V0/hCosPA", "V0 cosine of pointing angle", {HistType::kTH1F, {{100, 0.99, 1.f}}}},
      {"V0/hPCA", "V0 distance between 2 legs", {HistType::kTH1F, {{50, 0.f, 0.5f}}}},
      {"V0/hMassGamma", "V0 mass gamma", {HistType::kTH1F, {{100, 0, 0.1}}}},
      {"V0/hMassK0Short", "V0 mass K0S", {HistType::kTH1F, {{200, 0.4, 0.6}}}},
      {"V0/hMassLambda", "V0 mass Lambda", {HistType::kTH1F, {{100, 1.08, 1.18}}}},
      {"V0/hMassAntiLambda", "V0 mass AntiLambda", {HistType::kTH1F, {{100, 1.08, 1.18}}}},

      {"V0/hTPCdEdx_P_El", "TPC dEdx vs. p;p_{in} (GeV/c);TPC dE/dx", {HistType::kTH2F, {{500, 0, 5}, {200, 0, 200}}}},
      {"V0/hTPCdEdx_P_Pi", "TPC dEdx vs. p;p_{in} (GeV/c);TPC dE/dx", {HistType::kTH2F, {{500, 0, 5}, {200, 0, 200}}}},
      {"V0/hTPCdEdx_P_Ka", "TPC dEdx vs. p;p_{in} (GeV/c);TPC dE/dx", {HistType::kTH2F, {{500, 0, 5}, {200, 0, 200}}}},
      {"V0/hTPCdEdx_P_Pr", "TPC dEdx vs. p;p_{in} (GeV/c);TPC dE/dx", {HistType::kTH2F, {{500, 0, 5}, {200, 0, 200}}}},
      {"V0/hTOFbeta_P_El", "TOF beta vs. p;p_{in} (GeV/c);TOF #beta", {HistType::kTH2F, {{500, 0, 5}, {220, 0, 1.1}}}},
      {"V0/hTOFbeta_P_Pi", "TOF beta vs. p;p_{in} (GeV/c);TOF #beta", {HistType::kTH2F, {{500, 0, 5}, {220, 0, 1.1}}}},
      {"V0/hTOFbeta_P_Ka", "TOF beta vs. p;p_{in} (GeV/c);TOF #beta", {HistType::kTH2F, {{500, 0, 5}, {220, 0, 1.1}}}},
      {"V0/hTOFbeta_P_Pr", "TOF beta vs. p;p_{in} (GeV/c);TOF #beta", {HistType::kTH2F, {{500, 0, 5}, {220, 0, 1.1}}}},

      {"Cascade/hRxy_Xi", "R_{xy} of cascade vs. mass;m_{#Lambda#pi};R_{xy} (cm)", {HistType::kTH2F, {{200, 1.2, 1.4}, {200, 0, 20.f}}}},
      {"Cascade/hRxy_Omega", "R_{xy} of cascade vs. mass;m_{#LambdaK};R_{xy} (cm)", {HistType::kTH2F, {{200, 1.6, 1.8}, {200, 0, 20.f}}}},
      {"Cascade/hCTau_Xi", "c#tau vs. mass;m_{#Lambda#pi};c#tau (cm)", {HistType::kTH2F, {{200, 1.2, 1.4}, {200, 0, 20.f}}}},
      {"Cascade/hCTau_Omega", "c#tau vs. mass;m_{#LambdaK};c#tau (cm)", {HistType::kTH2F, {{200, 1.6, 1.8}, {200, 0, 20.f}}}},
      {"Cascade/hV0CosPA", "V0 cosine of pointing angle", {HistType::kTH1F, {{100, 0.99, 1.f}}}},
      {"Cascade/hV0PCA", "V0 distance between 2 legs", {HistType::kTH1F, {{50, 0.f, 0.5}}}},
      {"Cascade/hCosPA", "cascade cosine of pointing angle", {HistType::kTH1F, {{100, 0.99, 1.f}}}},
      {"Cascade/hPCA", "cascade distance between 2 legs", {HistType::kTH1F, {{50, 0.f, 0.5}}}},
      {"Cascade/hMassLambda", "V0 mass Lambda in cascade", {HistType::kTH1F, {{100, 1.08, 1.18}}}},
      {"Cascade/hMassXi", "cascade mass #Xi", {HistType::kTH1F, {{200, 1.2, 1.4}}}},
      {"Cascade/hMassOmega", "cascade mass #Omega", {HistType::kTH1F, {{200, 1.6, 1.8}}}},
      {"Cascade/hMassPt_Xi", "cascade mass #Xi^{#pm};m_{#Lambda#pi} (GeV/c^{2});p_{T,#Lambda#pi} (GeV/c)", {HistType::kTH2F, {{200, 1.2, 1.4}, {100, 0, 10}}}},
      {"Cascade/hMassPt_Omega", "cascade mass #Omega^{#pm};m_{#LambdaK} (GeV/c^{2});p_{T,#LambdaK} (GeV/c)", {HistType::kTH2F, {{200, 1.6, 1.8}, {100, 0, 10}}}},
    },
  };

  // CCDB options
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> mVtxPath{"mVtxPath", "GLO/Calib/MeanVertex", "Path of the mean vertex file"};
  Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};
  // Operation and minimisation criteria
  Configurable<double> d_bz_input{"d_bz_input", -999, "bz field, -999 is automatic"};
  Configurable<int> useMatCorrType{"useMatCorrType", 2, "0: none, 1: TGeo, 2: LUT"};
  Configurable<std::string> irSource{"irSource", "ZNC hadronic", "Estimator of the interaction rate (Recommended: pp --> T0VTX, Pb-Pb --> ZNC hadronic)"};

  // for zorro
  Configurable<std::string> cfg_swt_names{"cfg_swt_names", "fHighTrackMult,fHighFt0cFv0Mult", "comma-separated software trigger names"};
  o2::framework::Configurable<std::string> ccdbPathSoftwareTrigger{"ccdbPathSoftwareTrigger", "EventFiltering/Zorro/", "ccdb path for ZORRO objects"};
  Configurable<uint64_t> bcMarginForSoftwareTrigger{"bcMarginForSoftwareTrigger", 100, "Number of BCs of margin for software triggers"};
  Configurable<bool> cfgUseZorro{"cfgUseZorro", false, "flag to analyze software-triggered data"};

  Configurable<float> downscaling_electron_highP{"downscaling_electron_highP", 1.1, "down scaling factor to store electron at high p"};
  Configurable<float> downscaling_pion_highP{"downscaling_pion_highP", 1.1, "down scaling factor to store pion at high p"};
  Configurable<float> downscaling_kaon_highP{"downscaling_kaon_highP", 1.1, "down scaling factor to store kaon at high p"};
  Configurable<float> downscaling_proton_highP{"downscaling_proton_highP", 1.1, "down scaling factor to store proton at high p"};

  Configurable<float> downscaling_electron_midP{"downscaling_electron_midP", 0.1, "down scaling factor to store electron at intermediate p"};

  Configurable<float> downscaling_electron_lowP{"downscaling_electron_lowP", 0.01, "down scaling factor to store electron at low p"};
  Configurable<float> downscaling_pion_lowP{"downscaling_pion_lowP", 0.01, "down scaling factor to store pion at low p"};
  Configurable<float> downscaling_kaon_lowP{"downscaling_kaon_lowP", 1.1, "down scaling factor to store kaon at low p"};
  Configurable<float> downscaling_proton_lowP{"downscaling_proton_lowP", 0.01, "down scaling factor to store proton at low p"};

  Configurable<float> mid_p_for_downscaling_electron{"mid_p_for_downscaling_electron", 0.8, "intermediate p to apply down scaling factor to store electron"};

  Configurable<float> max_p_for_downscaling_electron{"max_p_for_downscaling_electron", 2.0, "max p to apply down scaling factor to store electron"};
  Configurable<float> max_p_for_downscaling_pion{"max_p_for_downscaling_pion", 2.0, "max p to apply down scaling factor to store pion"};
  Configurable<float> max_p_for_downscaling_kaon{"max_p_for_downscaling_kaon", 0.0, "max p to apply down scaling factor to store kaon"};
  Configurable<float> max_p_for_downscaling_proton{"max_p_for_downscaling_proton", 2.0, "max p to apply down scaling factor to store proton"};
  Configurable<bool> store_ele_band_only{"store_ele_band_only", false, "flag to store tracks around electron band only to reduce output size"};

  struct : ConfigurableGroup {
    std::string prefix = "eventcut_group";
    Configurable<float> cfgZvtxMin{"cfgZvtxMin", -10.f, "min. Zvtx"};
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", 10.f, "max. Zvtx"};
    Configurable<bool> cfgRequireSel8{"cfgRequireSel8", true, "require sel8 in event cut"};
    Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "require FT0AND in event cut"};
    Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", true, "require No time frame border in event cut"};
    Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", true, "require no ITS readout frame border in event cut"};
    Configurable<bool> cfgRequireVertexITSTPC{"cfgRequireVertexITSTPC", false, "require Vertex ITSTPC in event cut"};             // ITS-TPC matched track contributes PV.
    Configurable<bool> cfgRequireVertexTOFmatched{"cfgRequireVertexTOFmatched", false, "require Vertex TOFmatched in event cut"}; // ITS-TPC-TOF matched track contributes PV.
    Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "require no same bunch pileup in event cut"};
    Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", false, "require good Zvtx between FT0 vs. PV in event cut"};
    Configurable<int> cfgTrackOccupancyMin{"cfgTrackOccupancyMin", -2, "min. track occupancy"};
    Configurable<int> cfgTrackOccupancyMax{"cfgTrackOccupancyMax", 1000000000, "max. track occupancy"};
    Configurable<float> cfgFT0COccupancyMin{"cfgFT0COccupancyMin", -2, "min. FT0C occupancy"};
    Configurable<float> cfgFT0COccupancyMax{"cfgFT0COccupancyMax", 1000000000, "max. FT0C occupancy"};
    Configurable<bool> cfgRequireNoCollInTimeRangeStandard{"cfgRequireNoCollInTimeRangeStandard", false, "require no collision in time range standard"};
    Configurable<bool> cfgRequireNoCollInTimeRangeStrict{"cfgRequireNoCollInTimeRangeStrict", false, "require no collision in time range strict"};
    Configurable<bool> cfgRequirekNoCollInRofStandard{"cfgRequirekNoCollInRofStandard", false, "require no other collisions in this Readout Frame with per-collision multiplicity above threshold"};
    Configurable<bool> cfgRequirekNoCollInRofStrict{"cfgRequirekNoCollInRofStrict", false, "require no other collisions in this Readout Frame"};
    Configurable<bool> cfgRequirekNoHighMultCollInPrevRof{"cfgRequirekNoHighMultCollInPrevRof", false, "require no HM collision in previous ITS ROF"};
    Configurable<bool> cfgRequireGoodITSLayer3{"cfgRequireGoodITSLayer3", false, "number of inactive chips on ITS layer 3 are below threshold "};
    Configurable<bool> cfgRequireGoodITSLayer0123{"cfgRequireGoodITSLayer0123", false, "number of inactive chips on ITS layers 0-3 are below threshold "};
    Configurable<bool> cfgRequireGoodITSLayersAll{"cfgRequireGoodITSLayersAll", false, "number of inactive chips on all ITS layers are below threshold "};
    Configurable<uint16_t> cfgNumContribMin{"cfgNumContribMin", 0, "min. numContrib"};
    Configurable<uint16_t> cfgNumContribMax{"cfgNumContribMax", 65000, "max. numContrib"};
  } eventcuts;

  struct : ConfigurableGroup {
    std::string prefix = "v0cut_group";
    Configurable<float> cfg_min_pt{"cfg_min_pt", 0.05, "min pt for v0 legs"};
    Configurable<float> cfg_max_eta{"cfg_max_eta", 0.9, "max. eta for v0 legs"};
    Configurable<float> cfg_min_mass_photon{"cfg_min_mass_photon", 0.00, "min mass for photon conversion"};
    Configurable<float> cfg_max_mass_photon{"cfg_max_mass_photon", 0.02, "max mass for photon conversion"};
    Configurable<float> cfg_min_mass_k0s{"cfg_min_mass_k0s", 0.490, "min mass for K0S"};
    Configurable<float> cfg_max_mass_k0s{"cfg_max_mass_k0s", 0.505, "max mass for K0S"};
    Configurable<float> cfg_min_mass_lambda{"cfg_min_mass_lambda", 1.113, "min mass for Lambda"};
    Configurable<float> cfg_max_mass_lambda{"cfg_max_mass_lambda", 1.118, "max mass for Lambda"};

    Configurable<float> cfg_min_mass_k0s_veto{"cfg_min_mass_k0s_veto", 0.47, "min mass for K0S veto"};
    Configurable<float> cfg_max_mass_k0s_veto{"cfg_max_mass_k0s_veto", 0.52, "max mass for K0S veto"};
    Configurable<float> cfg_min_mass_lambda_veto{"cfg_min_mass_lambda_veto", 1.105, "min mass for Lambda veto"};
    Configurable<float> cfg_max_mass_lambda_veto{"cfg_max_mass_lambda_veto", 1.125, "max mass for Lambda veto"};

    Configurable<float> cfg_min_cospa{"cfg_min_cospa", 0.9998, "min cospa for v0"};
    Configurable<float> cfg_max_dcadau{"cfg_max_dcadau", 0.1, "max distance between 2 legs for v0"};
    Configurable<float> cfg_min_qt_strangeness{"cfg_min_qt_strangeness", 0.02, "min qt for Lambda and K0S"};
    Configurable<float> cfg_min_qt_k0s{"cfg_min_qt_k0s", 0.11, "min qt for K0S"};

    Configurable<float> cfg_min_cr2findable_ratio_tpc{"cfg_min_cr2findable_ratio_tpc", 0.8, "min. TPC Ncr/Nf ratio"};
    Configurable<float> cfg_max_frac_shared_clusters_tpc{"cfg_max_frac_shared_clusters_tpc", 0.7, "max fraction of shared clusters in TPC"};
    Configurable<int> cfg_min_ncrossedrows_tpc{"cfg_min_ncrossedrows_tpc", 70, "min ncrossed rows"};
    Configurable<int> cfg_min_ncluster_tpc_pid{"cfg_min_ncluster_tpc_pid", 0, "min ncluster tpc"};
    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<int> cfg_min_ncluster_its{"cfg_min_ncluster_its", 2, "min ncluster its"};
    Configurable<int> cfg_min_ncluster_itsib{"cfg_min_ncluster_itsib", 0, "min ncluster itsib"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 5.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 36.0, "max chi2/NclsITS"};
    Configurable<float> cfg_min_chi2its{"cfg_min_chi2its", 0.0, "min chi2/NclsITS"}; // remove ITS afterburner
    Configurable<float> cfg_min_dcaxy_v0leg{"cfg_min_dcaxy_v0leg", 0.1, "min dca XY to PV for v0 legs in cm"};
    Configurable<bool> cfg_includeITSsa{"cfg_includeITSsa", false, "Flag to include ITSsa tracks"};
    Configurable<float> cfg_max_pt_itssa{"cfg_max_pt_itssa", 0.15, "mix pt for ITSsa track"};

    Configurable<float> cfg_min_TPCNsigmaEl{"cfg_min_TPCNsigmaEl", -5, "min n sigma e in TPC"};
    Configurable<float> cfg_max_TPCNsigmaEl{"cfg_max_TPCNsigmaEl", +5, "max n sigma e in TPC"};
    Configurable<float> cfg_min_TPCNsigmaPi{"cfg_min_TPCNsigmaPi", -5, "min n sigma pi in TPC"};
    Configurable<float> cfg_max_TPCNsigmaPi{"cfg_max_TPCNsigmaPi", +5, "max n sigma pi in TPC"};
    Configurable<float> cfg_min_TPCNsigmaKa{"cfg_min_TPCNsigmaKa", -5, "min n sigma ka in TPC"};
    Configurable<float> cfg_max_TPCNsigmaKa{"cfg_max_TPCNsigmaKa", +5, "max n sigma ka in TPC"};
    Configurable<float> cfg_min_TPCNsigmaPr{"cfg_min_TPCNsigmaPr", -5, "min n sigma pr in TPC"};
    Configurable<float> cfg_max_TPCNsigmaPr{"cfg_max_TPCNsigmaPr", +5, "max n sigma pr in TPC"};

    Configurable<float> cfg_min_TOFNsigmaEl{"cfg_min_TOFNsigmaEl", -5, "min n sigma e in TOF"};
    Configurable<float> cfg_max_TOFNsigmaEl{"cfg_max_TOFNsigmaEl", +5, "max n sigma e in TOF"};
    Configurable<float> cfg_min_TOFNsigmaPi{"cfg_min_TOFNsigmaPi", -5, "min n sigma pi in TOF"};
    Configurable<float> cfg_max_TOFNsigmaPi{"cfg_max_TOFNsigmaPi", +5, "max n sigma pi in TOF"};
    Configurable<float> cfg_min_TOFNsigmaKa{"cfg_min_TOFNsigmaKa", -5, "min n sigma ka in TOF"};
    Configurable<float> cfg_max_TOFNsigmaKa{"cfg_max_TOFNsigmaKa", +5, "max n sigma ka in TOF"};
    Configurable<float> cfg_min_TOFNsigmaPr{"cfg_min_TOFNsigmaPr", -5, "min n sigma pr in TOF"};
    Configurable<float> cfg_max_TOFNsigmaPr{"cfg_max_TOFNsigmaPr", +5, "max n sigma pr in TOF"};
  } v0cuts;

  struct : ConfigurableGroup {
    std::string prefix = "tightv0cut_group";
    Configurable<float> cfg_min_cr2findable_ratio_tpc{"cfg_min_cr2findable_ratio_tpc", 0.8, "min. TPC Ncr/Nf ratio"};
    Configurable<float> cfg_max_frac_shared_clusters_tpc{"cfg_max_frac_shared_clusters_tpc", 0.7, "max fraction of shared clusters in TPC"};
    Configurable<int> cfg_min_ncrossedrows_tpc{"cfg_min_ncrossedrows_tpc", 80, "min ncrossed rows"};
    Configurable<int> cfg_min_ncluster_tpc_pid{"cfg_min_ncluster_tpc_pid", 60, "min ncluster tpc"};
    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<int> cfg_min_ncluster_its{"cfg_min_ncluster_its", 2, "min ncluster its"};
    Configurable<int> cfg_min_ncluster_itsib{"cfg_min_ncluster_itsib", 0, "min ncluster itsib"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 5.0, "max chi2/NclsITS"};
    Configurable<float> cfg_min_chi2its{"cfg_min_chi2its", -1e+10, "min chi2/NclsITS"}; // remove ITS afterburner
    Configurable<float> cfg_max_chi2tof{"cfg_max_chi2tof", 1e+10, "max chi2 TOF"};      // distance in cm

    Configurable<float> cfg_min_TPCNsigmaEl{"cfg_min_TPCNsigmaEl", -2, "min n sigma e in TPC for pc->ee"};
    Configurable<float> cfg_max_TPCNsigmaEl{"cfg_max_TPCNsigmaEl", +2, "max n sigma e in TPC for pc->ee"};
    Configurable<float> cfg_min_TOFNsigmaEl{"cfg_min_TOFNsigmaEl", -2, "min n sigma e in TOF for pc->ee"};
    Configurable<float> cfg_max_TOFNsigmaEl{"cfg_max_TOFNsigmaEl", +2, "max n sigma e in TOF for pc->ee"};

    Configurable<float> cfg_min_TPCNsigmaPi{"cfg_min_TPCNsigmaPi", -2, "min n sigma pi in TPC for Lambda and cascade"};
    Configurable<float> cfg_max_TPCNsigmaPi{"cfg_max_TPCNsigmaPi", +2, "max n sigma pi in TPC for Lambda and cascade"};
    Configurable<float> cfg_min_TPCNsigmaPr{"cfg_min_TPCNsigmaPr", -2, "min n sigma pr in TPC for cascade"};
    Configurable<float> cfg_max_TPCNsigmaPr{"cfg_max_TPCNsigmaPr", +2, "max n sigma pr in TPC for cascade"};

    Configurable<float> cfg_min_TOFNsigmaPi{"cfg_min_TOFNsigmaPi", -2, "min n sigma pi in TOF for Lambda and cascade"};
    Configurable<float> cfg_max_TOFNsigmaPi{"cfg_max_TOFNsigmaPi", +2, "max n sigma pi in TOF for Lambda and cascade"};
    Configurable<float> cfg_min_TOFNsigmaPr{"cfg_min_TOFNsigmaPr", -2, "min n sigma pr in TOF for cascade"};
    Configurable<float> cfg_max_TOFNsigmaPr{"cfg_max_TOFNsigmaPr", +2, "max n sigma pr in TOF for cascade"};
  } tightv0cuts;

  struct : ConfigurableGroup {
    std::string prefix = "cascadecut_group";
    Configurable<float> cfg_min_mass_lambda{"cfg_min_mass_lambda", 1.11, "min mass for lambda in cascade"};
    Configurable<float> cfg_max_mass_lambda{"cfg_max_mass_lambda", 1.12, "max mass for lambda in cascade"};
    Configurable<float> cfg_min_mass_Xi_veto{"cfg_min_mass_Xi_veto", 1.31, "min mass for Xi veto"};
    Configurable<float> cfg_max_mass_Xi_veto{"cfg_max_mass_Xi_veto", 1.33, "max mass for Xi veto"};
    Configurable<float> cfg_min_mass_Omega{"cfg_min_mass_Omega", 1.669, "min mass for Omega"};
    Configurable<float> cfg_max_mass_Omega{"cfg_max_mass_Omega", 1.675, "max mass for Omega"};
    Configurable<float> cfg_min_cospa_v0{"cfg_min_cospa_v0", 0.995, "minimum V0 CosPA in cascade"};
    Configurable<float> cfg_max_dcadau_v0{"cfg_max_dcadau_v0", 0.1, "max distance between V0 Daughters in cascade"};
    Configurable<float> cfg_min_cospa{"cfg_min_cospa", 0.9998, "minimum cascade CosPA"};
    Configurable<float> cfg_max_dcadau{"cfg_max_dcadau", 0.1, "max distance between bachelor and V0"};
    Configurable<float> cfg_min_rxy_v0{"cfg_min_rxy_v0", 1.2, "minimum V0 rxy in cascade"};
    Configurable<float> cfg_min_rxy{"cfg_min_rxy", 0.5, "minimum V0 rxy in cascade"};
    Configurable<float> cfg_min_dcaxy_v0leg{"cfg_min_dcaxy_v0leg", 0.1, "min dca XY for v0 legs in cm"};
    Configurable<float> cfg_min_dcaxy_bachelor{"cfg_min_dcaxy_bachelor", 0.05, "min dca XY for bachelor in cm"};
    Configurable<float> cfg_min_dcaxy_v0{"cfg_min_dcaxy_v0", 0.05, "min dca XY for V0 in cm"};
  } cascadecuts;

  // for RCT
  Configurable<bool> cfgRequireGoodRCT{"cfgRequireGoodRCT", false, "require good detector flag in run condtion table"};
  Configurable<std::string> cfgRCTLabel{"cfgRCTLabel", "CBT_hadronPID", "select 1 [CBT, CBT_hadronPID, CBT_muon_glo] see O2Physics/Common/CCDB/RCTSelectionFlags.h"};
  Configurable<bool> cfgCheckZDC{"cfgCheckZDC", false, "set ZDC flag for PbPb"};
  Configurable<bool> cfgTreatLimitedAcceptanceAsBad{"cfgTreatLimitedAcceptanceAsBad", false, "reject all events where the detectors relevant for the specified Runlist are flagged as LimitedAcceptance"};

  int mRunNumber;
  float d_bz;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  o2::dataformats::VertexBase mVtx;
  const o2::dataformats::MeanVertexObject* mMeanVtx = nullptr;
  o2::base::MatLayerCylSet* lut = nullptr;
  o2::dataformats::DCA mDcaInfoCov;
  o2::aod::rctsel::RCTFlagsChecker rctChecker;
  ctpRateFetcher mRateFetcher;
  Zorro zorro;

  std::mt19937 engine;
  std::uniform_real_distribution<float> dist01;

  void init(InitContext&)
  {
    mRunNumber = 0;
    d_bz = 0;

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    rctChecker.init(cfgRCTLabel.value, cfgCheckZDC.value, cfgTreatLimitedAcceptanceAsBad.value);

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

    std::random_device seed_gen;
    engine = std::mt19937(seed_gen());
    dist01 = std::uniform_real_distribution<float>(0.0f, 1.0f);
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    if (cfgUseZorro) {
      zorro.setCCDBpath(ccdbPathSoftwareTrigger);
      zorro.setBCtolerance(bcMarginForSoftwareTrigger); // this does nothing.
      zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), cfg_swt_names.value);
      zorro.populateHistRegistry(registry, bc.runNumber());
    }

    // load matLUT for this timestamp
    if (!lut) {
      LOG(info) << "Loading material look-up table for timestamp: " << bc.timestamp();
      lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->getForTimeStamp<o2::base::MatLayerCylSet>(lutPath, bc.timestamp()));
    } else {
      LOG(info) << "Material look-up table already in place. Not reloading.";
    }

    // In case override, don't proceed, please - no CCDB access required
    if (d_bz_input > -990) {
      d_bz = d_bz_input;
      o2::parameters::GRPMagField grpmag;
      if (std::fabs(d_bz) > 1e-5) {
        grpmag.setL3Current(30000.f / (d_bz / 5.0f));
      }
      o2::base::Propagator::initFieldFromGRP(&grpmag);
      o2::base::Propagator::Instance()->setMatLUT(lut);
      mMeanVtx = ccdb->getForTimeStamp<o2::dataformats::MeanVertexObject>(mVtxPath, bc.timestamp());
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
      o2::base::Propagator::Instance()->setMatLUT(lut);
      mMeanVtx = ccdb->getForTimeStamp<o2::dataformats::MeanVertexObject>(mVtxPath, bc.timestamp());
      // Fetch magnetic field from ccdb for current collision
      d_bz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      o2::base::Propagator::Instance()->setMatLUT(lut);
      mMeanVtx = ccdb->getForTimeStamp<o2::dataformats::MeanVertexObject>(mVtxPath, bc.timestamp());
      // Fetch magnetic field from ccdb for current collision
      d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    }
    mRunNumber = bc.runNumber();
  }

  template <int il0, int il1, typename TTrack>
  float meanClusterSizeITS(TTrack const& track)
  {
    int total_size = 0;
    int nl = 0;
    for (int il = il0; il < il1; il++) {
      if (track.itsClsSizeInLayer(il) > 0) {
        total_size += track.itsClsSizeInLayer(il);
        nl++;
      }
    }
    if (nl > 0) {
      return static_cast<float>(total_size) / static_cast<float>(nl);
    } else {
      return 0.f;
    }
  }

  template <typename TCollision, typename TTrack>
  bool isSelectedV0Leg(TCollision const& collision, TTrack const& track)
  {
    if (!track.hasITS()) {
      return false;
    }

    if (track.itsNCls() < v0cuts.cfg_min_ncluster_its) {
      return false;
    }
    if (track.itsNClsInnerBarrel() < v0cuts.cfg_min_ncluster_itsib) {
      return false;
    }
    if (track.itsChi2NCl() < v0cuts.cfg_min_chi2its || v0cuts.cfg_max_chi2its < track.itsChi2NCl()) {
      return false;
    }

    if (!v0cuts.cfg_includeITSsa && (!track.hasITS() || !track.hasTPC())) {
      return false;
    }

    if (track.hasTPC()) {
      if (track.tpcNClsCrossedRows() < v0cuts.cfg_min_ncrossedrows_tpc) {
        return false;
      }
      if (track.tpcNClsFound() < v0cuts.cfg_min_ncluster_tpc) {
        return false;
      }
      if (track.tpcChi2NCl() > v0cuts.cfg_max_chi2tpc) {
        return false;
      }
      if (track.tpcCrossedRowsOverFindableCls() < v0cuts.cfg_min_cr2findable_ratio_tpc) {
        return false;
      }
      if (track.tpcFractionSharedCls() > v0cuts.cfg_max_frac_shared_clusters_tpc) {
        return false;
      }
      if (track.tpcNClsPID() < v0cuts.cfg_min_ncluster_tpc_pid) {
        return false;
      }
    }

    mDcaInfoCov.set(999, 999, 999, 999, 999);
    auto trackParCov = getTrackParCov(track);
    trackParCov.setPID(track.pidForTracking());
    mVtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
    mVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
    o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, trackParCov, 2.f, matCorr, &mDcaInfoCov);
    // float dcaXY = mDcaInfoCov.getY();
    // float dcaZ = mDcaInfoCov.getZ();

    // if (std::fabs(dcaXY) < v0cuts.cfg_min_dcaxy_v0leg) { // this is applied in filter.
    //   return false;
    // }

    if (std::fabs(trackParCov.getEta()) > v0cuts.cfg_max_eta || trackParCov.getPt() < v0cuts.cfg_min_pt) {
      return false;
    }

    if ((track.hasITS() && !track.hasTPC() && !track.hasTRD() && !track.hasTOF()) && v0cuts.cfg_max_pt_itssa < track.pt()) {
      return true;
    }

    return true;
  }

  template <typename TCollision, typename TTrack>
  bool isSelectedV0LegTight(TCollision const& collision, TTrack const& track)
  {
    if (!track.hasITS() || !track.hasTPC()) {
      return false;
    }

    if (track.itsNCls() < tightv0cuts.cfg_min_ncluster_its) {
      return false;
    }
    if (track.itsNClsInnerBarrel() < tightv0cuts.cfg_min_ncluster_itsib) {
      return false;
    }
    if (track.itsChi2NCl() < tightv0cuts.cfg_min_chi2its || tightv0cuts.cfg_max_chi2its < track.itsChi2NCl()) {
      return false;
    }

    if (track.tpcNClsCrossedRows() < tightv0cuts.cfg_min_ncrossedrows_tpc) {
      return false;
    }
    if (track.tpcNClsFound() < tightv0cuts.cfg_min_ncluster_tpc) {
      return false;
    }
    if (track.tpcChi2NCl() < 0.f || tightv0cuts.cfg_max_chi2tpc < track.tpcChi2NCl()) {
      return false;
    }
    if (track.tpcCrossedRowsOverFindableCls() < tightv0cuts.cfg_min_cr2findable_ratio_tpc) {
      return false;
    }
    if (track.tpcFractionSharedCls() > tightv0cuts.cfg_max_frac_shared_clusters_tpc) {
      return false;
    }
    if (track.tpcNClsPID() < tightv0cuts.cfg_min_ncluster_tpc_pid) {
      return false;
    }

    mDcaInfoCov.set(999, 999, 999, 999, 999);
    auto trackParCov = getTrackParCov(track);
    trackParCov.setPID(track.pidForTracking());
    mVtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
    mVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
    o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, trackParCov, 2.f, matCorr, &mDcaInfoCov);
    // float dcaXY = mDcaInfoCov.getY();
    // float dcaZ = mDcaInfoCov.getZ();

    // if (std::fabs(dcaXY) < v0cuts.cfg_min_dcaxy_v0leg) { // this is applied in filter.
    //   return false;
    // }

    if (std::fabs(trackParCov.getEta()) > v0cuts.cfg_max_eta || trackParCov.getPt() < v0cuts.cfg_min_pt) {
      return false;
    }

    return true;
  }

  template <typename TTrack>
  bool isElectron(TTrack const& track)
  {
    if (v0cuts.cfg_includeITSsa && (track.hasITS() && !track.hasTPC() && !track.hasTRD() && !track.hasTOF())) {
      return true;
    }
    bool is_El_TPC = v0cuts.cfg_min_TPCNsigmaEl < track.tpcNSigmaEl() && track.tpcNSigmaEl() < v0cuts.cfg_max_TPCNsigmaEl;
    bool is_El_TOF = track.hasTOF() ? v0cuts.cfg_min_TOFNsigmaEl < track.tofNSigmaEl() && track.tofNSigmaEl() < v0cuts.cfg_max_TOFNsigmaEl : true; // TOFif
    return is_El_TPC && is_El_TOF;
  }

  template <typename TTrack>
  bool isPion(TTrack const& track)
  {
    if (v0cuts.cfg_includeITSsa && (track.hasITS() && !track.hasTPC() && !track.hasTRD() && !track.hasTOF())) {
      return true;
    }
    bool is_Pi_TPC = v0cuts.cfg_min_TPCNsigmaPi < track.tpcNSigmaPi() && track.tpcNSigmaPi() < v0cuts.cfg_max_TPCNsigmaPi;
    bool is_Pi_TOF = track.hasTOF() ? v0cuts.cfg_min_TOFNsigmaPi < track.tofNSigmaPi() && track.tofNSigmaPi() < v0cuts.cfg_max_TOFNsigmaPi : true; // TOFif
    return is_Pi_TPC && is_Pi_TOF;
  }

  template <typename TTrack>
  bool isKaon(TTrack const& track)
  {
    if (v0cuts.cfg_includeITSsa && (track.hasITS() && !track.hasTPC() && !track.hasTRD() && !track.hasTOF())) {
      return true;
    }
    bool is_Ka_TPC = v0cuts.cfg_min_TPCNsigmaKa < track.tpcNSigmaKa() && track.tpcNSigmaKa() < v0cuts.cfg_max_TPCNsigmaKa;
    bool is_Ka_TOF = track.hasTOF() ? v0cuts.cfg_min_TOFNsigmaKa < track.tofNSigmaKa() && track.tofNSigmaKa() < v0cuts.cfg_max_TOFNsigmaKa : true; // TOFif
    return is_Ka_TPC && is_Ka_TOF;
  }

  template <typename TTrack>
  bool isProton(TTrack const& track)
  {
    if (v0cuts.cfg_includeITSsa && (track.hasITS() && !track.hasTPC() && !track.hasTRD() && !track.hasTOF())) {
      return true;
    }
    bool is_Pr_TPC = v0cuts.cfg_min_TPCNsigmaPr < track.tpcNSigmaPr() && track.tpcNSigmaPr() < v0cuts.cfg_max_TPCNsigmaPr;
    bool is_Pr_TOF = track.hasTOF() ? v0cuts.cfg_min_TOFNsigmaPr < track.tofNSigmaPr() && track.tofNSigmaPr() < v0cuts.cfg_max_TOFNsigmaPr : true; // TOFif
    return is_Pr_TPC && is_Pr_TOF;
  }

  template <typename TTrack>
  bool isElectronTight(TTrack const& track)
  {
    bool is_El_TPC = tightv0cuts.cfg_min_TPCNsigmaEl < track.tpcNSigmaEl() && track.tpcNSigmaEl() < tightv0cuts.cfg_max_TPCNsigmaEl;
    bool is_El_TOF = track.hasTOF() ? tightv0cuts.cfg_min_TOFNsigmaEl < track.tofNSigmaEl() && track.tofNSigmaEl() < tightv0cuts.cfg_max_TOFNsigmaEl : true; // TOFif
    return is_El_TPC && is_El_TOF;
  }

  template <typename TTrack>
  bool isPionTight(TTrack const& track)
  {
    bool is_Pi_TPC = tightv0cuts.cfg_min_TPCNsigmaPi < track.tpcNSigmaPi() && track.tpcNSigmaPi() < tightv0cuts.cfg_max_TPCNsigmaPi;
    bool is_Pi_TOF = track.hasTOF() ? tightv0cuts.cfg_min_TOFNsigmaPi < track.tofNSigmaPi() && track.tofNSigmaPi() < tightv0cuts.cfg_max_TOFNsigmaPi : true; // TOFif
    return is_Pi_TPC && is_Pi_TOF;
  }

  template <typename TTrack>
  bool isProtonTight(TTrack const& track)
  {
    bool is_Pr_TPC = tightv0cuts.cfg_min_TPCNsigmaPr < track.tpcNSigmaPr() && track.tpcNSigmaPr() < tightv0cuts.cfg_max_TPCNsigmaPr;
    bool is_Pr_TOF = track.hasTOF() ? tightv0cuts.cfg_min_TOFNsigmaPr < track.tofNSigmaPr() && track.tofNSigmaPr() < tightv0cuts.cfg_max_TOFNsigmaPr : true; // TOFif
    return is_Pr_TPC && is_Pr_TOF;
  }

  template <typename TTrack>
  bool isPionTightTOFreq(TTrack const& track)
  {
    // only for K0S-> pi+ pi-
    bool is_Pi_TPC = tightv0cuts.cfg_min_TPCNsigmaPi < track.tpcNSigmaPi() && track.tpcNSigmaPi() < tightv0cuts.cfg_max_TPCNsigmaPi;
    bool is_Pi_TOF = tightv0cuts.cfg_min_TOFNsigmaPi < track.tofNSigmaPi() && track.tofNSigmaPi() < tightv0cuts.cfg_max_TOFNsigmaPi && std::fabs(track.tofChi2()) < tightv0cuts.cfg_max_chi2tof; // TOFreq
    return is_Pi_TPC && is_Pi_TOF;
  }

  template <typename TCollision, typename TTrack>
  void fillTrackTable(TCollision const& collision, TTrack const& track, const uint8_t pidlabel, const float hadronicRate)
  {
    if (store_ele_band_only && !isElectron(track)) {
      return;
    }
    mDcaInfoCov.set(999, 999, 999, 999, 999);
    auto trackParCov = getTrackParCov(track);
    // trackParCov.setPID(track.pidForTracking());
    trackParCov.setPID(o2::track::PID::Electron); // This is for eID in the end!
    mVtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
    mVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
    o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, trackParCov, 2.f, matCorr, &mDcaInfoCov);
    // float dcaXY = mDcaInfoCov.getY();
    // float dcaZ = mDcaInfoCov.getZ();

    if (pidlabel == static_cast<uint8_t>(o2::aod::pwgem::dilepton::ml::PID_Label::kElectron)) {
      if (track.tpcInnerParam() < mid_p_for_downscaling_electron) {
        if (dist01(engine) > downscaling_electron_lowP) {
          return;
        }
      } else if (track.tpcInnerParam() < max_p_for_downscaling_electron) {
        if (dist01(engine) > downscaling_electron_midP) {
          return;
        }
      } else {
        if (dist01(engine) > downscaling_electron_highP) {
          return;
        }
      }
    } else if (pidlabel == static_cast<uint8_t>(o2::aod::pwgem::dilepton::ml::PID_Label::kPion)) {
      if (track.tpcInnerParam() < max_p_for_downscaling_pion) {
        if (dist01(engine) > downscaling_pion_lowP) {
          return;
        }
      } else {
        if (dist01(engine) > downscaling_pion_highP) {
          return;
        }
      }
    } else if (pidlabel == static_cast<uint8_t>(o2::aod::pwgem::dilepton::ml::PID_Label::kKaon)) {
      if (track.tpcInnerParam() < max_p_for_downscaling_kaon) {
        if (dist01(engine) > downscaling_kaon_lowP) {
          return;
        }
      } else {
        if (dist01(engine) > downscaling_kaon_highP) {
          return;
        }
      }
    } else if (pidlabel == static_cast<uint8_t>(o2::aod::pwgem::dilepton::ml::PID_Label::kProton)) {
      if (track.tpcInnerParam() < max_p_for_downscaling_proton) {
        if (dist01(engine) > downscaling_proton_lowP) {
          return;
        }
      } else {
        if (dist01(engine) > downscaling_proton_highP) {
          return;
        }
      }
    }

    if (std::find(stored_trackIds.begin(), stored_trackIds.end(), track.globalIndex()) == stored_trackIds.end()) {
      emprimarytracks(collision.numContrib(), collision.trackOccupancyInTimeRange(), collision.ft0cOccupancyInTimeRange(), hadronicRate,
                      trackParCov.getP(), trackParCov.getTgl(), track.sign(),
                      track.tpcNClsFindable(), track.tpcNClsFound(), track.tpcNClsCrossedRows(), track.tpcNClsPID(),
                      track.tpcChi2NCl(), track.tpcInnerParam(),
                      track.tpcSignal(),
                      track.beta(),
                      track.itsClusterSizes(), track.itsChi2NCl(), track.tofChi2(), track.detectorMap(), pidlabel);

      empidel(track.tpcNSigmaEl(), track.tofNSigmaEl());
      empidpi(track.tpcNSigmaPi(), track.tofNSigmaPi());
      empidka(track.tpcNSigmaKa(), track.tofNSigmaKa());
      empidpr(track.tpcNSigmaPr(), track.tofNSigmaPr());

      stored_trackIds.emplace_back(track.globalIndex());
    }
  }

  template <typename TCollision>
  bool isSelectedEvent(TCollision const& collision)
  {
    if (eventcuts.cfgRequireSel8 && !collision.sel8()) {
      return false;
    }

    if (collision.posZ() < eventcuts.cfgZvtxMin || eventcuts.cfgZvtxMax < collision.posZ()) {
      return false;
    }

    if (eventcuts.cfgRequireFT0AND && !collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
      return false;
    }

    if (eventcuts.cfgRequireNoTFB && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }

    if (eventcuts.cfgRequireNoITSROFB && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }

    if (eventcuts.cfgRequireVertexITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return false;
    }

    if (eventcuts.cfgRequireVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
      return false;
    }

    if (eventcuts.cfgRequireNoSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }

    if (eventcuts.cfgRequireGoodZvtxFT0vsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }

    if (eventcuts.cfgRequireNoCollInTimeRangeStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return false;
    }

    if (eventcuts.cfgRequireNoCollInTimeRangeStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict)) {
      return false;
    }

    if (eventcuts.cfgRequirekNoCollInRofStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      return false;
    }

    if (eventcuts.cfgRequirekNoCollInRofStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStrict)) {
      return false;
    }

    if (eventcuts.cfgRequirekNoHighMultCollInPrevRof && !collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
      return false;
    }

    if (eventcuts.cfgRequireGoodITSLayer3 && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayer3)) {
      return false;
    }

    if (eventcuts.cfgRequireGoodITSLayer0123 && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayer0123)) {
      return false;
    }

    if (eventcuts.cfgRequireGoodITSLayersAll && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      return false;
    }

    if (!(eventcuts.cfgTrackOccupancyMin <= collision.trackOccupancyInTimeRange() && collision.trackOccupancyInTimeRange() < eventcuts.cfgTrackOccupancyMax)) {
      return false;
    }

    if (!(eventcuts.cfgFT0COccupancyMin <= collision.ft0cOccupancyInTimeRange() && collision.ft0cOccupancyInTimeRange() < eventcuts.cfgFT0COccupancyMax)) {
      return false;
    }

    return true;
  }

  //! type of V0. 0: built solely for cascades (does not pass standard V0 cuts), 1: standard 2, 3: photon-like with TPC-only use. Regular analysis should always use type 1.
  Filter v0Filter = o2::aod::v0data::v0Type == uint8_t(1) && o2::aod::v0data::v0cosPA > v0cuts.cfg_min_cospa&& o2::aod::v0data::dcaV0daughters<v0cuts.cfg_max_dcadau && nabs(o2::aod::v0data::dcanegtopv)> v0cuts.cfg_min_dcaxy_v0leg&& nabs(o2::aod::v0data::dcanegtopv) > v0cuts.cfg_min_dcaxy_v0leg;
  using filteredV0s = soa::Filtered<aod::V0Datas>;

  Filter cascadeFilter = o2::aod::cascdata::dcacascdaughters < cascadecuts.cfg_max_dcadau && nabs(o2::aod::cascdata::dcanegtopv) > cascadecuts.cfg_min_dcaxy_v0leg&& nabs(o2::aod::cascdata::dcanegtopv) > cascadecuts.cfg_min_dcaxy_v0leg&& nabs(o2::aod::cascdata::dcabachtopv) > cascadecuts.cfg_min_dcaxy_bachelor;
  using filteredCascades = soa::Filtered<aod::CascDatas>;

  Filter collisionFilter_track_occupancy = eventcuts.cfgTrackOccupancyMin <= o2::aod::evsel::trackOccupancyInTimeRange && o2::aod::evsel::trackOccupancyInTimeRange < eventcuts.cfgTrackOccupancyMax;
  Filter collisionFilter_ft0c_occupancy = eventcuts.cfgFT0COccupancyMin <= o2::aod::evsel::ft0cOccupancyInTimeRange && o2::aod::evsel::ft0cOccupancyInTimeRange < eventcuts.cfgFT0COccupancyMax;
  Filter collisionFilter_numContrib = eventcuts.cfgNumContribMin <= o2::aod::collision::numContrib && o2::aod::collision::numContrib < eventcuts.cfgNumContribMax;
  Filter collisionFilter_common = nabs(o2::aod::collision::posZ) < 10.f && o2::aod::evsel::sel8 == true;
  using filteredMyCollisions = soa::Filtered<MyCollisions>;

  Preslice<aod::V0Datas> perCollision_v0 = o2::aod::v0data::collisionId;
  Preslice<aod::CascDatas> perCollision_cascade = o2::aod::cascdata::collisionId;
  Preslice<MyTracks> perCollision_track = o2::aod::track::collisionId;

  Partition<MyTracks> posTracks = o2::aod::track::signed1Pt > 0.f && ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::ITS) == true;
  Partition<MyTracks> negTracks = o2::aod::track::signed1Pt < 0.f && ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::ITS) == true;
  std::vector<uint64_t> stored_trackIds;

  void processPID(filteredMyCollisions const& collisions, aod::BCsWithTimestamps const&, filteredV0s const& v0s, filteredCascades const& cascades, MyTracks const& tracks)
  {
    stored_trackIds.reserve(tracks.size());
    for (const auto& collision : collisions) {
      registry.fill(HIST("Event/hEventCounter"), 1.0); // all

      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!isSelectedEvent(collision)) {
        continue;
      }

      if (cfgUseZorro && !zorro.isSelected(bc.globalBC(), bcMarginForSoftwareTrigger)) {
        continue;
      }

      if (cfgRequireGoodRCT && !rctChecker.checkTable(collision)) {
        continue;
      }
      registry.fill(HIST("Event/hEventCounter"), 2.0); // selected
      registry.fill(HIST("Event/hNumContrib"), collision.numContrib());

      float hadronicRate = mRateFetcher.fetch(ccdb.service, bc.timestamp(), bc.runNumber(), irSource) * 1.e-3; // kHz

      auto v0s_coll = v0s.sliceBy(perCollision_v0, collision.globalIndex());
      for (const auto& v0 : v0s_coll) {
        // auto o2v0 = v0.template v0_as<aod::V0s>();
        auto pos = v0.template posTrack_as<MyTracks>();
        auto neg = v0.template negTrack_as<MyTracks>();
        // LOGF(info, "v0.globalIndex() = %d, v0.collisionId() = %d, v0.posTrackId() = %d, v0.negTrackId() = %d", v0.globalIndex(), v0.collisionId(), v0.posTrackId(), v0.negTrackId());
        // LOGF(info, "is pos ITSsa = %d", pos.hasITS() && !pos.hasTPC() && !pos.hasTRD() && !pos.hasTOF());
        // LOGF(info, "is neg ITSsa = %d", neg.hasITS() && !neg.hasTPC() && !neg.hasTRD() && !neg.hasTOF());

        if (v0.dcaV0daughters() > v0cuts.cfg_max_dcadau) {
          continue;
        }
        if (v0.v0cosPA() < v0cuts.cfg_min_cospa) {
          continue;
        }
        if (pos.sign() * neg.sign() > 0) {
          continue;
        }
        if (!isSelectedV0Leg(collision, pos) || !isSelectedV0Leg(collision, neg)) {
          continue;
        }

        registry.fill(HIST("V0/hPCA"), v0.dcaV0daughters());
        registry.fill(HIST("V0/hCosPA"), v0.v0cosPA());
        registry.fill(HIST("V0/hAP"), v0.alpha(), v0.qtarm());

        if (v0cuts.cfg_min_qt_strangeness < v0.qtarm()) {
          if (v0cuts.cfg_min_qt_k0s < v0.qtarm()) {
            if (!(v0cuts.cfg_min_mass_lambda_veto < v0.mLambda() && v0.mLambda() < v0cuts.cfg_max_mass_lambda_veto) && !(v0cuts.cfg_min_mass_lambda_veto < v0.mAntiLambda() && v0.mAntiLambda() < v0cuts.cfg_max_mass_lambda_veto)) {
              if ((isPionTightTOFreq(pos) && isSelectedV0LegTight(collision, pos)) && (isPion(neg) && isSelectedV0Leg(collision, neg))) {
                registry.fill(HIST("V0/hMassK0Short"), v0.mK0Short());
                if (v0cuts.cfg_min_mass_k0s < v0.mK0Short() && v0.mK0Short() < v0cuts.cfg_max_mass_k0s) {
                  registry.fill(HIST("V0/hTPCdEdx_P_Pi"), neg.tpcInnerParam(), neg.tpcSignal());
                  registry.fill(HIST("V0/hTOFbeta_P_Pi"), neg.tpcInnerParam(), neg.beta());
                  fillTrackTable(collision, neg, static_cast<uint8_t>(o2::aod::pwgem::dilepton::ml::PID_Label::kPion), hadronicRate);
                }
              }
              if (isPion(pos) && isSelectedV0Leg(collision, pos) && isPionTightTOFreq(neg) && isSelectedV0LegTight(collision, neg)) {
                registry.fill(HIST("V0/hMassK0Short"), v0.mK0Short());
                if (v0cuts.cfg_min_mass_k0s < v0.mK0Short() && v0.mK0Short() < v0cuts.cfg_max_mass_k0s) {
                  registry.fill(HIST("V0/hTPCdEdx_P_Pi"), pos.tpcInnerParam(), pos.tpcSignal());
                  registry.fill(HIST("V0/hTOFbeta_P_Pi"), pos.tpcInnerParam(), pos.beta());
                  fillTrackTable(collision, pos, static_cast<uint8_t>(o2::aod::pwgem::dilepton::ml::PID_Label::kPion), hadronicRate);
                }
              }
            } // end of K0S
          }

          if (!(v0cuts.cfg_min_mass_k0s_veto < v0.mK0Short() && v0.mK0Short() < v0cuts.cfg_max_mass_k0s_veto)) {
            if (isProton(pos) && isSelectedV0Leg(collision, pos) && isPionTight(neg) && isSelectedV0LegTight(collision, neg)) {
              registry.fill(HIST("V0/hMassLambda"), v0.mLambda());
              if (v0cuts.cfg_min_mass_lambda < v0.mLambda() && v0.mLambda() < v0cuts.cfg_max_mass_lambda) {
                fillTrackTable(collision, pos, static_cast<uint8_t>(o2::aod::pwgem::dilepton::ml::PID_Label::kProton), hadronicRate);
                registry.fill(HIST("V0/hTPCdEdx_P_Pr"), pos.tpcInnerParam(), pos.tpcSignal());
                registry.fill(HIST("V0/hTOFbeta_P_Pr"), pos.tpcInnerParam(), pos.beta());
              }
            } // end of Lambda
            if (isPionTight(pos) && isSelectedV0LegTight(collision, pos) && isProton(neg) && isSelectedV0Leg(collision, neg)) {
              registry.fill(HIST("V0/hMassAntiLambda"), v0.mAntiLambda());
              if (v0cuts.cfg_min_mass_lambda < v0.mAntiLambda() && v0.mAntiLambda() < v0cuts.cfg_max_mass_lambda) {
                fillTrackTable(collision, neg, static_cast<uint8_t>(o2::aod::pwgem::dilepton::ml::PID_Label::kProton), hadronicRate);
                registry.fill(HIST("V0/hTPCdEdx_P_Pr"), neg.tpcInnerParam(), neg.tpcSignal());
                registry.fill(HIST("V0/hTOFbeta_P_Pr"), neg.tpcInnerParam(), neg.beta());
              }
            } // end of AntiLambda
          }

        } // end of stangeness

        if (isElectronTight(pos) && isSelectedV0LegTight(collision, pos) && isElectron(neg) && isSelectedV0Leg(collision, neg)) {
          registry.fill(HIST("V0/hMassGamma"), v0.mGamma());
          registry.fill(HIST("V0/hMassGamma_Rxy"), v0.v0radius(), v0.mGamma());
          if (v0cuts.cfg_min_mass_photon < v0.mGamma() && v0.mGamma() < v0cuts.cfg_max_mass_photon) {
            registry.fill(HIST("V0/hXY_Gamma"), v0.x(), v0.y());
            fillTrackTable(collision, neg, static_cast<uint8_t>(o2::aod::pwgem::dilepton::ml::PID_Label::kElectron), hadronicRate);
            registry.fill(HIST("V0/hTPCdEdx_P_El"), neg.tpcInnerParam(), neg.tpcSignal());
            registry.fill(HIST("V0/hTOFbeta_P_El"), neg.tpcInnerParam(), neg.beta());
          }
        } // end of photon conversion

        if (isElectron(pos) && isSelectedV0Leg(collision, pos) && isElectronTight(neg) && isSelectedV0LegTight(collision, neg)) {
          registry.fill(HIST("V0/hMassGamma"), v0.mGamma());
          registry.fill(HIST("V0/hMassGamma_Rxy"), v0.v0radius(), v0.mGamma());
          if (v0cuts.cfg_min_mass_photon < v0.mGamma() && v0.mGamma() < v0cuts.cfg_max_mass_photon) {
            registry.fill(HIST("V0/hXY_Gamma"), v0.x(), v0.y());
            fillTrackTable(collision, pos, static_cast<uint8_t>(o2::aod::pwgem::dilepton::ml::PID_Label::kElectron), hadronicRate);
            registry.fill(HIST("V0/hTPCdEdx_P_El"), pos.tpcInnerParam(), pos.tpcSignal());
            registry.fill(HIST("V0/hTOFbeta_P_El"), pos.tpcInnerParam(), pos.beta());
          }
        } // end of photon conversion

      } // end of v0 loop

      auto cascades_coll = cascades.sliceBy(perCollision_cascade, collision.globalIndex());
      for (const auto& cascade : cascades_coll) {
        // Track casting
        auto bachelor = cascade.template bachelor_as<MyTracks>();
        auto pos = cascade.template posTrack_as<MyTracks>();
        auto neg = cascade.template negTrack_as<MyTracks>();
        if (pos.sign() * neg.sign() > 0) {
          continue;
        }

        if (bachelor.sign() < 0) { // Omega- -> L + K- -> p + pi- + K-
          if (!isProtonTight(pos) || !isPionTight(neg)) {
            continue;
          }
        } else { // Omegabar+ -> Lbar + K+ -> pbar + pi+ + K+
          if (!isProtonTight(neg) || !isPionTight(pos)) {
            continue;
          }
        }

        if (!(cascadecuts.cfg_min_mass_lambda < cascade.mLambda() && cascade.mLambda() < cascadecuts.cfg_max_mass_lambda)) {
          continue;
        }

        if (cascade.cascradius() > cascade.v0radius()) {
          continue;
        }

        if (cascade.dcaV0daughters() > cascadecuts.cfg_max_dcadau_v0) {
          continue;
        }
        if (cascade.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < cascadecuts.cfg_min_cospa_v0) {
          continue;
        }
        if (cascade.v0radius() < cascadecuts.cfg_min_rxy_v0) {
          continue;
        }

        if (cascade.cascradius() < cascadecuts.cfg_min_rxy) {
          continue;
        }

        if (cascade.dcacascdaughters() > cascadecuts.cfg_max_dcadau) {
          continue;
        }
        if (cascade.casccosPA(collision.posX(), collision.posY(), collision.posZ()) < cascadecuts.cfg_min_cospa) {
          continue;
        }

        if (cascade.dcav0topv(collision.posX(), collision.posY(), collision.posZ()) < cascadecuts.cfg_min_dcaxy_v0) {
          continue;
        }

        if (!isSelectedV0LegTight(collision, pos) || !isSelectedV0LegTight(collision, neg) || !isSelectedV0Leg(collision, bachelor)) {
          continue;
        }

        registry.fill(HIST("Cascade/hMassLambda"), cascade.mLambda());
        registry.fill(HIST("Cascade/hV0PCA"), cascade.dcaV0daughters());
        registry.fill(HIST("Cascade/hV0CosPA"), cascade.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
        registry.fill(HIST("Cascade/hPCA"), cascade.dcacascdaughters()); // distance between bachelor and V0.
        registry.fill(HIST("Cascade/hCosPA"), cascade.casccosPA(collision.posX(), collision.posY(), collision.posZ()));

        float length = std::sqrt(std::pow(cascade.x() - collision.posX(), 2) + std::pow(cascade.y() - collision.posY(), 2) + std::pow(cascade.z() - collision.posZ(), 2));
        float mom = cascade.p();
        float ctauXi = length / mom * o2::constants::physics::MassXiMinus;       // 4.91 cm in PDG
        float ctauOmega = length / mom * o2::constants::physics::MassOmegaMinus; // 2.46 cm in PDG

        if (isPion(bachelor)) {
          registry.fill(HIST("Cascade/hMassXi"), cascade.mXi());
          registry.fill(HIST("Cascade/hMassPt_Xi"), cascade.mXi(), cascade.pt());
          registry.fill(HIST("Cascade/hRxy_Xi"), cascade.mXi(), cascade.cascradius());
          registry.fill(HIST("Cascade/hCTau_Xi"), cascade.mXi(), ctauXi);
        }
        if (!(cascadecuts.cfg_min_mass_Xi_veto < cascade.mXi() && cascade.mXi() < cascadecuts.cfg_max_mass_Xi_veto) && isKaon(bachelor)) { // reject Xi candidates
          registry.fill(HIST("Cascade/hMassOmega"), cascade.mOmega());
          registry.fill(HIST("Cascade/hMassPt_Omega"), cascade.mOmega(), cascade.pt());
          registry.fill(HIST("Cascade/hRxy_Omega"), cascade.mOmega(), cascade.cascradius());
          registry.fill(HIST("Cascade/hCTau_Omega"), cascade.mOmega(), ctauOmega);
          if (cascadecuts.cfg_min_mass_Omega < cascade.mOmega() && cascade.mOmega() < cascadecuts.cfg_max_mass_Omega) { // select Omega candidates
            registry.fill(HIST("V0/hTPCdEdx_P_Ka"), bachelor.tpcInnerParam(), bachelor.tpcSignal());
            registry.fill(HIST("V0/hTOFbeta_P_Ka"), bachelor.tpcInnerParam(), bachelor.beta());
            fillTrackTable(collision, bachelor, static_cast<uint8_t>(o2::aod::pwgem::dilepton::ml::PID_Label::kKaon), hadronicRate);
          }
        }
      } // end of cascade loop
    } // end of collision loop
    stored_trackIds.clear();
    stored_trackIds.shrink_to_fit();
  } // end of process

  // please choose only 1 process function.
  void processDummy(filteredMyCollisions const&) {}

  PROCESS_SWITCH(TreeCreatorElectronMLDDA, processPID, "produce ML input for single track level", true); // this is for eID with ITSsa. e/pi/k/p are selected by TOF, and these can be used for ITS-TPC PID.
  PROCESS_SWITCH(TreeCreatorElectronMLDDA, processDummy, "process dummy", false);
};

struct MLTrackQC {
  // Basic checks
  HistogramRegistry registry{
    "registry",
    {
      {"hTPCdEdx_P_All", "TPC dE/dx vs. p;p_{in} (GeV/c);TPC dE/dx", {HistType::kTH2F, {{1000, 0.f, 10.f}, {200, 0, 200}}}},
      {"hTOFbeta_P_All", "TOF beta vs. p;p_{in} (GeV/c);TOF #beta", {HistType::kTH2F, {{1000, 0.f, 10.f}, {220, 0.0, 1.1}}}},
      {"hITSobClusterSize_P_All", "mean ITSob cluster size vs. p;p_{in} (GeV/c);<ITSob cluster size> #times cos(#lambda)", {HistType::kTH2F, {{1000, 0.f, 10.f}, {150, 0.0, 15}}}},
      {"hTPCdEdx_P_Electron", "TPC dE/dx vs. p;p_{in} (GeV/c);TPC dE/dx", {HistType::kTH2F, {{1000, 0.f, 10.f}, {200, 0, 200}}}},
      {"hTOFbeta_P_Electron", "TOF beta vs. p;p_{in} (GeV/c);TOF #beta", {HistType::kTH2F, {{1000, 0.f, 10.f}, {220, 0.0, 1.1}}}},
      {"hITSobClusterSize_P_Electron", "mean ITSob cluster size vs. p;p_{in} (GeV/c);<ITSob cluster size> #times cos(#lambda)", {HistType::kTH2F, {{1000, 0.f, 10.f}, {150, 0.0, 15}}}},
      {"hTPCdEdx_P_Pion", "TPC dE/dx vs. p;p_{in} (GeV/c);TPC dE/dx", {HistType::kTH2F, {{1000, 0.f, 10.f}, {200, 0, 200}}}},
      {"hTOFbeta_P_Pion", "TOF beta vs. p;p_{in} (GeV/c);TOF #beta", {HistType::kTH2F, {{1000, 0.f, 10.f}, {220, 0.0, 1.1}}}},
      {"hITSobClusterSize_P_Pion", "mean ITSob cluster size vs. p;p_{in} (GeV/c);<ITSob cluster size> #times cos(#lambda)", {HistType::kTH2F, {{1000, 0.f, 10.f}, {150, 0.0, 15}}}},
      {"hTPCdEdx_P_Kaon", "TPC dE/dx vs. p;p_{in} (GeV/c);TPC dE/dx", {HistType::kTH2F, {{1000, 0.f, 10.f}, {200, 0, 200}}}},
      {"hTOFbeta_P_Kaon", "TOF beta vs. p;p_{in} (GeV/c);TOF #beta", {HistType::kTH2F, {{1000, 0.f, 10.f}, {220, 0.0, 1.1}}}},
      {"hITSobClusterSize_P_Kaon", "mean ITSob cluster size vs. p;p_{in} (GeV/c);<ITSob cluster size> #times cos(#lambda)", {HistType::kTH2F, {{1000, 0.f, 10.f}, {150, 0.0, 15}}}},
      {"hTPCdEdx_P_Proton", "TPC dE/dx vs. p;p_{in} (GeV/c);TPC dE/dx", {HistType::kTH2F, {{1000, 0.f, 10.f}, {200, 0, 200}}}},
      {"hTOFbeta_P_Proton", "TOF beta vs. p;p_{in} (GeV/c);TOF #beta", {HistType::kTH2F, {{1000, 0.f, 10.f}, {220, 0.0, 1.1}}}},
      {"hITSobClusterSize_P_Proton", "mean ITSob cluster size vs. p;p_{in} (GeV/c);<ITSob cluster size> #times cos(#lambda)", {HistType::kTH2F, {{1000, 0.f, 10.f}, {150, 0.0, 15}}}},

      {"hTPCNsigmaEl_P", "TPC n#sigma_{e} vs. p;p_{in} (GeV/c);n #sigma_{e}^{TPC}", {HistType::kTH2F, {{1000, 0.f, 10.f}, {100, -5, +5}}}},
      {"hTPCNsigmaPi_P", "TPC n#sigma_{#pi} vs. p;p_{in} (GeV/c);n #sigma_{#pi}^{TPC}", {HistType::kTH2F, {{1000, 0.f, 10.f}, {100, -5, +5}}}},
      {"hTPCNsigmaKa_P", "TPC n#sigma_{K} vs. p;p_{in} (GeV/c);n #sigma_{K}^{TPC}", {HistType::kTH2F, {{1000, 0.f, 10.f}, {100, -5, +5}}}},
      {"hTPCNsigmaPr_P", "TPC n#sigma_{p} vs. p;p_{in} (GeV/c);n #sigma_{p}^{TPC}", {HistType::kTH2F, {{1000, 0.f, 10.f}, {100, -5, +5}}}},
      {"hTOFNsigmaEl_P", "TOF n#sigma_{e} vs. p;p_{in} (GeV/c);n #sigma_{e}^{TOF}", {HistType::kTH2F, {{1000, 0.f, 10.f}, {100, -5, +5}}}},
      {"hTOFNsigmaPi_P", "TOF n#sigma_{#pi} vs. p;p_{in} (GeV/c);n #sigma_{#pi}^{TOF}", {HistType::kTH2F, {{1000, 0.f, 10.f}, {100, -5, +5}}}},
      {"hTOFNsigmaKa_P", "TOF n#sigma_{K} vs. p;p_{in} (GeV/c);n #sigma_{K}^{TOF}", {HistType::kTH2F, {{1000, 0.f, 10.f}, {100, -5, +5}}}},
      {"hTOFNsigmaPr_P", "TOF n#sigma_{p} vs. p;p_{in} (GeV/c);n #sigma_{p}^{TOF}", {HistType::kTH2F, {{1000, 0.f, 10.f}, {100, -5, +5}}}},
    },
  };

  using MyPIDTracks = soa::Join<aod::EMTracksForMLPID, aod::EMPIDsEl, aod::EMPIDsPi, aod::EMPIDsKa, aod::EMPIDsPr>;

  void processQC(MyPIDTracks const& tracks)
  {
    for (const auto& track : tracks) {
      registry.fill(HIST("hTPCdEdx_P_All"), track.tpcInnerParam(), track.tpcSignal());
      registry.fill(HIST("hTOFbeta_P_All"), track.tpcInnerParam(), track.beta());
      registry.fill(HIST("hITSobClusterSize_P_All"), track.tpcInnerParam(), track.meanClusterSizeITSob() * std::cos(std::atan(track.tgl())));
      if (track.pidlabel() == static_cast<uint8_t>(o2::aod::pwgem::dilepton::ml::PID_Label::kElectron)) {
        registry.fill(HIST("hTPCdEdx_P_Electron"), track.tpcInnerParam(), track.tpcSignal());
        registry.fill(HIST("hTOFbeta_P_Electron"), track.tpcInnerParam(), track.beta());
        registry.fill(HIST("hITSobClusterSize_P_Electron"), track.tpcInnerParam(), track.meanClusterSizeITSob() * std::cos(std::atan(track.tgl())));
        registry.fill(HIST("hTPCNsigmaEl_P"), track.tpcInnerParam(), track.tpcNSigmaEl());
        registry.fill(HIST("hTOFNsigmaEl_P"), track.tpcInnerParam(), track.tofNSigmaEl());
      } else if (track.pidlabel() == static_cast<uint8_t>(o2::aod::pwgem::dilepton::ml::PID_Label::kPion)) {
        registry.fill(HIST("hTPCdEdx_P_Pion"), track.tpcInnerParam(), track.tpcSignal());
        registry.fill(HIST("hTOFbeta_P_Pion"), track.tpcInnerParam(), track.beta());
        registry.fill(HIST("hITSobClusterSize_P_Pion"), track.tpcInnerParam(), track.meanClusterSizeITSob() * std::cos(std::atan(track.tgl())));
        registry.fill(HIST("hTPCNsigmaPi_P"), track.tpcInnerParam(), track.tpcNSigmaPi());
        registry.fill(HIST("hTOFNsigmaPi_P"), track.tpcInnerParam(), track.tofNSigmaPi());
      } else if (track.pidlabel() == static_cast<uint8_t>(o2::aod::pwgem::dilepton::ml::PID_Label::kKaon)) {
        registry.fill(HIST("hTPCdEdx_P_Kaon"), track.tpcInnerParam(), track.tpcSignal());
        registry.fill(HIST("hTOFbeta_P_Kaon"), track.tpcInnerParam(), track.beta());
        registry.fill(HIST("hITSobClusterSize_P_Kaon"), track.tpcInnerParam(), track.meanClusterSizeITSob() * std::cos(std::atan(track.tgl())));
        registry.fill(HIST("hTPCNsigmaKa_P"), track.tpcInnerParam(), track.tpcNSigmaKa());
        registry.fill(HIST("hTOFNsigmaKa_P"), track.tpcInnerParam(), track.tofNSigmaKa());
      } else if (track.pidlabel() == static_cast<uint8_t>(o2::aod::pwgem::dilepton::ml::PID_Label::kProton)) {
        registry.fill(HIST("hTPCdEdx_P_Proton"), track.tpcInnerParam(), track.tpcSignal());
        registry.fill(HIST("hTOFbeta_P_Proton"), track.tpcInnerParam(), track.beta());
        registry.fill(HIST("hITSobClusterSize_P_Proton"), track.tpcInnerParam(), track.meanClusterSizeITSob() * std::cos(std::atan(track.tgl())));
        registry.fill(HIST("hTPCNsigmaPr_P"), track.tpcInnerParam(), track.tpcNSigmaPr());
        registry.fill(HIST("hTOFNsigmaPr_P"), track.tpcInnerParam(), track.tofNSigmaPr());
      }
    } // end of track loop
  }
  PROCESS_SWITCH(MLTrackQC, processQC, "process QC for single track level", false);

  void processDummy(aod::EMTracksForMLPID const&) {}
  PROCESS_SWITCH(MLTrackQC, processDummy, "process dummy", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TreeCreatorElectronMLDDA>(cfgc, TaskName{"tree-creator-ele-ml-dda"}),
    adaptAnalysisTask<MLTrackQC>(cfgc, TaskName{"ml-track-qc"})};
}
