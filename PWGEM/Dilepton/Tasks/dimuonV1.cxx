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
// This code runs loop over leptons.
//    Please write to: daiki.sekihata@cern.ch

#include "PWGEM/Dilepton/Core/DimuonCut.h"
#include "PWGEM/Dilepton/Core/EMEventCut.h"
#include "PWGEM/Dilepton/DataModel/dileptonTables.h"
#include "PWGEM/Dilepton/Utils/EMFwdTrack.h"
#include "PWGEM/Dilepton/Utils/EMTrack.h"
#include "PWGEM/Dilepton/Utils/EMTrackUtilities.h"
#include "PWGEM/Dilepton/Utils/EventHistograms.h"
#include "PWGEM/Dilepton/Utils/EventMixingHandler.h"
#include "PWGEM/Dilepton/Utils/PairUtilities.h"

#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"

#include <CommonConstants/LHCConstants.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/SliceCache.h>
#include <Framework/runDataProcessing.h>
#include <MathUtils/Utils.h>

#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>

#include <sys/types.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <random>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include <math.h>

struct dimuonV1 {

  // // Configurables
  o2::framework::Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  // o2::framework::Configurable<int> cfgOccupancyEstimator{"cfgOccupancyEstimator", 0, "FT0C:0, Track:1"};
  // o2::framework::Configurable<bool> cfgDoMix{"cfgDoMix", false, "flag for event mixing"};
  // o2::framework::Configurable<int> ndepth{"ndepth", 1000, "depth for event mixing"};
  // o2::framework::Configurable<uint64_t> ndiff_bc_mix{"ndiff_bc_mix", 594, "difference in global BC required in mixed events"};
  // o2::framework::ConfigurableAxis ConfVtxBins{"ConfVtxBins", {o2::framework::VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  // o2::framework::ConfigurableAxis ConfCentBins{"ConfCentBins", {o2::framework::VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.f, 999.f}, "Mixing bins - centrality"};
  // o2::framework::ConfigurableAxis ConfEPBins{"ConfEPBins", {16, -M_PI / 2, +M_PI / 2}, "Mixing bins - event plane angle"};
  // o2::framework::ConfigurableAxis ConfOccupancyBins{"ConfOccupancyBins", {o2::framework::VARIABLE_WIDTH, -1, 1e+10}, "Mixing bins - occupancy"};

  o2::framework::Configurable<bool> cfgApplyWeightTTCA{"cfgApplyWeightTTCA", false, "flag to apply weighting by 1/N"};

  o2::framework::ConfigurableAxis ConfMllBins{"ConfMllBins", {o2::framework::VARIABLE_WIDTH, 0.50, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.70, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.00, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.10, 1.11, 1.12, 1.13, 1.14, 1.15, 1.16, 1.17, 1.18, 1.19, 1.20, 1.3, 1.4, 1.5, 2.00, 2.01, 2.02, 2.03, 2.04, 2.05, 2.06, 2.07, 2.08, 2.09, 2.1, 2.11, 2.12, 2.13, 2.14, 2.15, 2.16, 2.17, 2.18, 2.19, 2.2, 2.21, 2.22, 2.23, 2.24, 2.25, 2.26, 2.27, 2.28, 2.29, 2.3, 2.31, 2.32, 2.33, 2.34, 2.35, 2.36, 2.37, 2.38, 2.39, 2.4, 2.41, 2.42, 2.43, 2.44, 2.45, 2.46, 2.47, 2.48, 2.49, 2.5, 2.51, 2.52, 2.53, 2.54, 2.55, 2.56, 2.57, 2.58, 2.59, 2.6, 2.61, 2.62, 2.63, 2.64, 2.65, 2.66, 2.67, 2.68, 2.69, 2.7, 2.71, 2.72, 2.73, 2.74, 2.75, 2.76, 2.77, 2.78, 2.79, 2.8, 2.81, 2.82, 2.83, 2.84, 2.85, 2.86, 2.87, 2.88, 2.89, 2.9, 2.91, 2.92, 2.93, 2.94, 2.95, 2.96, 2.97, 2.98, 2.99, 3., 3.01, 3.02, 3.03, 3.04, 3.05, 3.06, 3.07, 3.08, 3.09, 3.1, 3.11, 3.12, 3.13, 3.14, 3.15, 3.16, 3.17, 3.18, 3.19, 3.2, 3.21, 3.22, 3.23, 3.24, 3.25, 3.26, 3.27, 3.28, 3.29, 3.3, 3.31, 3.32, 3.33, 3.34, 3.35, 3.36, 3.37, 3.38, 3.39, 3.4, 3.41, 3.42, 3.43, 3.44, 3.45, 3.46, 3.47, 3.48, 3.49, 3.5, 3.51, 3.52, 3.53, 3.54, 3.55, 3.56, 3.57, 3.58, 3.59, 3.6, 3.61, 3.62, 3.63, 3.64, 3.65, 3.66, 3.67, 3.68, 3.69, 3.7, 3.71, 3.72, 3.73, 3.74, 3.75, 3.76, 3.77, 3.78, 3.79, 3.8, 3.81, 3.82, 3.83, 3.84, 3.85, 3.86, 3.87, 3.88, 3.89, 3.9, 3.91, 3.92, 3.93, 3.94, 3.95, 3.96, 3.97, 3.98, 3.99, 4.0, 6.00, 6.10, 6.20, 6.30, 6.40, 6.50, 6.60, 6.70, 6.80, 6.90, 7.00, 7.10, 7.20, 7.30, 7.40, 7.50, 7.60, 7.70, 7.80, 7.90, 8.00, 8.10, 8.20, 8.30, 8.40, 8.50, 8.60, 8.70, 8.80, 8.90, 9.00, 9.10, 9.20, 9.30, 9.40, 9.50, 9.60, 9.70, 9.80, 9.90, 10.00, 10.10, 10.20, 10.30, 10.40, 10.50, 10.60, 10.70, 10.80, 10.90, 11.00, 11.1, 11.2, 11.3, 11.4, 11.5, 11.6, 11.7, 11.8, 11.9, 12.0}, "mmumu bins for output histograms"};

  o2::framework::ConfigurableAxis ConfPtllBins{"ConfPtllBins", {10, 0, 10}, "pTll bins for output histograms"};
  o2::framework::ConfigurableAxis ConfYllBins{"ConfYllBins", {3, -4.0, -2.5}, "yll bins for output histograms"};
  o2::framework::ConfigurableAxis ConfUQBins{"ConfUQBins", {200, -1, 1}, "uQ bins for output histograms"};

  EMEventCut fEMEventCut;
  struct : o2::framework::ConfigurableGroup {
    std::string prefix = "eventcut_group";
    o2::framework::Configurable<float> cfgZvtxMin{"cfgZvtxMin", -10.f, "min. Zvtx"};
    o2::framework::Configurable<float> cfgZvtxMax{"cfgZvtxMax", +10.f, "max. Zvtx"};
    o2::framework::Configurable<bool> cfgRequireSel8{"cfgRequireSel8", false, "require sel8 in event cut"};
    o2::framework::Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "require FT0AND in event cut"};
    o2::framework::Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", true, "require No time frame border in event cut"};
    o2::framework::Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", true, "require no ITS readout frame border in event cut"};
    o2::framework::Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", true, "require no same bunch pileup in event cut"};
    o2::framework::Configurable<bool> cfgRequireVertexITSTPC{"cfgRequireVertexITSTPC", false, "require Vertex ITSTPC in event cut"};             // ITS-TPC matched track contributes PV.
    o2::framework::Configurable<bool> cfgRequireVertexTOFmatched{"cfgRequireVertexTOFmatched", false, "require Vertex TOFmatched in event cut"}; // ITS-TPC-TOF matched track contributes PV.
    o2::framework::Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", true, "require good Zvtx between FT0 vs. PV in event cut"};
    o2::framework::Configurable<int> cfgTrackOccupancyMin{"cfgTrackOccupancyMin", -2, "min. occupancy"};
    o2::framework::Configurable<int> cfgTrackOccupancyMax{"cfgTrackOccupancyMax", 1000000000, "max. occupancy"};
    o2::framework::Configurable<float> cfgFT0COccupancyMin{"cfgFT0COccupancyMin", -2, "min. FT0C occupancy"};
    o2::framework::Configurable<float> cfgFT0COccupancyMax{"cfgFT0COccupancyMax", 1000000000, "max. FT0C occupancy"};

    // o2::framework::Configurable<bool> cfgRequireNoCollInTimeRangeStandard{"cfgRequireNoCollInTimeRangeStandard", false, "require no collision in time range standard"};
    // o2::framework::Configurable<bool> cfgRequireNoCollInTimeRangeStrict{"cfgRequireNoCollInTimeRangeStrict", false, "require no collision in time range strict"};
    // o2::framework::Configurable<bool> cfgRequireNoCollInITSROFStandard{"cfgRequireNoCollInITSROFStandard", false, "require no collision in time range standard"};
    // o2::framework::Configurable<bool> cfgRequireNoCollInITSROFStrict{"cfgRequireNoCollInITSROFStrict", false, "require no collision in time range strict"};
    // o2::framework::Configurable<bool> cfgRequireNoHighMultCollInPrevRof{"cfgRequireNoHighMultCollInPrevRof", false, "require no HM collision in previous ITS ROF"};
    // o2::framework::Configurable<bool> cfgRequireGoodITSLayer3{"cfgRequireGoodITSLayer3", false, "number of inactive chips on ITS layer 3 are below threshold "};
    // o2::framework::Configurable<bool> cfgRequireGoodITSLayer0123{"cfgRequireGoodITSLayer0123", false, "number of inactive chips on ITS layers 0-3 are below threshold "};
    // o2::framework::Configurable<bool> cfgRequireGoodITSLayersAll{"cfgRequireGoodITSLayersAll", false, "number of inactive chips on all ITS layers are below threshold "};
    // for RCT
    o2::framework::Configurable<bool> cfgRequireGoodRCT{"cfgRequireGoodRCT", true, "require good detector flag in run condtion table"};
    o2::framework::Configurable<std::string> cfgRCTLabel{"cfgRCTLabel", "CBT_muon", "select 1 [CBT, CBT_hadronPID, CBT_muon_glo] see O2Physics/Common/CCDB/RCTSelectionFlags.h"};
    o2::framework::Configurable<bool> cfgCheckZDC{"cfgCheckZDC", true, "set ZDC flag for PbPb"};
    o2::framework::Configurable<bool> cfgTreatLimitedAcceptanceAsBad{"cfgTreatLimitedAcceptanceAsBad", false, "reject all events where the detectors relevant for the specified Runlist are flagged as LimitedAcceptance"};

    o2::framework::Configurable<float> cfgCentMin{"cfgCentMin", -1, "min. centrality"};
    o2::framework::Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};
    o2::framework::Configurable<uint16_t> cfgNumContribMin{"cfgNumContribMin", 0, "min. numContrib"};
    o2::framework::Configurable<uint16_t> cfgNumContribMax{"cfgNumContribMax", 65000, "max. numContrib"};
  } eventcuts;

  DimuonCut fDimuonCut;
  struct : o2::framework::ConfigurableGroup {
    std::string prefix = "dimuoncut_group";
    o2::framework::Configurable<float> cfg_min_mass{"cfg_min_mass", 0.0, "min mass"};
    o2::framework::Configurable<float> cfg_max_mass{"cfg_max_mass", 1e+10, "max mass"};
    o2::framework::Configurable<float> cfg_min_pair_pt{"cfg_min_pair_pt", 0.0, "min pair pt"};
    o2::framework::Configurable<float> cfg_max_pair_pt{"cfg_max_pair_pt", 1e+10, "max pair pt"};
    o2::framework::Configurable<float> cfg_min_pair_y{"cfg_min_pair_y", -4.0, "min pair rapidity"};
    o2::framework::Configurable<float> cfg_max_pair_y{"cfg_max_pair_y", -2.5, "max pair rapidity"};
    o2::framework::Configurable<float> cfg_min_pair_dcaxy{"cfg_min_pair_dcaxy", 0.0, "min pair dca3d in sigma"};
    o2::framework::Configurable<float> cfg_max_pair_dcaxy{"cfg_max_pair_dcaxy", 1e+10, "max pair dca3d in sigma"};
    o2::framework::Configurable<bool> cfg_apply_detadphi{"cfg_apply_detadphi", false, "flag to apply deta-dphi elliptic cut"};
    o2::framework::Configurable<float> cfg_min_deta{"cfg_min_deta", 0.02, "min deta between 2 muons (elliptic cut)"};
    o2::framework::Configurable<float> cfg_min_dphi{"cfg_min_dphi", 0.02, "min dphi between 2 muons (elliptic cut)"};

    o2::framework::Configurable<uint8_t> cfg_track_type{"cfg_track_type", 3, "muon track type [0: MFT-MCH-MID, 3: MCH-MID]"};
    o2::framework::Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.8, "min pT for single track"};
    o2::framework::Configurable<float> cfg_max_pt_track{"cfg_max_pt_track", 1e+10, "max pT for single track"};
    o2::framework::Configurable<float> cfg_min_eta_track{"cfg_min_eta_track", -4.0, "min eta for single track"};
    o2::framework::Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", -2.5, "max eta for single track"};
    o2::framework::Configurable<float> cfg_min_phi_track{"cfg_min_phi_track", 0.f, "min phi for single track"};
    o2::framework::Configurable<float> cfg_max_phi_track{"cfg_max_phi_track", 6.3, "max phi for single track"};
    o2::framework::Configurable<int> cfg_min_ncluster_mft{"cfg_min_ncluster_mft", 5, "min ncluster MFT"};
    o2::framework::Configurable<int> cfg_min_ncluster_mch{"cfg_min_ncluster_mch", 5, "min ncluster MCH"};
    o2::framework::Configurable<float> cfg_max_chi2{"cfg_max_chi2", 1e+6, "max chi2/ndf"};
    o2::framework::Configurable<float> cfg_max_chi2mft{"cfg_max_chi2mft", 1e+6, "max chi2/ndf"};
    o2::framework::Configurable<float> cfg_border_pt_for_chi2mchmft{"cfg_border_pt_for_chi2mchmft", 0, "border pt for different max chi2 for MFT-MCH matching"};
    o2::framework::Configurable<float> cfg_max_matching_chi2_mftmch_lowPt{"cfg_max_matching_chi2_mftmch_lowPt", 8, "max chi2 for MFT-MCH matching for low pT"};
    o2::framework::Configurable<float> cfg_max_matching_chi2_mftmch_highPt{"cfg_max_matching_chi2_mftmch_highPt", 40, "max chi2 for MFT-MCH matching for high pT"};
    o2::framework::Configurable<float> cfg_max_matching_chi2_mchmid{"cfg_max_matching_chi2_mchmid", 1e+10, "max chi2 for MCH-MID matching"};
    o2::framework::Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 1e+10, "max dca XY for single track in cm"};
    o2::framework::Configurable<float> cfg_min_rabs{"cfg_min_rabs", 17.6, "min Radius at the absorber end"};
    o2::framework::Configurable<float> cfg_max_rabs{"cfg_max_rabs", 89.5, "max Radius at the absorber end"};
    o2::framework::Configurable<float> cfg_max_diff_chi2_mftmch{"cfg_max_diff_chi2_mftmch", -1.f, "max. diff chi2MatchingMCHMFT between the best and the 2nd best matched candidates"};
    o2::framework::Configurable<bool> enableTTCA{"enableTTCA", false, "Flag to enable or disable TTCA"};
    o2::framework::Configurable<float> cfg_max_relDPt_wrt_matchedMCHMID{"cfg_max_relDPt_wrt_matchedMCHMID", 1e+10f, "max. relative dpt between MFT-MCH-MID and MCH-MID"};
    o2::framework::Configurable<float> cfg_max_DEta_wrt_matchedMCHMID{"cfg_max_DEta_wrt_matchedMCHMID", 1e+10f, "max. deta between MFT-MCH-MID and MCH-MID"};
    o2::framework::Configurable<float> cfg_max_DPhi_wrt_matchedMCHMID{"cfg_max_DPhi_wrt_matchedMCHMID", 1e+10f, "max. dphi between MFT-MCH-MID and MCH-MID"};
    o2::framework::Configurable<bool> requireMFTHitMap{"requireMFTHitMap", false, "flag to apply MFT hit map"};
    o2::framework::Configurable<std::vector<int>> requiredMFTDisks{"requiredMFTDisks", std::vector<int>{0}, "hit map on MFT disks [0,1,2,3,4]. logical-OR of each double-sided disk"};
  } dimuoncuts;

  o2::aod::rctsel::RCTFlagsChecker rctChecker;
  int mRunNumber{0};

  o2::framework::HistogramRegistry fRegistry{"output", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject, false, false};
  static constexpr std::string_view event_pair_types[2] = {"same/", "mix/"};

  void init(o2::framework::InitContext& /*context*/)
  {
    mRunNumber = 0;

    rctChecker.init(eventcuts.cfgRCTLabel.value, eventcuts.cfgCheckZDC.value, eventcuts.cfgTreatLimitedAcceptanceAsBad.value);

    // emh_pos = new MyEMH_muon(ndepth);
    // emh_neg = new MyEMH_muon(ndepth);

    addhistograms();
    DefineEMEventCut();
    DefineDimuonCut();
  }

  template <typename TCollision>
  void initCCDB(TCollision const& collision)
  {
    if (mRunNumber == collision.runNumber()) {
      return;
    }
    mRunNumber = collision.runNumber();
  }

  ~dimuonV1()
  {
    // delete emh_pos;
    // emh_pos = 0x0;
    // delete emh_neg;
    // emh_neg = 0x0;
  }

  void addhistograms()
  {
    // event info
    o2::aod::pwgem::dilepton::utils::eventhistogram::addEventHistograms<-1>(&fRegistry);

    fRegistry.add("Event/before/ZDC/hQxtQxp", "Q_{x}^{t} #upoint Q_{x}^{p} vs. centrality;centrality FT0C (%);Q_{x}^{t} #upoint Q_{x}^{p}", o2::framework::HistType::kTH2D, {{110, 0, 110}, {2000, -1, +1}}, false);
    fRegistry.add("Event/before/ZDC/hQytQyp", "Q_{y}^{t} #upoint Q_{y}^{p} vs. centrality;centrality FT0C (%);Q_{y}^{t} #upoint Q_{y}^{p}", o2::framework::HistType::kTH2D, {{110, 0, 110}, {2000, -1, +1}}, false);
    fRegistry.add("Event/before/ZDC/hQxtQyp", "Q_{x}^{t} #upoint Q_{y}^{p} vs. centrality;centrality FT0C (%);Q_{x}^{t} #upoint Q_{y}^{p}", o2::framework::HistType::kTH2D, {{110, 0, 110}, {2000, -1, +1}}, false);
    fRegistry.add("Event/before/ZDC/hQytQxp", "Q_{y}^{t} #upoint Q_{x}^{p} vs. centrality;centrality FT0C (%);Q_{y}^{t} #upoint Q_{x}^{p}", o2::framework::HistType::kTH2D, {{110, 0, 110}, {2000, -1, +1}}, false);
    fRegistry.addClone("Event/before/ZDC/", "Event/after/ZDC/");

    // pair info
    const o2::framework::AxisSpec axis_mass{ConfMllBins, "m_{#mu#mu} (GeV/c^{2})"};
    const o2::framework::AxisSpec axis_pt{ConfPtllBins, "p_{T,#mu#mu} (GeV/c)"};
    const o2::framework::AxisSpec axis_y{ConfYllBins, "y_{#mu#mu}"};
    const o2::framework::AxisSpec axis_uxQxt{ConfUQBins, "u_{x}Q_{x}^{t}"};
    const o2::framework::AxisSpec axis_uxQxp{ConfUQBins, "u_{x}Q_{x}^{p}"};
    const o2::framework::AxisSpec axis_uyQyt{ConfUQBins, "u_{y}Q_{y}^{t}"};
    const o2::framework::AxisSpec axis_uyQyp{ConfUQBins, "u_{y}Q_{y}^{p}"};

    fRegistry.add("Pair/same/uls/hs", "dilepton", o2::framework::HistType::kTHnSparseD, {axis_mass, axis_pt, axis_y, axis_uxQxt, axis_uxQxp, axis_uyQyt, axis_uyQyp}, true);
    fRegistry.addClone("Pair/same/uls/", "Pair/same/lspp/");
    fRegistry.addClone("Pair/same/uls/", "Pair/same/lsmm/");
    fRegistry.addClone("Pair/same/", "Pair/mix/");
  }

  void DefineEMEventCut()
  {
    fEMEventCut = EMEventCut("fEMEventCut", "fEMEventCut");
    fEMEventCut.SetRequireSel8(eventcuts.cfgRequireSel8);
    fEMEventCut.SetRequireFT0AND(eventcuts.cfgRequireFT0AND);
    fEMEventCut.SetZvtxRange(eventcuts.cfgZvtxMin, eventcuts.cfgZvtxMax);
    fEMEventCut.SetRequireNoTFB(eventcuts.cfgRequireNoTFB);
    fEMEventCut.SetRequireNoITSROFB(eventcuts.cfgRequireNoITSROFB);
    fEMEventCut.SetRequireNoSameBunchPileup(eventcuts.cfgRequireNoSameBunchPileup);
    fEMEventCut.SetRequireVertexITSTPC(eventcuts.cfgRequireVertexITSTPC);
    fEMEventCut.SetRequireVertexTOFmatched(eventcuts.cfgRequireVertexTOFmatched);
    fEMEventCut.SetRequireGoodZvtxFT0vsPV(eventcuts.cfgRequireGoodZvtxFT0vsPV);
    // fEMEventCut.SetRequireNoCollInTimeRangeStandard(eventcuts.cfgRequireNoCollInTimeRangeStandard);
    // fEMEventCut.SetRequireNoCollInTimeRangeStrict(eventcuts.cfgRequireNoCollInTimeRangeStrict);
    // fEMEventCut.SetRequireNoCollInITSROFStandard(eventcuts.cfgRequireNoCollInITSROFStandard);
    // fEMEventCut.SetRequireNoCollInITSROFStrict(eventcuts.cfgRequireNoCollInITSROFStrict);
    // fEMEventCut.SetRequireNoHighMultCollInPrevRof(eventcuts.cfgRequireNoHighMultCollInPrevRof);
    // fEMEventCut.SetRequireGoodITSLayer3(eventcuts.cfgRequireGoodITSLayer3);
    // fEMEventCut.SetRequireGoodITSLayer0123(eventcuts.cfgRequireGoodITSLayer0123);
    // fEMEventCut.SetRequireGoodITSLayersAll(eventcuts.cfgRequireGoodITSLayersAll);
  }

  void DefineDimuonCut()
  {
    fDimuonCut = DimuonCut("fDimuonCut", "fDimuonCut");

    // for pair
    fDimuonCut.SetMassRange(dimuoncuts.cfg_min_mass, dimuoncuts.cfg_max_mass);
    fDimuonCut.SetPairPtRange(dimuoncuts.cfg_min_pair_pt, dimuoncuts.cfg_max_pair_pt);
    fDimuonCut.SetPairYRange(dimuoncuts.cfg_min_pair_y, dimuoncuts.cfg_max_pair_y);
    fDimuonCut.SetPairDCAxyRange(dimuoncuts.cfg_min_pair_dcaxy, dimuoncuts.cfg_max_pair_dcaxy);
    fDimuonCut.SetMindEtadPhi(dimuoncuts.cfg_apply_detadphi, dimuoncuts.cfg_min_deta, dimuoncuts.cfg_min_dphi);

    // for track
    fDimuonCut.SetTrackType(dimuoncuts.cfg_track_type);
    fDimuonCut.SetTrackPtRange(dimuoncuts.cfg_min_pt_track, dimuoncuts.cfg_max_pt_track);
    fDimuonCut.SetTrackEtaRange(dimuoncuts.cfg_min_eta_track, dimuoncuts.cfg_max_eta_track);
    fDimuonCut.SetTrackPhiRange(dimuoncuts.cfg_min_phi_track, dimuoncuts.cfg_max_phi_track);
    fDimuonCut.SetNClustersMFT(dimuoncuts.cfg_min_ncluster_mft, 10);
    fDimuonCut.SetNClustersMCHMID(dimuoncuts.cfg_min_ncluster_mch, 20);
    fDimuonCut.SetChi2(0.f, dimuoncuts.cfg_max_chi2);
    fDimuonCut.SetChi2MFT(0.f, dimuoncuts.cfg_max_chi2mft);
    fDimuonCut.SetMaxDiffMatchingChi2MCHMFT(dimuoncuts.cfg_max_diff_chi2_mftmch);
    fDimuonCut.SetMaxMatchingChi2MCHMFTPtDep([&](float pt) { return (pt < dimuoncuts.cfg_border_pt_for_chi2mchmft ? dimuoncuts.cfg_max_matching_chi2_mftmch_lowPt : dimuoncuts.cfg_max_matching_chi2_mftmch_highPt); });
    fDimuonCut.SetMatchingChi2MCHMID(0.f, dimuoncuts.cfg_max_matching_chi2_mchmid);
    fDimuonCut.SetDCAxy(0.f, dimuoncuts.cfg_max_dcaxy);
    fDimuonCut.SetRabs(dimuoncuts.cfg_min_rabs, dimuoncuts.cfg_max_rabs);
    fDimuonCut.SetMaxPDCARabsDep([&](float rabs) { return (rabs < 26.5 ? 594.f : 324.f); });
    fDimuonCut.SetMaxdPtdEtadPhiwrtMCHMID(dimuoncuts.cfg_max_relDPt_wrt_matchedMCHMID, dimuoncuts.cfg_max_DEta_wrt_matchedMCHMID, dimuoncuts.cfg_max_DPhi_wrt_matchedMCHMID); // this is relevant for global muons
    fDimuonCut.SetMFTHitMap(dimuoncuts.requireMFTHitMap, dimuoncuts.requiredMFTDisks);
    fDimuonCut.EnableTTCA(dimuoncuts.enableTTCA);
  }

  template <int ev_id, typename TCollision, typename TTrack1, typename TTrack2, typename TCut>
  bool fillPairInfo(TCollision const& collision, TTrack1 const& t1, TTrack2 const& t2, TCut const& cut)
  {
    if (!cut.template IsSelectedTrack<false>(t1) || !cut.template IsSelectedTrack<false>(t2)) {
      return false;
    }

    if (!cut.IsSelectedPair(t1, t2)) {
      return false;
    }

    float weight = 1.f;
    if constexpr (ev_id == 0) {
      if (cfgApplyWeightTTCA) {
        weight = map_weight[std::make_pair(t1.globalIndex(), t2.globalIndex())];
      }
    }

    ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassMuon);
    ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassMuon);
    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
    float phi = RecoDecay::constrainAngle(v12.Phi(), 0, 1U);

    // TODO : implement correction for NUA
    float uxQxt = std::cos(1.f * phi) * collision.qxZDCC();
    float uxQxp = std::cos(1.f * phi) * collision.qxZDCA();
    float uyQyt = std::sin(1.f * phi) * collision.qyZDCC();
    float uyQyp = std::sin(1.f * phi) * collision.qyZDCA();

    if (t1.sign() * t2.sign() < 0) { // ULS
      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("uls/hs"), v12.M(), v12.Pt(), v12.Rapidity(), uxQxt, uxQxp, uyQyt, uyQyp, weight);
    } else if (t1.sign() > 0 && t2.sign() > 0) { // LS++
      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hs"), v12.M(), v12.Pt(), v12.Rapidity(), uxQxt, uxQxp, uyQyt, uyQyp, weight);
    } else if (t1.sign() < 0 && t2.sign() < 0) { // LS--
      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hs"), v12.M(), v12.Pt(), v12.Rapidity(), uxQxt, uxQxp, uyQyt, uyQyp, weight);
    }
    return true;
  }

  template <typename TCollisions, typename TLeptons, typename TPresilce, typename TCut>
  void runPairing(TCollisions const& collisions, TLeptons const& posTracks, TLeptons const& negTracks, TPresilce const& perCollision, TCut const& cut)
  {
    for (const auto& collision : collisions) {
      initCCDB(collision);

      float centrality = std::array<float, 3>{collision.centFT0M(), collision.centFT0A(), collision.centFT0C()}[cfgCentEstimator];
      if (centrality < eventcuts.cfgCentMin || eventcuts.cfgCentMax < centrality) {
        continue;
      }
      o2::aod::pwgem::dilepton::utils::eventhistogram::fillEventInfo<0, -1>(&fRegistry, collision);

      float QxtQxp = collision.qxZDCA() * collision.qxZDCC();
      float QytQyp = collision.qyZDCA() * collision.qyZDCC();
      float QytQxp = collision.qxZDCA() * collision.qyZDCC();
      float QxtQyp = collision.qxZDCC() * collision.qyZDCA();

      fRegistry.fill(HIST("Event/before/ZDC/hQxtQxp"), centrality, QxtQxp);
      fRegistry.fill(HIST("Event/before/ZDC/hQytQyp"), centrality, QytQyp);
      fRegistry.fill(HIST("Event/before/ZDC/hQxtQyp"), centrality, QxtQyp);
      fRegistry.fill(HIST("Event/before/ZDC/hQytQxp"), centrality, QytQxp);

      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      if (eventcuts.cfgRequireGoodRCT && !rctChecker.checkTable(collision)) {
        continue;
      }

      o2::aod::pwgem::dilepton::utils::eventhistogram::fillEventInfo<1, -1>(&fRegistry, collision);

      fRegistry.fill(HIST("Event/before/hCollisionCounter"), o2::aod::pwgem::dilepton::utils::eventhistogram::nbin_ev); // accepted
      fRegistry.fill(HIST("Event/after/hCollisionCounter"), o2::aod::pwgem::dilepton::utils::eventhistogram::nbin_ev);  // accepted

      fRegistry.fill(HIST("Event/after/ZDC/hQxtQxp"), centrality, QxtQxp);
      fRegistry.fill(HIST("Event/after/ZDC/hQytQyp"), centrality, QytQyp);
      fRegistry.fill(HIST("Event/after/ZDC/hQxtQyp"), centrality, QxtQyp);
      fRegistry.fill(HIST("Event/after/ZDC/hQytQxp"), centrality, QytQxp);

      auto posTracks_per_coll = posTracks.sliceByCached(perCollision, collision.globalIndex(), cache);
      auto negTracks_per_coll = negTracks.sliceByCached(perCollision, collision.globalIndex(), cache);

      for (const auto& [pos, neg] : combinations(o2::soa::CombinationsFullIndexPolicy(posTracks_per_coll, negTracks_per_coll))) { // ULS
        fillPairInfo<0>(collision, pos, neg, cut);
      }
      // for (const auto& [pos1, pos2] : combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(posTracks_per_coll, posTracks_per_coll))) { // LS++
      //   bool is_pair_ok = fillPairInfo<0>(collision, pos1, pos2, cut);
      //   if (is_pair_ok) {
      //     nlspp++;
      //   }
      // }
      // for (const auto& [neg1, neg2] : combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(negTracks_per_coll, negTracks_per_coll))) { // LS--
      //   bool is_pair_ok = fillPairInfo<0>(collision, neg1, neg2, cut);
      //   if (is_pair_ok) {
      //     nlsmm++;
      //   }
      // }

    } // end of collision loop

  } // end of DF

  template <typename TTrack1, typename TTrack2, typename TCut>
  bool isPairOK(TTrack1 const& t1, TTrack2 const& t2, TCut const& cut)
  {
    if (!cut.IsSelectedTrack(t1) || !cut.IsSelectedTrack(t2)) {
      return false;
    }

    if (!cut.IsSelectedPair(t1, t2)) {
      return false;
    }
    return true;
  }

  std::map<std::pair<int, int>, float> map_weight; // <posId, negId> -> float
  template <typename TCollisions, typename TLeptons, typename TPresilce, typename TCut, typename TAllTracks>
  void fillPairWeightMap(TCollisions const& collisions, TLeptons const& posTracks, TLeptons const& negTracks, TPresilce const& perCollision, TCut const& cut, TAllTracks const& tracks)
  {
    std::vector<std::pair<int, int>> passed_pairIds;
    passed_pairIds.reserve(posTracks.size() * negTracks.size());

    for (const auto& collision : collisions) {
      initCCDB(collision);
      float centrality = std::array<float, 3>{collision.centFT0M(), collision.centFT0A(), collision.centFT0C()}[cfgCentEstimator];
      if (centrality < eventcuts.cfgCentMin || eventcuts.cfgCentMax < centrality) {
        continue;
      }

      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      if (eventcuts.cfgRequireGoodRCT && !rctChecker.checkTable(collision)) {
        continue;
      }

      auto posTracks_per_coll = posTracks.sliceByCached(perCollision, collision.globalIndex(), cache);
      auto negTracks_per_coll = negTracks.sliceByCached(perCollision, collision.globalIndex(), cache);

      for (const auto& [pos, neg] : combinations(o2::soa::CombinationsFullIndexPolicy(posTracks_per_coll, negTracks_per_coll))) { // ULS
        if (isPairOK(pos, neg, cut)) {
          passed_pairIds.emplace_back(std::make_pair(pos.globalIndex(), neg.globalIndex()));
        }
      }
      for (const auto& [pos1, pos2] : combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(posTracks_per_coll, posTracks_per_coll))) { // LS++
        if (isPairOK(pos1, pos2, cut)) {
          passed_pairIds.emplace_back(std::make_pair(pos1.globalIndex(), pos2.globalIndex()));
        }
      }
      for (const auto& [neg1, neg2] : combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(negTracks_per_coll, negTracks_per_coll))) { // LS--
        if (isPairOK(neg1, neg2, cut)) {
          passed_pairIds.emplace_back(std::make_pair(neg1.globalIndex(), neg2.globalIndex()));
        }
      }
    } // end of collision loop

    for (const auto& pairId : passed_pairIds) {
      auto t1 = tracks.rawIteratorAt(std::get<0>(pairId));
      auto t2 = tracks.rawIteratorAt(std::get<1>(pairId));

      float n = 1.f; // include myself.
      for (const auto& ambId1 : t1.ambiguousMuonsIds()) {
        for (const auto& ambId2 : t2.ambiguousMuonsIds()) {
          if (std::find(passed_pairIds.begin(), passed_pairIds.end(), std::make_pair(ambId1, ambId2)) != passed_pairIds.end()) {
            n += 1.f;
          }
        }
      }
      map_weight[pairId] = 1.f / n;
    } // end of passed_pairIds loop
    passed_pairIds.clear();
    passed_pairIds.shrink_to_fit();
  }

  using MyCollisions = o2::soa::Join<o2::aod::EMEvents, o2::aod::EMEventsMult, o2::aod::EMEventsCent, o2::aod::EMEventsZDC>;
  using MyCollision = MyCollisions::iterator;

  using MyMuons = o2::soa::Join<o2::aod::EMPrimaryMuons, o2::aod::EMPrimaryMuonEMEventIds, o2::aod::EMAmbiguousMuonSelfIds>;
  using MyMuon = MyMuons::iterator;
  using FilteredMyMuons = o2::soa::Filtered<MyMuons>;
  using FilteredMyMuon = FilteredMyMuons::iterator;

  using MyEMH_muon = o2::aod::pwgem::dilepton::utils::EventMixingHandler<std::tuple<int, int, int, int>, std::pair<int, int>, o2::aod::pwgem::dilepton::utils::EMFwdTrack>;

  o2::framework::expressions::Filter collisionFilter_centrality = (eventcuts.cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < eventcuts.cfgCentMax) || (eventcuts.cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < eventcuts.cfgCentMax) || (eventcuts.cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < eventcuts.cfgCentMax);
  o2::framework::expressions::Filter collisionFilter_numContrib = eventcuts.cfgNumContribMin <= o2::aod::collision::numContrib && o2::aod::collision::numContrib < eventcuts.cfgNumContribMax;
  o2::framework::expressions::Filter collisionFilter_occupancy_track = eventcuts.cfgTrackOccupancyMin <= o2::aod::evsel::trackOccupancyInTimeRange && o2::aod::evsel::trackOccupancyInTimeRange < eventcuts.cfgTrackOccupancyMax;
  o2::framework::expressions::Filter collisionFilter_occupancy_ft0c = eventcuts.cfgFT0COccupancyMin <= o2::aod::evsel::ft0cOccupancyInTimeRange && o2::aod::evsel::ft0cOccupancyInTimeRange < eventcuts.cfgFT0COccupancyMax;
  using FilteredMyCollisions = o2::soa::Filtered<MyCollisions>;

  o2::framework::SliceCache cache;
  o2::framework::Preslice<o2::aod::EMPrimaryMuonEMEventIds> perCollision_muon = o2::aod::emprimarymuon::emeventId;
  o2::framework::expressions::Filter trackFilter_muon = o2::aod::fwdtrack::trackType == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) && dimuoncuts.cfg_min_pt_track < o2::aod::fwdtrack::pt && o2::aod::fwdtrack::pt < dimuoncuts.cfg_max_pt_track && dimuoncuts.cfg_min_eta_track < o2::aod::fwdtrack::eta && o2::aod::fwdtrack::eta < dimuoncuts.cfg_max_eta_track;
  o2::framework::expressions::Filter ttcaFilter_muon = ifnode(dimuoncuts.enableTTCA.node(), true, o2::aod::emprimarymuon::isAssociatedToMPC == true);

  o2::framework::Partition<FilteredMyMuons> positive_muons = o2::aod::emprimarymuon::sign > int8_t(0);
  o2::framework::Partition<FilteredMyMuons> negative_muons = o2::aod::emprimarymuon::sign < int8_t(0);

  // MyEMH_muon* emh_pos = nullptr;
  // MyEMH_muon* emh_neg = nullptr;

  void processAnalysis(FilteredMyCollisions const& collisions, FilteredMyMuons const& muons)
  {
    if (cfgApplyWeightTTCA) {
      fillPairWeightMap(collisions, positive_muons, negative_muons, o2::aod::emprimarymuon::emeventId, fDimuonCut, muons);
    }
    runPairing(collisions, positive_muons, negative_muons, o2::aod::emprimarymuon::emeventId, fDimuonCut);
  }
  PROCESS_SWITCH(dimuonV1, processAnalysis, "run dimuon v1 analysis", true);

  void processDummy(MyCollisions const&) {}
  PROCESS_SWITCH(dimuonV1, processDummy, "Dummy function", false);
};
o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  return o2::framework::WorkflowSpec{adaptAnalysisTask<dimuonV1>(cfgc, o2::framework::TaskName{"dimuon-v1"})};
}
