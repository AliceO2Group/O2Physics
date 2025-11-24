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

#include "PWGEM/Dilepton/DataModel/dileptonTables.h"
#include "PWGEM/Dilepton/Utils/EMFwdTrack.h"
#include "PWGEM/Dilepton/Utils/EMTrack.h"
#include "PWGEM/Dilepton/Utils/EMTrackUtilities.h"
#include "PWGEM/Dilepton/Utils/EventHistograms.h"
#include "PWGEM/Dilepton/Utils/EventMixingHandler.h"
#include "PWGEM/Dilepton/Utils/PairUtilities.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/LHCConstants.h"
#include "DataFormatsParameters/GRPECSObject.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "MathUtils/Utils.h"

#include "Math/Vector4D.h"

#include <algorithm>
#include <array>
#include <iterator>
#include <map>
#include <random>
#include <ranges>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::dilepton::utils;
using namespace o2::aod::pwgem::dilepton::utils::emtrackutil;
using namespace o2::aod::pwgem::dilepton::utils::pairutil;

using MyEMH_electron = o2::aod::pwgem::dilepton::utils::EventMixingHandler<std::tuple<int, int, int, int>, std::pair<int, int>, EMTrack>;
using MyEMH_muon = o2::aod::pwgem::dilepton::utils::EventMixingHandler<std::tuple<int, int, int, int>, std::pair<int, int>, EMFwdTrack>;
using MyEMH_pair = o2::aod::pwgem::dilepton::utils::EventMixingHandler<std::tuple<int, int, int, int>, std::pair<int, int>, std::tuple<int, int, int, int, EMPair>>;

struct DileptonPolarization {
  // Configurables
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};
  Configurable<float> d_bz_input{"d_bz_input", -999, "bz field in kG, -999 is automatic"};

  Configurable<int> cfgPairType{"cfgPairType", 0, "0:dielectron:0, 1:dimuon"};
  Configurable<int> cfgOccupancyEstimator{"cfgOccupancyEstimator", 0, "FT0C:0, Track:1"};
  Configurable<bool> cfgDoMix{"cfgDoMix", true, "flag for event mixing"};
  Configurable<int> ndepth{"ndepth", 100, "depth for event mixing"};
  Configurable<uint64_t> ndiff_bc_mix{"ndiff_bc_mix", 594, "difference in global BC required in mixed events"};
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis ConfCentBins{"ConfCentBins", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.f, 999.f}, "Mixing bins - centrality"};
  ConfigurableAxis ConfEPBins{"ConfEPBins", {16, -M_PI / 2, +M_PI / 2}, "Mixing bins - event plane angle"};
  ConfigurableAxis ConfOccupancyBins{"ConfOccupancyBins", {VARIABLE_WIDTH, -1, 1e+10}, "Mixing bins - occupancy"};
  Configurable<int> cfgPolarizationFrame{"cfgPolarizationFrame", 0, "frame of polarization. 0:CS, 1:HX, else:FATAL"};

  ConfigurableAxis ConfMllBins{"ConfMllBins", {VARIABLE_WIDTH, 0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.70, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.00, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.10, 1.11, 1.12, 1.13, 1.14, 1.15, 1.16, 1.17, 1.18, 1.19, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.75, 2.80, 2.85, 2.90, 2.95, 3.00, 3.05, 3.10, 3.15, 3.20, 3.25, 3.30, 3.35, 3.40, 3.45, 3.50, 3.55, 3.60, 3.65, 3.70, 3.75, 3.80, 3.85, 3.90, 3.95, 4.00, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.00, 8.10, 8.20, 8.30, 8.40, 8.50, 8.60, 8.70, 8.80, 8.90, 9.00, 9.10, 9.20, 9.30, 9.40, 9.50, 9.60, 9.70, 9.80, 9.90, 10.00, 10.10, 10.20, 10.30, 10.40, 10.50, 10.60, 10.70, 10.80, 10.90, 11.00, 11.1, 11.2, 11.3, 11.4, 11.50, 11.6, 11.7, 11.8, 11.9, 12.0}, "mll bins for output histograms"};
  ConfigurableAxis ConfPtllBins{"ConfPtllBins", {VARIABLE_WIDTH, 0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.50, 3.00, 3.50, 4.00, 4.50, 5.00, 6.00, 7.00, 8.00, 9.00, 10.00}, "pTll bins for output histograms"};
  ConfigurableAxis ConfDCAllBins{"ConfDCAllBins", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0}, "DCAll bins for output histograms"};
  ConfigurableAxis ConfYllBins{"ConYllBins", {1, -1.f, 1.f}, "yll bins for output histograms"}; // pair rapidity

  ConfigurableAxis ConfPolarizationCosThetaBins{"ConfPolarizationCosThetaBins", {20, -1.f, 1.f}, "cos(theta) bins for polarization analysis"};
  ConfigurableAxis ConfPolarizationPhiBins{"ConfPolarizationPhiBins", {1, -M_PI, M_PI}, "phi bins for polarization analysis"};
  ConfigurableAxis ConfPolarizationQuadMomBins{"ConfPolarizationQuadMomBins", {15, -0.5, 1}, "quadrupole moment bins for polarization analysis"}; // quardrupole moment <(3 x cos^2(theta) -1)/2>

  struct : ConfigurableGroup {
    std::string prefix = "eventcut_group";
    Configurable<float> cfgZvtxMin{"cfgZvtxMin", -10.f, "min. Zvtx"};
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", +10.f, "max. Zvtx"};
    // Configurable<bool> cfgRequireSel8{"cfgRequireSel8", true, "require sel8 in event cut"};
    // Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "require FT0AND in event cut"};
    // Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", false, "require No time frame border in event cut"};
    // Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", false, "require no ITS readout frame border in event cut"};
    // Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "require no same bunch pileup in event cut"};
    // Configurable<bool> cfgRequireVertexITSTPC{"cfgRequireVertexITSTPC", false, "require Vertex ITSTPC in event cut"};             // ITS-TPC matched track contributes PV.
    // Configurable<bool> cfgRequireVertexTOFmatched{"cfgRequireVertexTOFmatched", false, "require Vertex TOFmatched in event cut"}; // ITS-TPC-TOF matched track contributes PV.
    // Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", false, "require good Zvtx between FT0 vs. PV in event cut"};
    Configurable<float> cfgCentMin{"cfgCentMin", -1, "min. centrality"};
    Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};
    Configurable<int> cfgTrackOccupancyMin{"cfgTrackOccupancyMin", -2, "min. occupancy"};
    Configurable<int> cfgTrackOccupancyMax{"cfgTrackOccupancyMax", 1000000000, "max. occupancy"};
    Configurable<float> cfgFT0COccupancyMin{"cfgFT0COccupancyMin", -2, "min. FT0C occupancy"};
    Configurable<float> cfgFT0COccupancyMax{"cfgFT0COccupancyMax", 1000000000, "max. FT0C occupancy"};
  } eventcuts;

  struct : ConfigurableGroup {
    std::string prefix = "dileptoncut_group";
    Configurable<float> cfg_min_mass{"cfg_min_mass", 0.0, "min mass"};
    Configurable<float> cfg_max_mass{"cfg_max_mass", 1e+10, "max mass"};
    Configurable<float> cfg_min_pair_pt{"cfg_min_pair_pt", 0.0, "min pair pT"};
    Configurable<float> cfg_max_pair_pt{"cfg_max_pair_pt", 1e+10, "max pair pT"};
    Configurable<float> cfg_min_pair_y{"cfg_min_pair_y", -0.9, "min pair rapidity"};
    Configurable<float> cfg_max_pair_y{"cfg_max_pair_y", +0.9, "max pair rapidity"};

    Configurable<float> cfg_min_track_pt{"cfg_min_track_pt", 0.2, "min pT for single track"};
    Configurable<float> cfg_max_track_pt{"cfg_max_track_pt", 1e+10, "max pT for single track"};
    Configurable<float> cfg_min_track_eta{"cfg_min_track_eta", -0.9, "min eta for single track"};
    Configurable<float> cfg_max_track_eta{"cfg_max_track_eta", +0.9, "max eta for single track"};
  } dileptoncuts;

  struct : ConfigurableGroup {
    std::string prefix = "accBins";
    ConfigurableAxis ConfMllAccBins{"ConfMllAccBins", {40, 0, 4}, "mll bins for acceptance for plarization"};
    ConfigurableAxis ConfPtllAccBins{"ConfPtllAccBins", {100, 0, 10}, "pTll bins for acceptance for plarization"};
    ConfigurableAxis ConfEtallAccBins{"ConEtallAccBins", {30, -1.5f, 1.5f}, "etall bins for acceptance for plarization"};   // pair pseudo-rapidity
    ConfigurableAxis ConfPhillAccBins{"ConPhillAccBins", {36, 0.f, 2 * M_PI}, "phill bins for acceptance for plarization"}; // pair pseudo-rapidity
  } accBins;

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber;
  float d_bz;

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  // static constexpr std::string_view event_cut_types[2] = {"before/", "after/"};
  static constexpr std::string_view event_pair_types[2] = {"same/", "mix/"};

  std::mt19937 engine;
  std::vector<float> cent_bin_edges;
  std::vector<float> zvtx_bin_edges;
  std::vector<float> ep_bin_edges;
  std::vector<float> occ_bin_edges;
  std::vector<float> mll_bin_edges;
  std::vector<float> ptll_bin_edges;
  std::vector<float> etall_bin_edges;
  std::vector<float> phill_bin_edges;

  float leptonM1 = 0.f;
  float leptonM2 = 0.f;

  float beamM1 = o2::constants::physics::MassProton; // mass of beam
  float beamM2 = o2::constants::physics::MassProton; // mass of beam
  float beamE1 = 0.f;                                // beam energy
  float beamE2 = 0.f;                                // beam energy
  float beamP1 = 0.f;                                // beam momentum
  float beamP2 = 0.f;                                // beam momentum

  void init(InitContext& /*context*/)
  {
    mRunNumber = 0;
    d_bz = 0;

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    if (cfgPairType.value == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron)) {
      leptonM1 = o2::constants::physics::MassElectron;
      leptonM2 = o2::constants::physics::MassElectron;
    } else if (cfgPairType.value == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon)) {
      leptonM1 = o2::constants::physics::MassMuon;
      leptonM2 = o2::constants::physics::MassMuon;
    }

    if (ConfVtxBins.value[0] == VARIABLE_WIDTH) {
      zvtx_bin_edges = std::vector<float>(ConfVtxBins.value.begin(), ConfVtxBins.value.end());
      zvtx_bin_edges.erase(zvtx_bin_edges.begin());
      for (const auto& edge : zvtx_bin_edges) {
        LOGF(info, "VARIABLE_WIDTH: zvtx_bin_edges = %f", edge);
      }
    } else {
      int nbins = static_cast<int>(ConfVtxBins.value[0]);
      float xmin = static_cast<float>(ConfVtxBins.value[1]);
      float xmax = static_cast<float>(ConfVtxBins.value[2]);
      zvtx_bin_edges.resize(nbins + 1);
      for (int i = 0; i < nbins + 1; i++) {
        zvtx_bin_edges[i] = (xmax - xmin) / (nbins)*i + xmin;
        LOGF(info, "FIXED_WIDTH: zvtx_bin_edges[%d] = %f", i, zvtx_bin_edges[i]);
      }
    }

    if (ConfCentBins.value[0] == VARIABLE_WIDTH) {
      cent_bin_edges = std::vector<float>(ConfCentBins.value.begin(), ConfCentBins.value.end());
      cent_bin_edges.erase(cent_bin_edges.begin());
      for (const auto& edge : cent_bin_edges) {
        LOGF(info, "VARIABLE_WIDTH: cent_bin_edges = %f", edge);
      }
    } else {
      int nbins = static_cast<int>(ConfCentBins.value[0]);
      float xmin = static_cast<float>(ConfCentBins.value[1]);
      float xmax = static_cast<float>(ConfCentBins.value[2]);
      cent_bin_edges.resize(nbins + 1);
      for (int i = 0; i < nbins + 1; i++) {
        cent_bin_edges[i] = (xmax - xmin) / (nbins)*i + xmin;
        LOGF(info, "FIXED_WIDTH: cent_bin_edges[%d] = %f", i, cent_bin_edges[i]);
      }
    }

    if (ConfEPBins.value[0] == VARIABLE_WIDTH) {
      ep_bin_edges = std::vector<float>(ConfEPBins.value.begin(), ConfEPBins.value.end());
      ep_bin_edges.erase(ep_bin_edges.begin());
      for (const auto& edge : ep_bin_edges) {
        LOGF(info, "VARIABLE_WIDTH: ep_bin_edges = %f", edge);
      }
    } else {
      int nbins = static_cast<int>(ConfEPBins.value[0]);
      float xmin = static_cast<float>(ConfEPBins.value[1]);
      float xmax = static_cast<float>(ConfEPBins.value[2]);
      ep_bin_edges.resize(nbins + 1);
      for (int i = 0; i < nbins + 1; i++) {
        ep_bin_edges[i] = (xmax - xmin) / (nbins)*i + xmin;
        LOGF(info, "FIXED_WIDTH: ep_bin_edges[%d] = %f", i, ep_bin_edges[i]);
      }
    }

    LOGF(info, "cfgOccupancyEstimator = %d", cfgOccupancyEstimator.value);
    if (ConfOccupancyBins.value[0] == VARIABLE_WIDTH) {
      occ_bin_edges = std::vector<float>(ConfOccupancyBins.value.begin(), ConfOccupancyBins.value.end());
      occ_bin_edges.erase(occ_bin_edges.begin());
      for (const auto& edge : occ_bin_edges) {
        LOGF(info, "VARIABLE_WIDTH: occ_bin_edges = %f", edge);
      }
    } else {
      int nbins = static_cast<int>(ConfOccupancyBins.value[0]);
      float xmin = static_cast<float>(ConfOccupancyBins.value[1]);
      float xmax = static_cast<float>(ConfOccupancyBins.value[2]);
      occ_bin_edges.resize(nbins + 1);
      for (int i = 0; i < nbins + 1; i++) {
        occ_bin_edges[i] = (xmax - xmin) / (nbins)*i + xmin;
        LOGF(info, "FIXED_WIDTH: occ_bin_edges[%d] = %f", i, occ_bin_edges[i]);
      }
    }

    emh_pair_uls = new MyEMH_pair(ndepth);
    emh_pair_lspp = new MyEMH_pair(ndepth);
    emh_pair_lsmm = new MyEMH_pair(ndepth);

    if (accBins.ConfMllAccBins.value[0] == VARIABLE_WIDTH) {
      mll_bin_edges = std::vector<float>(accBins.ConfMllAccBins.value.begin(), accBins.ConfMllAccBins.value.end());
      mll_bin_edges.erase(mll_bin_edges.begin());
      for (const auto& edge : mll_bin_edges) {
        LOGF(info, "VARIABLE_WIDTH: mll_bin_edges = %f", edge);
      }
    } else {
      int nbins = static_cast<int>(accBins.ConfMllAccBins.value[0]);
      float xmin = static_cast<float>(accBins.ConfMllAccBins.value[1]);
      float xmax = static_cast<float>(accBins.ConfMllAccBins.value[2]);
      mll_bin_edges.resize(nbins + 1);
      for (int i = 0; i < nbins + 1; i++) {
        mll_bin_edges[i] = (xmax - xmin) / (nbins)*i + xmin;
        LOGF(info, "FIXED_WIDTH: mll_bin_edges[%d] = %f", i, mll_bin_edges[i]);
      }
    }

    if (accBins.ConfPtllAccBins.value[0] == VARIABLE_WIDTH) {
      ptll_bin_edges = std::vector<float>(accBins.ConfPtllAccBins.value.begin(), accBins.ConfPtllAccBins.value.end());
      ptll_bin_edges.erase(ptll_bin_edges.begin());
      for (const auto& edge : ptll_bin_edges) {
        LOGF(info, "VARIABLE_WIDTH: ptll_bin_edges = %f", edge);
      }
    } else {
      int nbins = static_cast<int>(accBins.ConfPtllAccBins.value[0]);
      float xmin = static_cast<float>(accBins.ConfPtllAccBins.value[1]);
      float xmax = static_cast<float>(accBins.ConfPtllAccBins.value[2]);
      ptll_bin_edges.resize(nbins + 1);
      for (int i = 0; i < nbins + 1; i++) {
        ptll_bin_edges[i] = (xmax - xmin) / (nbins)*i + xmin;
        LOGF(info, "FIXED_WIDTH: ptll_bin_edges[%d] = %f", i, ptll_bin_edges[i]);
      }
    }

    if (accBins.ConfEtallAccBins.value[0] == VARIABLE_WIDTH) {
      etall_bin_edges = std::vector<float>(accBins.ConfEtallAccBins.value.begin(), accBins.ConfEtallAccBins.value.end());
      etall_bin_edges.erase(etall_bin_edges.begin());
      for (const auto& edge : etall_bin_edges) {
        LOGF(info, "VARIABLE_WIDTH: etall_bin_edges = %f", edge);
      }
    } else {
      int nbins = static_cast<int>(accBins.ConfEtallAccBins.value[0]);
      float xmin = static_cast<float>(accBins.ConfEtallAccBins.value[1]);
      float xmax = static_cast<float>(accBins.ConfEtallAccBins.value[2]);
      etall_bin_edges.resize(nbins + 1);
      for (int i = 0; i < nbins + 1; i++) {
        etall_bin_edges[i] = (xmax - xmin) / (nbins)*i + xmin;
        LOGF(info, "FIXED_WIDTH: etall_bin_edges[%d] = %f", i, etall_bin_edges[i]);
      }
    }

    if (accBins.ConfPhillAccBins.value[0] == VARIABLE_WIDTH) {
      phill_bin_edges = std::vector<float>(accBins.ConfPhillAccBins.value.begin(), accBins.ConfPhillAccBins.value.end());
      phill_bin_edges.erase(phill_bin_edges.begin());
      for (const auto& edge : phill_bin_edges) {
        LOGF(info, "VARIABLE_WIDTH: phill_bin_edges = %f", edge);
      }
    } else {
      int nbins = static_cast<int>(accBins.ConfPhillAccBins.value[0]);
      float xmin = static_cast<float>(accBins.ConfPhillAccBins.value[1]);
      float xmax = static_cast<float>(accBins.ConfPhillAccBins.value[2]);
      phill_bin_edges.resize(nbins + 1);
      for (int i = 0; i < nbins + 1; i++) {
        phill_bin_edges[i] = (xmax - xmin) / (nbins)*i + xmin;
        LOGF(info, "FIXED_WIDTH: phill_bin_edges[%d] = %f", i, phill_bin_edges[i]);
      }
    }

    std::random_device seed_gen;
    engine = std::mt19937(seed_gen());

    addhistograms();

    fRegistry.add("Pair/mix/hDiffBC", "diff. global BC in mixed event;|BC_{current} - BC_{mixed}|", kTH1D, {{10001, -0.5, 10000.5}}, true);
  }

  template <typename TCollision>
  void initCCDB(TCollision const& collision)
  {
    if (mRunNumber == collision.runNumber()) {
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
      mRunNumber = collision.runNumber();
      return;
    }

    auto run3grp_timestamp = collision.timestamp();
    o2::parameters::GRPObject* grpo = 0x0;
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (!skipGRPOquery)
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      // Fetch magnetic field from ccdb for current collision
      d_bz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      // Fetch magnetic field from ccdb for current collision
      d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kG";
    }
    mRunNumber = collision.runNumber();

    auto grplhcif = ccdb->getForTimeStamp<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", collision.timestamp());
    int beamZ1 = grplhcif->getBeamZ(o2::constants::lhc::BeamC);
    int beamZ2 = grplhcif->getBeamZ(o2::constants::lhc::BeamA);
    int beamA1 = grplhcif->getBeamA(o2::constants::lhc::BeamC);
    int beamA2 = grplhcif->getBeamA(o2::constants::lhc::BeamA);
    beamE1 = grplhcif->getBeamEnergyPerNucleonInGeV(o2::constants::lhc::BeamC);
    beamE2 = grplhcif->getBeamEnergyPerNucleonInGeV(o2::constants::lhc::BeamA);
    beamM1 = o2::constants::physics::MassProton * beamA1;
    beamM2 = o2::constants::physics::MassProton * beamA2;
    beamP1 = std::sqrt(std::pow(beamE1, 2) - std::pow(beamM1, 2));
    beamP2 = std::sqrt(std::pow(beamE2, 2) - std::pow(beamM2, 2));
    LOGF(info, "beamZ1 = %d, beamZ2 = %d, beamA1 = %d, beamA2 = %d, beamE1 = %f (GeV), beamE2 = %f (GeV), beamM1 = %f (GeV), beamM2 = %f (GeV), beamP1 = %f (GeV), beamP2 = %f (GeV)", beamZ1, beamZ2, beamA1, beamA2, beamE1, beamE2, beamM1, beamM2, beamP1, beamP2);
  }

  ~DileptonPolarization()
  {
    delete emh_pair_uls;
    emh_pair_uls = 0x0;
    delete emh_pair_lspp;
    emh_pair_lspp = 0x0;
    delete emh_pair_lsmm;
    emh_pair_lsmm = 0x0;

    map_mixed_eventId_to_globalBC.clear();
  }

  void addhistograms()
  {
    auto hCollisionCounter = fRegistry.add<TH1>("Event/before/hCollisionCounter", "collision counter;;Number of events", kTH1D, {{9, 0.5, 9 + 0.5}}, false);
    hCollisionCounter->GetXaxis()->SetBinLabel(1, "all");
    hCollisionCounter->GetXaxis()->SetBinLabel(2, "FT0AND");
    hCollisionCounter->GetXaxis()->SetBinLabel(3, "No TF border");
    hCollisionCounter->GetXaxis()->SetBinLabel(4, "No ITS ROF border");
    hCollisionCounter->GetXaxis()->SetBinLabel(5, "No Same Bunch Pileup");
    hCollisionCounter->GetXaxis()->SetBinLabel(6, "Is Good Zvtx FT0vsPV");
    hCollisionCounter->GetXaxis()->SetBinLabel(7, "sel8");
    hCollisionCounter->GetXaxis()->SetBinLabel(8, "|Z_{vtx}| < 10 cm");
    hCollisionCounter->GetXaxis()->SetBinLabel(9, "accepted");

    fRegistry.add("Event/before/hZvtx", "vertex z; Z_{vtx} (cm)", kTH1D, {{100, -50, +50}}, false);
    fRegistry.add("Event/before/hCentFT0C", "hCentFT0C;centrality FT0C (%)", kTH1D, {{110, 0, 110}}, false);
    fRegistry.add("Event/before/hCorrOccupancy", "occupancy correlation;FT0C occupancy;track occupancy", kTH2D, {{200, 0, 200000}, {200, 0, 20000}}, false);
    fRegistry.add("Event/before/hEP2_CentFT0C_forMix", "2nd harmonics event plane for mix;centrality FT0C (%);#Psi_{2} (rad.)", kTH2F, {{110, 0, 110}, {180, -M_PI_2, +M_PI_2}}, false);
    fRegistry.addClone("Event/before/", "Event/after/");

    std::string mass_axis_title = "m_{ll} (GeV/c^{2})";
    std::string pair_pt_axis_title = "p_{T,ll} (GeV/c)";
    std::string pair_dca_axis_title = "DCA_{ll} (#sigma)";
    std::string pair_y_axis_title = "y_{ll}";

    if (cfgPairType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron)) {
      mass_axis_title = "m_{ee} (GeV/c^{2})";
      pair_pt_axis_title = "p_{T,ee} (GeV/c)";
      pair_dca_axis_title = "DCA_{ee} (#sigma)";
      pair_y_axis_title = "y_{ee}";
    } else if (cfgPairType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon)) {
      mass_axis_title = "m_{#mu#mu} (GeV/c^{2})";
      pair_pt_axis_title = "p_{T,#mu#mu} (GeV/c)";
      pair_dca_axis_title = "DCA_{#mu#mu} (#sigma)";
      pair_y_axis_title = "y_{#mu#mu}";
    }

    // pair info
    const AxisSpec axis_mass{ConfMllBins, mass_axis_title};
    const AxisSpec axis_pt{ConfPtllBins, pair_pt_axis_title};
    const AxisSpec axis_dca{ConfDCAllBins, pair_dca_axis_title};
    const AxisSpec axis_y{ConfYllBins, pair_y_axis_title};

    std::string frameName = "CS";
    if (cfgPolarizationFrame == 0) {
      frameName = "CS";
    } else if (cfgPolarizationFrame == 1) {
      frameName = "HX";
    } else {
      LOG(fatal) << "set 0 or 1 to cfgPolarizationFrame!";
    }

    const AxisSpec axis_cos_theta{ConfPolarizationCosThetaBins, Form("cos(#theta^{%s})", frameName.data())};
    const AxisSpec axis_phi{ConfPolarizationPhiBins, Form("#varphi^{%s} (rad.)", frameName.data())};
    const AxisSpec axis_quadmom{ConfPolarizationQuadMomBins, Form("#frac{3 cos^{2}(#theta^{%s}) -1}{2}", frameName.data())};
    fRegistry.add("Pair/same/uls/hs", "dilepton", kTHnSparseD, {axis_mass, axis_pt, axis_dca, axis_y, axis_cos_theta, axis_phi, axis_quadmom}, true);

    fRegistry.addClone("Pair/same/uls/", "Pair/same/lspp/");
    fRegistry.addClone("Pair/same/uls/", "Pair/same/lsmm/");
    fRegistry.addClone("Pair/same/", "Pair/mix/");
  }

  template <int ev_id, typename TCollision, typename TDilepton>
  bool fillPairInfo(TCollision const& collision, TDilepton const& dilepton)
  {
    float weight = 1.f;
    ROOT::Math::PtEtaPhiMVector v1(dilepton.pt1(), dilepton.eta1(), dilepton.phi1(), leptonM1);
    ROOT::Math::PtEtaPhiMVector v2(dilepton.pt2(), dilepton.eta2(), dilepton.phi2(), leptonM2);
    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

    if (v12.Rapidity() < dileptoncuts.cfg_min_pair_y || dileptoncuts.cfg_max_pair_y < v12.Rapidity()) {
      return false;
    }

    float pair_dca = pairDCAQuadSum(dilepton.dca1(), dilepton.dca2());
    float cos_thetaPol = 999, phiPol = 999.f;
    auto arrM = std::array<float, 4>{static_cast<float>(v12.Px()), static_cast<float>(v12.Py()), static_cast<float>(v12.Pz()), static_cast<float>(v12.M())};
    auto random_sign = std::pow(-1, engine() % 2); // -1^0 = +1 or -1^1 = -1;
    auto arrD = dilepton.sign1() * dilepton.sign2() < 0 ? (dilepton.sign1() > 0 ? std::array<float, 4>{static_cast<float>(v1.Px()), static_cast<float>(v1.Py()), static_cast<float>(v1.Pz()), leptonM1} : std::array<float, 4>{static_cast<float>(v2.Px()), static_cast<float>(v2.Py()), static_cast<float>(v2.Pz()), leptonM2}) : (random_sign > 0 ? std::array<float, 4>{static_cast<float>(v1.Px()), static_cast<float>(v1.Py()), static_cast<float>(v1.Pz()), leptonM1} : std::array<float, 4>{static_cast<float>(v2.Px()), static_cast<float>(v2.Py()), static_cast<float>(v2.Pz()), leptonM2});
    if (cfgPolarizationFrame == 0) {
      o2::aod::pwgem::dilepton::utils::pairutil::getAngleCS(arrM, arrD, beamE1, beamE2, beamP1, beamP2, cos_thetaPol, phiPol);
    } else if (cfgPolarizationFrame == 1) {
      o2::aod::pwgem::dilepton::utils::pairutil::getAngleHX(arrM, arrD, beamE1, beamE2, beamP1, beamP2, cos_thetaPol, phiPol);
    }
    o2::math_utils::bringToPMPi(phiPol);
    float quadmom = (3.f * std::pow(cos_thetaPol, 2) - 1.f) / 2.f;

    if (dilepton.sign1() * dilepton.sign2() < 0) { // ULS
      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("uls/hs"), v12.M(), v12.Pt(), pair_dca, v12.Rapidity(), cos_thetaPol, phiPol, quadmom, weight);
    } else if (dilepton.sign1() > 0 && dilepton.sign2() > 0) { // LS++
      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hs"), v12.M(), v12.Pt(), pair_dca, v12.Rapidity(), cos_thetaPol, phiPol, quadmom, weight);
    } else if (dilepton.sign1() < 0 && dilepton.sign2() < 0) { // LS--
      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hs"), v12.M(), v12.Pt(), pair_dca, v12.Rapidity(), cos_thetaPol, phiPol, quadmom, weight);
    }

    if constexpr (ev_id == 0) { // same event
      int mbin = lower_bound(mll_bin_edges.begin(), mll_bin_edges.end(), v12.M()) - mll_bin_edges.begin() - 1;
      if (mbin < 0) {
        mbin = 0;
      } else if (static_cast<int>(mll_bin_edges.size()) - 2 < mbin) {
        mbin = static_cast<int>(mll_bin_edges.size()) - 2;
      }

      int ptbin = lower_bound(ptll_bin_edges.begin(), ptll_bin_edges.end(), v12.Pt()) - ptll_bin_edges.begin() - 1;
      if (ptbin < 0) {
        ptbin = 0;
      } else if (static_cast<int>(ptll_bin_edges.size()) - 2 < ptbin) {
        ptbin = static_cast<int>(ptll_bin_edges.size()) - 2;
      }

      int etabin = lower_bound(etall_bin_edges.begin(), etall_bin_edges.end(), v12.Eta()) - etall_bin_edges.begin() - 1;
      if (etabin < 0) {
        etabin = 0;
      } else if (static_cast<int>(etall_bin_edges.size()) - 2 < etabin) {
        etabin = static_cast<int>(etall_bin_edges.size()) - 2;
      }

      float phi12 = v12.Phi();
      o2::math_utils::bringTo02Pi(phi12);
      int phibin = lower_bound(phill_bin_edges.begin(), phill_bin_edges.end(), phi12) - phill_bin_edges.begin() - 1;
      if (phibin < 0) {
        phibin = 0;
      } else if (static_cast<int>(phill_bin_edges.size()) - 2 < phibin) {
        phibin = static_cast<int>(phill_bin_edges.size()) - 2;
      }

      auto key_df_collision = std::make_pair(ndf, collision.globalIndex());
      float phi12_tmp = v12.Phi();
      o2::math_utils::bringTo02Pi(phi12_tmp);
      EMPair empair = EMPair(v12.Pt(), v12.Eta(), phi12_tmp, v12.M(), 0);
      empair.setPositiveLegPxPyPzM(arrD[0], arrD[1], arrD[2], leptonM1);
      // empair.setNegativeLegPtEtaPhiM(t2.pt(), t2.eta(), t2.phi(), leptonM2);
      empair.setPairDCA(pair_dca);
      auto pair_tmp = std::make_tuple(mbin, ptbin, etabin, phibin, empair);
      if (dilepton.sign1() * dilepton.sign2() < 0) { // ULS
        emh_pair_uls->AddTrackToEventPool(key_df_collision, pair_tmp);
      } else if (dilepton.sign1() > 0 && dilepton.sign2() > 0) { // LS++
        emh_pair_lspp->AddTrackToEventPool(key_df_collision, pair_tmp);
      } else if (dilepton.sign1() < 0 && dilepton.sign2() < 0) { // LS--
        emh_pair_lsmm->AddTrackToEventPool(key_df_collision, pair_tmp);
      }
    }
    return true;
  }

  Filter collisionFilter_centrality = eventcuts.cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < eventcuts.cfgCentMax;
  Filter collisionFilter_occupancy_track = eventcuts.cfgTrackOccupancyMin <= o2::aod::evsel::trackOccupancyInTimeRange && o2::aod::evsel::trackOccupancyInTimeRange < eventcuts.cfgTrackOccupancyMax;
  Filter collisionFilter_occupancy_ft0c = eventcuts.cfgFT0COccupancyMin <= o2::aod::evsel::ft0cOccupancyInTimeRange && o2::aod::evsel::ft0cOccupancyInTimeRange < eventcuts.cfgFT0COccupancyMax;
  using filteredCollisions = soa::Filtered<aod::EMThinEvents>;

  Filter dileptonFilter_track1 = (dileptoncuts.cfg_min_track_pt < o2::aod::emdilepton::pt1 && o2::aod::emdilepton::pt1 < dileptoncuts.cfg_max_track_pt) && (dileptoncuts.cfg_min_track_eta < o2::aod::emdilepton::eta1 && o2::aod::emdilepton::eta1 < dileptoncuts.cfg_max_track_eta);
  Filter dileptonFilter_track2 = (dileptoncuts.cfg_min_track_pt < o2::aod::emdilepton::pt2 && o2::aod::emdilepton::pt2 < dileptoncuts.cfg_max_track_pt) && (dileptoncuts.cfg_min_track_eta < o2::aod::emdilepton::eta2 && o2::aod::emdilepton::eta2 < dileptoncuts.cfg_max_track_eta);
  using filteredDileptons = soa::Filtered<aod::EMDileptons>;

  SliceCache cache;
  Preslice<aod::EMDileptons> perCollision = aod::emdilepton::emthineventId;
  Partition<filteredDileptons> dileptonsULS = (o2::aod::emdilepton::sign1 > int16_t(0) && o2::aod::emdilepton::sign2 < int16_t(0)) || (o2::aod::emdilepton::sign1 < int16_t(0) && o2::aod::emdilepton::sign2 > int16_t(0));
  Partition<filteredDileptons> dileptonsLSPP = o2::aod::emdilepton::sign1 > int16_t(0) && o2::aod::emdilepton::sign2 > int16_t(0);
  Partition<filteredDileptons> dileptonsLSMM = o2::aod::emdilepton::sign1 < int16_t(0) && o2::aod::emdilepton::sign2 < int16_t(0);

  MyEMH_pair* emh_pair_uls = nullptr;
  MyEMH_pair* emh_pair_lspp = nullptr;
  MyEMH_pair* emh_pair_lsmm = nullptr;

  std::map<std::pair<int, int>, uint64_t> map_mixed_eventId_to_globalBC;

  int ndf = 0;

  template <typename TCollisions, typename TDileptons>
  void runPairing(TCollisions const& collisions, TDileptons const&)
  {
    for (const auto& collision : collisions) {
      initCCDB(collision);
      float centrality = collision.centFT0C();
      if (centrality < eventcuts.cfgCentMin || eventcuts.cfgCentMax < centrality) {
        continue;
      }

      float ep2 = collision.ep2();
      fRegistry.fill(HIST("Event/after/hZvtx"), collision.posZ());
      fRegistry.fill(HIST("Event/after/hCollisionCounter"), 9); // is qvector calibarated
      fRegistry.fill(HIST("Event/after/hCorrOccupancy"), collision.ft0cOccupancyInTimeRange(), collision.trackOccupancyInTimeRange());
      fRegistry.fill(HIST("Event/after/hEP2_CentFT0C_forMix"), collision.centFT0C(), ep2);

      auto dileptons_uls_per_coll = dileptonsULS->sliceByCached(aod::emdilepton::emthineventId, collision.globalIndex(), cache);
      auto dileptons_lspp_per_coll = dileptonsLSPP->sliceByCached(aod::emdilepton::emthineventId, collision.globalIndex(), cache);
      auto dileptons_lsmm_per_coll = dileptonsLSMM->sliceByCached(aod::emdilepton::emthineventId, collision.globalIndex(), cache);
      // LOGF(info, "collision.globalIndex() = %d, dileptons_uls_per_coll.size() = %d, dileptons_lspp_per_coll.size() = %d, dileptons_lsmm_per_coll.size() = %d", collision.globalIndex(), dileptons_uls_per_coll.size(), dileptons_lspp_per_coll.size(), dileptons_lsmm_per_coll.size());

      int nuls = 0, nlspp = 0, nlsmm = 0;
      for (const auto& dilepton : dileptons_uls_per_coll) { // ULS
        bool is_pair_ok = fillPairInfo<0>(collision, dilepton);
        if (is_pair_ok) {
          nuls++;
        }
      }
      for (const auto& dilepton : dileptons_lspp_per_coll) { // LS++
        bool is_pair_ok = fillPairInfo<0>(collision, dilepton);
        if (is_pair_ok) {
          nlspp++;
        }
      }
      for (const auto& dilepton : dileptons_lsmm_per_coll) { // LS--
        bool is_pair_ok = fillPairInfo<0>(collision, dilepton);
        if (is_pair_ok) {
          nlsmm++;
        }
      }

      if (!cfgDoMix || !(nuls > 0 || nlspp > 0 || nlsmm > 0)) {
        continue;
      }

      // event mixing
      int zbin = lower_bound(zvtx_bin_edges.begin(), zvtx_bin_edges.end(), collision.posZ()) - zvtx_bin_edges.begin() - 1;
      if (zbin < 0) {
        zbin = 0;
      } else if (static_cast<int>(zvtx_bin_edges.size()) - 2 < zbin) {
        zbin = static_cast<int>(zvtx_bin_edges.size()) - 2;
      }

      int centbin = lower_bound(cent_bin_edges.begin(), cent_bin_edges.end(), centrality) - cent_bin_edges.begin() - 1;
      if (centbin < 0) {
        centbin = 0;
      } else if (static_cast<int>(cent_bin_edges.size()) - 2 < centbin) {
        centbin = static_cast<int>(cent_bin_edges.size()) - 2;
      }

      int epbin = lower_bound(ep_bin_edges.begin(), ep_bin_edges.end(), ep2) - ep_bin_edges.begin() - 1;
      if (epbin < 0) {
        epbin = 0;
      } else if (static_cast<int>(ep_bin_edges.size()) - 2 < epbin) {
        epbin = static_cast<int>(ep_bin_edges.size()) - 2;
      }

      int occbin = -1;
      if (cfgOccupancyEstimator == 0) {
        occbin = lower_bound(occ_bin_edges.begin(), occ_bin_edges.end(), collision.ft0cOccupancyInTimeRange()) - occ_bin_edges.begin() - 1;
      } else if (cfgOccupancyEstimator == 1) {
        occbin = lower_bound(occ_bin_edges.begin(), occ_bin_edges.end(), collision.trackOccupancyInTimeRange()) - occ_bin_edges.begin() - 1;
      } else {
        occbin = lower_bound(occ_bin_edges.begin(), occ_bin_edges.end(), collision.ft0cOccupancyInTimeRange()) - occ_bin_edges.begin() - 1;
      }

      if (occbin < 0) {
        occbin = 0;
      } else if (static_cast<int>(occ_bin_edges.size()) - 2 < occbin) {
        occbin = static_cast<int>(occ_bin_edges.size()) - 2;
      }

      // LOGF(info, "collision.globalIndex() = %d, collision.posZ() = %f, centrality = %f, ep2 = %f, collision.ft0cOccupancyInTimeRange() = %f, zbin = %d, centbin = %d, epbin = %d, occbin = %d", collision.globalIndex(), collision.posZ(), centrality, ep2, collision.ft0cOccupancyInTimeRange(), zbin, centbin, epbin, occbin);

      std::tuple<int, int, int, int> key_bin = std::make_tuple(zbin, centbin, epbin, occbin);
      std::pair<int, int> key_df_collision = std::make_pair(ndf, collision.globalIndex()); // this gives the current event.

      auto selected_pairs_uls_in_this_event = emh_pair_uls->GetTracksPerCollision(key_df_collision);
      auto selected_pairs_lspp_in_this_event = emh_pair_lspp->GetTracksPerCollision(key_df_collision);
      auto selected_pairs_lsmm_in_this_event = emh_pair_lsmm->GetTracksPerCollision(key_df_collision);
      auto collisionIds_in_mixing_pool = emh_pair_uls->GetCollisionIdsFromEventPool(key_bin);
      float weight = 1.f;

      for (const auto& mix_dfId_collisionId : collisionIds_in_mixing_pool) {
        auto pairs_uls_from_event_pool = emh_pair_uls->GetTracksPerCollision(mix_dfId_collisionId);
        auto pairs_lspp_from_event_pool = emh_pair_lspp->GetTracksPerCollision(mix_dfId_collisionId);
        auto pairs_lsmm_from_event_pool = emh_pair_lsmm->GetTracksPerCollision(mix_dfId_collisionId);

        auto globalBC_mix = map_mixed_eventId_to_globalBC[mix_dfId_collisionId];
        uint64_t diffBC = std::max(collision.globalBC(), globalBC_mix) - std::min(collision.globalBC(), globalBC_mix);
        fRegistry.fill(HIST("Pair/mix/hDiffBC"), diffBC);
        if (diffBC < ndiff_bc_mix) {
          continue;
        }

        for (const auto& pair1 : selected_pairs_uls_in_this_event) { // ULS mix
          auto empair1 = std::get<4>(pair1);
          auto v_pos = empair1.getPositiveLeg(); // pt, eta, phi, M
          auto arrD = std::array<float, 4>{static_cast<float>(v_pos.Px()), static_cast<float>(v_pos.Py()), static_cast<float>(v_pos.Pz()), leptonM1};

          int mbin = lower_bound(mll_bin_edges.begin(), mll_bin_edges.end(), empair1.mass()) - mll_bin_edges.begin() - 1;
          if (mbin < 0) {
            mbin = 0;
          } else if (static_cast<int>(mll_bin_edges.size()) - 2 < mbin) {
            mbin = static_cast<int>(mll_bin_edges.size()) - 2;
          }

          int ptbin = lower_bound(ptll_bin_edges.begin(), ptll_bin_edges.end(), empair1.pt()) - ptll_bin_edges.begin() - 1;
          if (ptbin < 0) {
            ptbin = 0;
          } else if (static_cast<int>(ptll_bin_edges.size()) - 2 < ptbin) {
            ptbin = static_cast<int>(ptll_bin_edges.size()) - 2;
          }

          int etabin = lower_bound(etall_bin_edges.begin(), etall_bin_edges.end(), empair1.eta()) - etall_bin_edges.begin() - 1;
          if (etabin < 0) {
            etabin = 0;
          } else if (static_cast<int>(etall_bin_edges.size()) - 2 < etabin) {
            etabin = static_cast<int>(etall_bin_edges.size()) - 2;
          }

          int phibin = lower_bound(phill_bin_edges.begin(), phill_bin_edges.end(), empair1.phi()) - phill_bin_edges.begin() - 1;
          if (phibin < 0) {
            phibin = 0;
          } else if (static_cast<int>(phill_bin_edges.size()) - 2 < phibin) {
            phibin = static_cast<int>(phill_bin_edges.size()) - 2;
          }

          for (const auto& pair2 : std::views::filter(pairs_uls_from_event_pool, [&mbin, &ptbin, &etabin, &phibin](std::tuple<int, int, int, int, EMPair> t) { return std::get<0>(t) == mbin && std::get<1>(t) == ptbin && std::get<2>(t) == etabin && std::get<3>(t) == phibin; })) {
            auto empair2 = std::get<4>(pair2);
            auto arrM = std::array<float, 4>{static_cast<float>(empair2.px()), static_cast<float>(empair2.py()), static_cast<float>(empair2.pz()), static_cast<float>(empair2.mass())};

            float cos_thetaPol = 999, phiPol = 999.f;
            if (cfgPolarizationFrame == 0) {
              o2::aod::pwgem::dilepton::utils::pairutil::getAngleCS(arrM, arrD, beamE1, beamE2, beamP1, beamP2, cos_thetaPol, phiPol);
            } else if (cfgPolarizationFrame == 1) {
              o2::aod::pwgem::dilepton::utils::pairutil::getAngleHX(arrM, arrD, beamE1, beamE2, beamP1, beamP2, cos_thetaPol, phiPol);
            }
            o2::math_utils::bringToPMPi(phiPol);
            float quadmom = (3.f * std::pow(cos_thetaPol, 2) - 1.f) / 2.f;
            fRegistry.fill(HIST("Pair/mix/uls/hs"), empair1.mass(), empair1.pt(), empair1.getPairDCA(), empair1.rapidity(), cos_thetaPol, phiPol, quadmom, weight);
          }
        } // end of ULS

        for (const auto& pair1 : selected_pairs_lspp_in_this_event) { // LS++
          auto empair1 = std::get<4>(pair1);
          auto v_pos = empair1.getPositiveLeg(); // pt, eta, phi, M
          auto arrD = std::array<float, 4>{static_cast<float>(v_pos.Px()), static_cast<float>(v_pos.Py()), static_cast<float>(v_pos.Pz()), leptonM1};

          int mbin = lower_bound(mll_bin_edges.begin(), mll_bin_edges.end(), empair1.mass()) - mll_bin_edges.begin() - 1;
          if (mbin < 0) {
            mbin = 0;
          } else if (static_cast<int>(mll_bin_edges.size()) - 2 < mbin) {
            mbin = static_cast<int>(mll_bin_edges.size()) - 2;
          }

          int ptbin = lower_bound(ptll_bin_edges.begin(), ptll_bin_edges.end(), empair1.pt()) - ptll_bin_edges.begin() - 1;
          if (ptbin < 0) {
            ptbin = 0;
          } else if (static_cast<int>(ptll_bin_edges.size()) - 2 < ptbin) {
            ptbin = static_cast<int>(ptll_bin_edges.size()) - 2;
          }

          int etabin = lower_bound(etall_bin_edges.begin(), etall_bin_edges.end(), empair1.eta()) - etall_bin_edges.begin() - 1;
          if (etabin < 0) {
            etabin = 0;
          } else if (static_cast<int>(etall_bin_edges.size()) - 2 < etabin) {
            etabin = static_cast<int>(etall_bin_edges.size()) - 2;
          }

          int phibin = lower_bound(phill_bin_edges.begin(), phill_bin_edges.end(), empair1.phi()) - phill_bin_edges.begin() - 1;
          if (phibin < 0) {
            phibin = 0;
          } else if (static_cast<int>(phill_bin_edges.size()) - 2 < phibin) {
            phibin = static_cast<int>(phill_bin_edges.size()) - 2;
          }

          for (const auto& pair2 : std::views::filter(pairs_lspp_from_event_pool, [&mbin, &ptbin, &etabin, &phibin](std::tuple<int, int, int, int, EMPair> t) { return std::get<0>(t) == mbin && std::get<1>(t) == ptbin && std::get<2>(t) == etabin && std::get<3>(t) == phibin; })) {
            auto empair2 = std::get<4>(pair2);
            auto arrM = std::array<float, 4>{static_cast<float>(empair2.px()), static_cast<float>(empair2.py()), static_cast<float>(empair2.pz()), static_cast<float>(empair2.mass())};

            float cos_thetaPol = 999, phiPol = 999.f;
            if (cfgPolarizationFrame == 0) {
              o2::aod::pwgem::dilepton::utils::pairutil::getAngleCS(arrM, arrD, beamE1, beamE2, beamP1, beamP2, cos_thetaPol, phiPol);
            } else if (cfgPolarizationFrame == 1) {
              o2::aod::pwgem::dilepton::utils::pairutil::getAngleHX(arrM, arrD, beamE1, beamE2, beamP1, beamP2, cos_thetaPol, phiPol);
            }
            o2::math_utils::bringToPMPi(phiPol);
            float quadmom = (3.f * std::pow(cos_thetaPol, 2) - 1.f) / 2.f;
            fRegistry.fill(HIST("Pair/mix/lspp/hs"), empair1.mass(), empair1.pt(), empair1.getPairDCA(), empair1.rapidity(), cos_thetaPol, phiPol, quadmom, weight);
          }
        } // end of LS++

        for (const auto& pair1 : selected_pairs_lsmm_in_this_event) { // LS--
          auto empair1 = std::get<4>(pair1);
          auto v_pos = empair1.getPositiveLeg(); // pt, eta, phi, M
          auto arrD = std::array<float, 4>{static_cast<float>(v_pos.Px()), static_cast<float>(v_pos.Py()), static_cast<float>(v_pos.Pz()), leptonM1};

          int mbin = lower_bound(mll_bin_edges.begin(), mll_bin_edges.end(), empair1.mass()) - mll_bin_edges.begin() - 1;
          if (mbin < 0) {
            mbin = 0;
          } else if (static_cast<int>(mll_bin_edges.size()) - 2 < mbin) {
            mbin = static_cast<int>(mll_bin_edges.size()) - 2;
          }

          int ptbin = lower_bound(ptll_bin_edges.begin(), ptll_bin_edges.end(), empair1.pt()) - ptll_bin_edges.begin() - 1;
          if (ptbin < 0) {
            ptbin = 0;
          } else if (static_cast<int>(ptll_bin_edges.size()) - 2 < ptbin) {
            ptbin = static_cast<int>(ptll_bin_edges.size()) - 2;
          }

          int etabin = lower_bound(etall_bin_edges.begin(), etall_bin_edges.end(), empair1.eta()) - etall_bin_edges.begin() - 1;
          if (etabin < 0) {
            etabin = 0;
          } else if (static_cast<int>(etall_bin_edges.size()) - 2 < etabin) {
            etabin = static_cast<int>(etall_bin_edges.size()) - 2;
          }

          int phibin = lower_bound(phill_bin_edges.begin(), phill_bin_edges.end(), empair1.phi()) - phill_bin_edges.begin() - 1;
          if (phibin < 0) {
            phibin = 0;
          } else if (static_cast<int>(phill_bin_edges.size()) - 2 < phibin) {
            phibin = static_cast<int>(phill_bin_edges.size()) - 2;
          }

          for (const auto& pair2 : std::views::filter(pairs_lsmm_from_event_pool, [&mbin, &ptbin, &etabin, &phibin](std::tuple<int, int, int, int, EMPair> t) { return std::get<0>(t) == mbin && std::get<1>(t) == ptbin && std::get<2>(t) == etabin && std::get<3>(t) == phibin; })) {
            auto empair2 = std::get<4>(pair2);
            auto arrM = std::array<float, 4>{static_cast<float>(empair2.px()), static_cast<float>(empair2.py()), static_cast<float>(empair2.pz()), static_cast<float>(empair2.mass())};

            float cos_thetaPol = 999, phiPol = 999.f;
            if (cfgPolarizationFrame == 0) {
              o2::aod::pwgem::dilepton::utils::pairutil::getAngleCS(arrM, arrD, beamE1, beamE2, beamP1, beamP2, cos_thetaPol, phiPol);
            } else if (cfgPolarizationFrame == 1) {
              o2::aod::pwgem::dilepton::utils::pairutil::getAngleHX(arrM, arrD, beamE1, beamE2, beamP1, beamP2, cos_thetaPol, phiPol);
            }
            o2::math_utils::bringToPMPi(phiPol);
            float quadmom = (3.f * std::pow(cos_thetaPol, 2) - 1.f) / 2.f;
            fRegistry.fill(HIST("Pair/mix/lsmm/hs"), empair1.mass(), empair1.pt(), empair1.getPairDCA(), empair1.rapidity(), cos_thetaPol, phiPol, quadmom, weight);
          }
        } // end of LS--

      } // end of loop over mixed event pool

      if (nuls > 0 || nlspp > 0 || nlsmm > 0) {
        map_mixed_eventId_to_globalBC[key_df_collision] = collision.globalBC();
        emh_pair_uls->AddCollisionIdAtLast(key_bin, key_df_collision);
        emh_pair_lspp->AddCollisionIdAtLast(key_bin, key_df_collision);
        emh_pair_lsmm->AddCollisionIdAtLast(key_bin, key_df_collision);
      } // end of if pair exist

    } // end of collision loop

  } // end of DF

  void processAnalysis(filteredCollisions const& collisions, filteredDileptons const& dileptons)
  {
    runPairing(collisions, dileptons);
    ndf++;
  }
  PROCESS_SWITCH(DileptonPolarization, processAnalysis, "run dilepton analysis", true);

  void processDummy(aod::EMThinEvents const&) {}
  PROCESS_SWITCH(DileptonPolarization, processDummy, "Dummy function", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<DileptonPolarization>(cfgc, TaskName{"dilepton-polarization"})};
}
