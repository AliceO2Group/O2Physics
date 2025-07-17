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

/// \file   zdcQVectors.cxx
/// \author Noor Koster
/// \since  11/2024
/// \brief  In this task the energy calibration and recentring of Q-vectors constructed in the ZDCs will be done

#include <stdlib.h>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <vector>
#include <typeinfo>
#include <memory>
#include <string>

#include "CCDB/BasicCCDBManager.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StaticFor.h"

#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "ReconstructionDataFormats/Track.h"
#include "PWGCF/DataModel/SPTableZDC.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TROOT.h"

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::track;
using namespace o2::aod::evsel;

namespace o2::analysis::qvectortask
{
int counter = 0;

// Define histogrm names here to use same names for creating and later uploading and retrieving data from ccdb
// Energy calibration:
std::vector<TString> namesEcal(10, "");
std::vector<std::vector<TString>> names(5, std::vector<TString>()); //(1x 4d 4x 1d)
std::vector<TString> vnames = {"hvertex_vx", "hvertex_vy"};

// https://alice-notes.web.cern.ch/system/files/notes/analysis/620/017-May-31-analysis_note-ALICE_analysis_note_v2.pdf
std::vector<double> pxZDC = {-1.75, 1.75, -1.75, 1.75};
std::vector<double> pyZDC = {-1.75, -1.75, 1.75, 1.75};
double alphaZDC = 0.395;

// q-vectors before (q) and after (qRec) recentering.
std::vector<double> q(4); // start values of [QxA, QyA, QxC, QyC]
std::vector<double> qNoEq(4); // start values of [QxA, QyA, QxC, QyC]

// for energy calibration
std::vector<double> eZN(8);      // uncalibrated energy for the 2x4 towers (a1, a2, a3, a4, c1, c2, c3, c4)
std::vector<double> meanEZN(10); // mean energies from calibration histos (common A, t1-4 A,common C, t1-4C)
std::vector<double> e(8, 0.);    // calibrated energies (a1, a2, a3, a4, c1, c2, c3, c4))

//  Define variables needed to do the recentring steps.
double centrality = 0;
int runnumber = 0;
int lastRunNumber = 0;
std::vector<double> v(3, 0); // vx, vy, vz
bool isSelected = true;

} // namespace o2::analysis::qvectortask

using namespace o2::analysis::qvectortask;

struct ZdcQVectors {

  Produces<aod::SPTableZDC> spTableZDC;

  ConfigurableAxis axisCent{"axisCent", {90, 0, 90}, "Centrality axis in 1% bins"};
  ConfigurableAxis axisCent10{"axisCent10", {9, 0, 90}, "Centrality axis in 10% bins"};
  ConfigurableAxis axisQ{"axisQ", {100, -2, 2}, "Q vector (xy) in ZDC"};
  ConfigurableAxis axisVxBig{"axisVxBig", {3, -0.01, 0.01}, "for Pos X of collision"};
  ConfigurableAxis axisVyBig{"axisVyBig", {3, -0.01, 0.01}, "for Pos Y of collision"};
  ConfigurableAxis axisVzBig{"axisVzBig", {3, -10, 10}, "for Pos Z of collision"};
  ConfigurableAxis axisVx{"axisVx", {100, -0.01, 0.01}, "for Pos X of collision"};
  ConfigurableAxis axisVy{"axisVy", {100, -0.01, 0.01}, "for Pos Y of collision"};
  ConfigurableAxis axisVz{"axisVz", {100, -10, 10}, "for vz of collision"};

  // Centrality Estimators -> standard is FT0C
  O2_DEFINE_CONFIGURABLE(cfgFT0Cvariant1, bool, false, "Set centrality estimator to cfgFT0Cvariant1");
  O2_DEFINE_CONFIGURABLE(cfgFT0M, bool, false, "Set centrality estimator to cfgFT0M");
  O2_DEFINE_CONFIGURABLE(cfgFV0A, bool, false, "Set centrality estimator to cfgFV0A");
  O2_DEFINE_CONFIGURABLE(cfgNGlobal, bool, false, "Set centrality estimator to cfgNGlobal");

  O2_DEFINE_CONFIGURABLE(cfgVtxZ, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgMagField, float, 99999, "Configurable magnetic field; default CCDB will be queried")
  O2_DEFINE_CONFIGURABLE(cfgEnergyCal, std::string, "Users/c/ckoster/ZDC/LHC23_zzh_pass4/Energy", "ccdb path for energy calibration histos")
  O2_DEFINE_CONFIGURABLE(cfgMeanv, std::string, "Users/c/ckoster/ZDC/LHC23_zzh_pass4/vmean", "ccdb path for mean v histos")
  O2_DEFINE_CONFIGURABLE(cfgMinEntriesSparseBin, int, 100, "Minimal number of entries allowed in 4D recentering histogram to use for recentering.")
  O2_DEFINE_CONFIGURABLE(cfgRec, std::string, "Users/c/ckoster/ZDC/LHC23_PbPb_pass4", "ccdb path for recentering histos");
  O2_DEFINE_CONFIGURABLE(cfgFillCommonRegistry, bool, true, "Fill common registry with histograms");

  // Additional event selections
  O2_DEFINE_CONFIGURABLE(cfgEvSelsMaxOccupancy, int, 10000, "Maximum occupancy of selected events");
  O2_DEFINE_CONFIGURABLE(cfgEvSelsNoSameBunchPileupCut, bool, true, "kNoSameBunchPileupCut");
  O2_DEFINE_CONFIGURABLE(cfgEvSelsIsGoodZvtxFT0vsPV, bool, true, "kIsGoodZvtxFT0vsPV");
  O2_DEFINE_CONFIGURABLE(cfgEvSelsNoCollInTimeRangeStandard, bool, true, "kNoCollInTimeRangeStandard");
  O2_DEFINE_CONFIGURABLE(cfgEvSelsDoOccupancySel, bool, true, "Bool for event selection on detector occupancy");
  O2_DEFINE_CONFIGURABLE(cfgEvSelsIsVertexITSTPC, bool, true, "Selects collisions with at least one ITS-TPC track");
  O2_DEFINE_CONFIGURABLE(cfgEvSelsIsGoodITSLayersAll, bool, true, "Cut time intervals with dead ITS staves");
  O2_DEFINE_CONFIGURABLE(cfgEvSelsCentMin, float, 0, "Minimum cenrality for selected events");
  O2_DEFINE_CONFIGURABLE(cfgEvSelsCentMax, float, 90, "Maximum cenrality for selected events");

  // define my.....
  // Filter collisionFilter = nabs(aod::collision::posZ) < cfgVtxZ;
  using UsedCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::CentFT0Ms, aod::CentFV0As, aod::CentNGlobals>;
  using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;

  enum SelectionCriteria {
    evSel_FilteredEvent,
    evSel_Zvtx,
    evSel_sel8,
    evSel_occupancy,
    evSel_kNoSameBunchPileup,
    evSel_kIsGoodZvtxFT0vsPV,
    evSel_kNoCollInTimeRangeStandard,
    evSel_kIsVertexITSTPC,
    evSel_CentCuts,
    evSel_kIsGoodITSLayersAll,
    evSel_isSelectedZDC,
    nEventSelections
  };

  enum CalibModes {
    kEnergyCal,
    kMeanv,
    kRec
  };

  //  Define output
  HistogramRegistry registry{"Registry"};

  Service<ccdb::BasicCCDBManager> ccdb;

  // keep track of calibration histos for each given step and iteration
  struct Calib {
    std::vector<TList*> calibList = std::vector<TList*>(3, nullptr); // [0] Enerfy cal, [1] vmean, [2] recentering
    std::vector<bool> calibfilesLoaded = std::vector<bool>(3, false);
    int atStep = 0;
    int atIteration = 0;
  } cal;

  enum FillType {
    kBefore,
    kAfter
  };

  void init(InitContext const&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    std::vector<const char*> sides = {"A", "C"};
    std::vector<const char*> capCOORDS = {"X", "Y"};

    if (cfgFillCommonRegistry) {
      registry.add<TH2>(Form("QA/before/hSPplaneA"), "hSPplaneA", kTH2D, {{100, -4, 4}, axisCent10});
      registry.add<TH2>(Form("QA/before/hSPplaneC"), "hSPplaneC", kTH2D, {{100, -4, 4}, axisCent10});
      registry.add<TH2>(Form("QA/before/hSPplaneFull"), "hSPplaneFull", kTH2D, {{100, -4, 4}, axisCent10});
      for (const auto& side : sides) {
        registry.add<TH2>(Form("QA/before/hZN%s_Qx_vs_Qy", side), Form("hZN%s_Qx_vs_Qy", side), kTH2F, {axisQ, axisQ});
      }

      for (const auto& COORD1 : capCOORDS) {
        for (const auto& COORD2 : capCOORDS) {
          // Now we get: <XX> <XY> & <YX> <YY> vs. Centrality
          registry.add<TProfile>(Form("QA/before/hQ%sA_Q%sC_vs_cent", COORD1, COORD2), Form("hQ%sA_Q%sC_vs_cent", COORD1, COORD2), kTProfile, {axisCent});
          registry.add<TProfile>(Form("QA/before/hQ%sA_Q%sC_vs_vx", COORD1, COORD2), Form("hQ%sA_Q%sC_vs_vx", COORD1, COORD2), kTProfile, {axisVx});
          registry.add<TProfile>(Form("QA/before/hQ%sA_Q%sC_vs_vy", COORD1, COORD2), Form("hQ%sA_Q%sC_vs_vy", COORD1, COORD2), kTProfile, {axisVy});
          registry.add<TProfile>(Form("QA/before/hQ%sA_Q%sC_vs_vz", COORD1, COORD2), Form("hQ%sA_Q%sC_vs_vz", COORD1, COORD2), kTProfile, {axisVz});
        }
      }

      // Add histograms for each step in the calibration process.
      // Sides is {A,C} and capcoords is {X,Y}
      for (const auto& side : sides) {
        for (const auto& coord : capCOORDS) {
          registry.add(Form("QA/before/hQ%s%s_vs_cent", coord, side), Form("hQ%s%s_vs_cent", coord, side), {HistType::kTProfile, {axisCent10}});
          registry.add(Form("QA/before/hQ%s%s_vs_vx", coord, side), Form("hQ%s%s_vs_vx", coord, side), {HistType::kTProfile, {axisVx}});
          registry.add(Form("QA/before/hQ%s%s_vs_vy", coord, side), Form("hQ%s%s_vs_vy", coord, side), {HistType::kTProfile, {axisVy}});
          registry.add(Form("QA/before/hQ%s%s_vs_vz", coord, side), Form("hQ%s%s_vs_vz", coord, side), {HistType::kTProfile, {axisVz}});

          names[0].push_back(TString::Format("hQ%s%s_mean_Cent_V_run", coord, side));
          names[1].push_back(TString::Format("hQ%s%s_mean_cent_run", coord, side));
          names[2].push_back(TString::Format("hQ%s%s_mean_vx_run", coord, side));
          names[3].push_back(TString::Format("hQ%s%s_mean_vy_run", coord, side));
          names[4].push_back(TString::Format("hQ%s%s_mean_vz_run", coord, side));
        } // end of capCOORDS
      } // end of sides

      registry.add<TH1>("QA/centrality_before", "centrality_before", kTH1D, {{100, 0, 100}});
      registry.add<TH1>("QA/centrality_after", "centrality_after", kTH1D, {{100, 0, 100}});

      registry.add<TProfile>("QA/ZNA_Energy", "ZNA_Energy", kTProfile, {{8, 0, 8}});
      registry.add<TProfile>("QA/ZNC_Energy", "ZNC_Energy", kTProfile, {{8, 0, 8}});

      registry.add<TProfile>("QA/before/ZNA_pmC", "ZNA_pmC", kTProfile, {{1, 0, 1.}});
      registry.add<TProfile>("QA/before/ZNA_pm1", "ZNA_pm1", kTProfile, {{1, 0, 1.}});
      registry.add<TProfile>("QA/before/ZNA_pm2", "ZNA_pm2", kTProfile, {{1, 0, 1.}});
      registry.add<TProfile>("QA/before/ZNA_pm3", "ZNA_pm3", kTProfile, {{1, 0, 1.}});
      registry.add<TProfile>("QA/before/ZNA_pm4", "ZNA_pm4", kTProfile, {{1, 0, 1.}});

      registry.add<TProfile>("QA/before/ZNC_pmC", "ZNC_pmC", kTProfile, {{1, 0, 1.}});
      registry.add<TProfile>("QA/before/ZNC_pm1", "ZNC_pm1", kTProfile, {{1, 0, 1.}});
      registry.add<TProfile>("QA/before/ZNC_pm2", "ZNC_pm2", kTProfile, {{1, 0, 1.}});
      registry.add<TProfile>("QA/before/ZNC_pm3", "ZNC_pm3", kTProfile, {{1, 0, 1.}});
      registry.add<TProfile>("QA/before/ZNC_pm4", "ZNC_pm4", kTProfile, {{1, 0, 1.}});

      registry.add<TProfile>("QA/before/ZNA_Qx", "ZNA_Qx", kTProfile, {{1, 0, 1.}});
      registry.add<TProfile>("QA/before/ZNA_Qy", "ZNA_Qy", kTProfile, {{1, 0, 1.}});
      registry.add<TProfile>("QA/before/ZNC_Qx", "ZNC_Qx", kTProfile, {{1, 0, 1.}});
      registry.add<TProfile>("QA/before/ZNC_Qy", "ZNC_Qy", kTProfile, {{1, 0, 1.}});

      registry.add<TH2>("QA/before/ZNA_Qx_vs_Centrality", "ZNA_Qx_vs_Centrality", kTH2D, {{100, 0, 100}, {200, -2, 2}});
      registry.add<TH2>("QA/before/ZNA_Qy_vs_Centrality", "ZNA_Qy_vs_Centrality", kTH2D, {{100, 0, 100}, {200, -2, 2}});
      registry.add<TH2>("QA/before/ZNC_Qx_vs_Centrality", "ZNC_Qx_vs_Centrality", kTH2D, {{100, 0, 100}, {200, -2, 2}});
      registry.add<TH2>("QA/before/ZNC_Qy_vs_Centrality", "ZNC_Qy_vs_Centrality", kTH2D, {{100, 0, 100}, {200, -2, 2}});

      registry.add<TH2>("QA/before/ZNA_pmC_vs_Centrality", "ZNA_pmC_vs_Centrality", kTH2D, {{100, 0, 100}, {300, 0, 300}});
      registry.add<TH2>("QA/before/ZNA_pmSUM_vs_Centrality", "ZNA_pmSUM_vs_Centrality", kTH2D, {{100, 0, 100}, {300, 0, 300}});
      registry.add<TH2>("QA/before/ZNA_pm1_vs_Centrality", "ZNA_pm1_vs_Centrality", kTH2D, {{100, 0, 100}, {100, 0, 1}});
      registry.add<TH2>("QA/before/ZNA_pm2_vs_Centrality", "ZNA_pm2_vs_Centrality", kTH2D, {{100, 0, 100}, {100, 0, 1}});
      registry.add<TH2>("QA/before/ZNA_pm3_vs_Centrality", "ZNA_pm3_vs_Centrality", kTH2D, {{100, 0, 100}, {100, 0, 1}});
      registry.add<TH2>("QA/before/ZNA_pm4_vs_Centrality", "ZNA_pm4_vs_Centrality", kTH2D, {{100, 0, 100}, {100, 0, 1}});

      registry.add<TH2>("QA/before/ZNC_pmC_vs_Centrality", "ZNC_pmC_vs_Centrality", kTH2D, {{100, 0, 100}, {300, 0, 300}});
      registry.add<TH2>("QA/before/ZNC_pmSUM_vs_Centrality", "ZNC_pmSUM_vs_Centrality", kTH2D, {{100, 0, 100}, {300, 0, 300}});
      registry.add<TH2>("QA/before/ZNC_pm1_vs_Centrality", "ZNC_pm1_vs_Centrality", kTH2D, {{100, 0, 100}, {100, 0, 1}});
      registry.add<TH2>("QA/before/ZNC_pm2_vs_Centrality", "ZNC_pm2_vs_Centrality", kTH2D, {{100, 0, 100}, {100, 0, 1}});
      registry.add<TH2>("QA/before/ZNC_pm3_vs_Centrality", "ZNC_pm3_vs_Centrality", kTH2D, {{100, 0, 100}, {100, 0, 1}});
      registry.add<TH2>("QA/before/ZNC_pm4_vs_Centrality", "ZNC_pm4_vs_Centrality", kTH2D, {{100, 0, 100}, {100, 0, 1}});

      registry.addClone("QA/before/", "QA/after/");

      registry.add<TProfile>("QA/before/ZNA_Qx_noEq", "ZNA_Qx_noEq", kTProfile, {{1, 0, 1.}});
      registry.add<TProfile>("QA/before/ZNA_Qy_noEq", "ZNA_Qy_noEq", kTProfile, {{1, 0, 1.}});
      registry.add<TProfile>("QA/before/ZNC_Qx_noEq", "ZNC_Qx_noEq", kTProfile, {{1, 0, 1.}});
      registry.add<TProfile>("QA/before/ZNC_Qy_noEq", "ZNC_Qy_noEq", kTProfile, {{1, 0, 1.}});
    }

    // Tower mean energies vs. centrality used for tower gain equalisation
    int totalTowers = 10;
    int totalTowersPerSide = 5;
    for (int tower = 0; tower < totalTowers; tower++) {
      namesEcal[tower] = TString::Format("hZN%s_mean_t%i_cent", sides[(tower < totalTowersPerSide) ? 0 : 1], tower % 5);
      registry.add<TProfile2D>(Form("Energy/%s", namesEcal[tower].Data()), Form("%s", namesEcal[tower].Data()), kTProfile2D, {{1, 0, 1}, axisCent});
    }

    // recentered q-vectors (to check what steps are finished in the end)
    registry.add<TProfile>("vmean/hvertex_vx", "hvertex_vx", kTProfile, {{1, 0., 1.}});
    registry.add<TProfile>("vmean/hvertex_vy", "hvertex_vy", kTProfile, {{1, 0., 1.}});
    registry.add<TProfile>("vmean/hvertex_vz", "hvertex_vz", kTProfile, {{1, 0., 1.}});

    registry.add("hEventCount", "Number of Event; Cut; #Events Passed Cut", {HistType::kTH1D, {{nEventSelections, 0, nEventSelections}}});
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_FilteredEvent + 1, "Filtered events");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_Zvtx + 1, "Z vertex cut event");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_sel8 + 1, "Sel8");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_occupancy + 1, "kOccupancy");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_kNoSameBunchPileup + 1, "kNoSameBunchPileup");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_kIsGoodZvtxFT0vsPV + 1, "kIsGoodZvtxFT0vsPV");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_kNoCollInTimeRangeStandard + 1, "kNoCollInTimeRangeStandard");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_kIsVertexITSTPC + 1, "kIsVertexITSTPC");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_CentCuts + 1, "Cenrality range");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_kIsGoodITSLayersAll + 1, "kkIsGoodITSLayersAll");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_isSelectedZDC + 1, "isSelected");
  }

  template <typename TCollision>
  bool eventSelected(TCollision collision, const float& centrality)
  {
    if (std::fabs(collision.posZ()) > cfgVtxZ)
      return 0;
    registry.fill(HIST("hEventCount"), evSel_Zvtx);

    if (!collision.sel8())
      return 0;
    registry.fill(HIST("hEventCount"), evSel_sel8);

    // Occupancy
    if (cfgEvSelsDoOccupancySel) {
      auto occupancy = collision.trackOccupancyInTimeRange();
      if (occupancy > cfgEvSelsMaxOccupancy) {
        return 0;
      }
      registry.fill(HIST("hEventCount"), evSel_occupancy);
    }

    if (cfgEvSelsNoSameBunchPileupCut) {
      if (!collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
        // rejects collisions which are associated with the same "found-by-T0" bunch crossing
        // https://indico.cern.ch/event/1396220/#1-event-selection-with-its-rof
        return 0;
      }
      registry.fill(HIST("hEventCount"), evSel_kNoSameBunchPileup);
    }
    if (cfgEvSelsIsGoodZvtxFT0vsPV) {
      if (!collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
        // removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference
        // use this cut at low multiplicities with caution
        return 0;
      }
      registry.fill(HIST("hEventCount"), evSel_kIsGoodZvtxFT0vsPV);
    }
    if (cfgEvSelsNoCollInTimeRangeStandard) {
      if (!collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
        //  Rejection of the collisions which have other events nearby
        return 0;
      }
      registry.fill(HIST("hEventCount"), evSel_kNoCollInTimeRangeStandard);
    }

    if (cfgEvSelsIsVertexITSTPC) {
      if (!collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
        // selects collisions with at least one ITS-TPC track, and thus rejects vertices built from ITS-only tracks
        return 0;
      }
      registry.fill(HIST("hEventCount"), evSel_kIsVertexITSTPC);
    }

    if (centrality > cfgEvSelsCentMax || centrality < cfgEvSelsCentMin)
      return 0;
    registry.fill(HIST("hEventCount"), evSel_CentCuts);

    if (cfgEvSelsIsGoodITSLayersAll) {
      if (!collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
        // New event selection bits to cut time intervals with dead ITS staves
        // https://indico.cern.ch/event/1493023/ (09-01-2025)
        return 0;
      }
      registry.fill(HIST("hEventCount"), evSel_kIsGoodITSLayersAll);
    }

    return 1;
  }

  template <FillType ft>
  inline void fillCommonRegistry(double qxa, double qya, double qxc, double qyc, std::vector<double> v, double centrality)
  {
    // loop for filling multiple histograms with different naming patterns
    //  Always fill the uncentered "raw" Q-vector histos!
    static constexpr std::string_view Time[] = {"before", "after"};

    registry.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/hZNA_Qx_vs_Qy"), qxa, qya);
    registry.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/hZNC_Qx_vs_Qy"), qxc, qyc);

    registry.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/hQXA_QXC_vs_cent"), centrality, qxa * qxc);
    registry.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/hQYA_QYC_vs_cent"), centrality, qya * qyc);
    registry.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/hQYA_QXC_vs_cent"), centrality, qya * qxc);
    registry.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/hQXA_QYC_vs_cent"), centrality, qxa * qyc);

    registry.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/hQXA_vs_cent"), centrality, qxa);
    registry.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/hQYA_vs_cent"), centrality, qya);
    registry.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/hQXC_vs_cent"), centrality, qxc);
    registry.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/hQYC_vs_cent"), centrality, qyc);

    registry.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/hQXA_vs_vx"), v[0], qxa);
    registry.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/hQYA_vs_vx"), v[0], qya);
    registry.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/hQXC_vs_vx"), v[0], qxc);
    registry.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/hQYC_vs_vx"), v[0], qyc);

    registry.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/hQXA_QXC_vs_vx"), v[0], qxa * qxc);
    registry.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/hQYA_QYC_vs_vx"), v[0], qya * qyc);
    registry.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/hQYA_QXC_vs_vx"), v[0], qya * qxc);
    registry.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/hQXA_QYC_vs_vx"), v[0], qxa * qyc);

    registry.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/hQXA_vs_vy"), v[1], qxa);
    registry.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/hQYA_vs_vy"), v[1], qya);
    registry.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/hQXC_vs_vy"), v[1], qxc);
    registry.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/hQYC_vs_vy"), v[1], qyc);

    registry.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/hQXA_QXC_vs_vy"), v[1], qxa * qxc);
    registry.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/hQYA_QYC_vs_vy"), v[1], qya * qyc);
    registry.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/hQYA_QXC_vs_vy"), v[1], qya * qxc);
    registry.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/hQXA_QYC_vs_vy"), v[1], qxa * qyc);

    registry.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/hQXA_vs_vz"), v[2], qxa);
    registry.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/hQYA_vs_vz"), v[2], qya);
    registry.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/hQXC_vs_vz"), v[2], qxc);
    registry.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/hQYC_vs_vz"), v[2], qyc);

    registry.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/hQXA_QXC_vs_vz"), v[2], qxa * qxc);
    registry.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/hQYA_QYC_vs_vz"), v[2], qya * qyc);
    registry.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/hQYA_QXC_vs_vz"), v[2], qya * qxc);
    registry.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/hQXA_QYC_vs_vz"), v[2], qxa * qyc);

    registry.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/ZNA_Qx_vs_Centrality"), centrality, qxa);
    registry.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/ZNA_Qy_vs_Centrality"), centrality, qya);
    registry.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/ZNC_Qx_vs_Centrality"), centrality, qxc);
    registry.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/ZNC_Qy_vs_Centrality"), centrality, qyc);

    // add psi!!
    double psiA = 1.0 * std::atan2(qxc, qxa);
    registry.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/hSPplaneA"), psiA, centrality, 1);
    double psiC = 1.0 * std::atan2(qyc, qya);
    registry.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/hSPplaneC"), psiC, centrality, 1);
    double psiFull = 1.0 * std::atan2(qxc + qyc, qxa + qya);
    registry.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/hSPplaneFull"), psiFull, centrality, 1);
  }

  template <CalibModes cm>
  void loadCalibrations(uint64_t timestamp, std::string ccdb_dir)
  {
    // iteration = 0 (Energy calibration) -> step 0 only
    // iteration 1,2,3,4,5 = recentering -> 5 steps per iteration (1x 4D + 4x 1D)

    if (cal.calibfilesLoaded[cm])
      return;

    if (ccdb_dir.empty() == false) {
      cal.calibList[cm] = ccdb->getForTimeStamp<TList>(ccdb_dir, timestamp);
      cal.calibfilesLoaded[cm] = true;
      LOGF(info, "Loaded calibration histos from %s", ccdb_dir.c_str());
      if (cm == kRec) {
        cal.atStep = 5;
        cal.atIteration = 5;
      }
    }
  }

  template <typename T, CalibModes cm>
  double getCorrection(const char* objName, int iteration = 0, int step = 0)
  {
    T* hist = nullptr;
    double calibConstant{0};

    if (cm == kEnergyCal) {
      TList* list = cal.calibList[cm];
      hist = reinterpret_cast<T*>(list->FindObject(Form("%s", objName)));
    } else if (cm == kMeanv) {
      TList* list = cal.calibList[cm];
      hist = reinterpret_cast<T*>(list->FindObject(Form("%s", objName)));
    } else if (cm == kRec) {
      TList* list = reinterpret_cast<TList*>(cal.calibList[cm]->FindObject(Form("it%i_step%i", iteration, step)));
      if (!list) {
        LOGF(fatal, "No calibration list for iteration %i and step %i", iteration, step);
      }
      hist = reinterpret_cast<T*>(list->FindObject(Form("%s", objName)));
      if (!hist) {
        LOGF(fatal, "No calibration histo for iteration %i and step %i -> %s", iteration, step, objName);
      }
      cal.atStep = step;
      cal.atIteration = iteration;
    }

    if (!hist) {
      LOGF(fatal, "%s not available.. Abort..", objName);
    }

    if (hist->InheritsFrom("TProfile2D")) {
      // needed for energy calibration!
      TProfile2D* h = reinterpret_cast<TProfile2D*>(hist);
      TString name = h->GetName();
      int binrunnumber = h->GetXaxis()->FindBin(TString::Format("%d", runnumber));
      int bin = h->GetYaxis()->FindBin(centrality);
      calibConstant = h->GetBinContent(binrunnumber, bin);
    } else if (hist->InheritsFrom("TProfile")) {
      TProfile* h = reinterpret_cast<TProfile*>(hist);
      TString name = h->GetName();
      int bin{};
      if (name.Contains("mean_vx"))
        bin = h->GetXaxis()->FindBin(v[0]);
      if (name.Contains("mean_vy"))
        bin = h->GetXaxis()->FindBin(v[1]);
      if (name.Contains("mean_vz"))
        bin = h->GetXaxis()->FindBin(v[2]);
      if (name.Contains("mean_cent"))
        bin = h->GetXaxis()->FindBin(centrality);
      if (name.Contains("vertex"))
        bin = h->GetXaxis()->FindBin(TString::Format("%i", runnumber));
      calibConstant = h->GetBinContent(bin);
    } else if (hist->InheritsFrom("THnSparse")) {
      std::vector<int> sparsePars;
      THnSparseD* h = reinterpret_cast<THnSparseD*>(hist);
      sparsePars.push_back(h->GetAxis(0)->FindBin(centrality));
      sparsePars.push_back(h->GetAxis(1)->FindBin(v[0]));
      sparsePars.push_back(h->GetAxis(2)->FindBin(v[1]));
      sparsePars.push_back(h->GetAxis(3)->FindBin(v[2]));

      for (std::size_t i = 0; i < sparsePars.size(); i++) {
        h->GetAxis(i)->SetRange(sparsePars[i], sparsePars[i]);
      }

      TH1D* tempProj = h->Projection(4);
      calibConstant = tempProj->GetMean();

      if (tempProj->GetEntries() < cfgMinEntriesSparseBin) {
        LOGF(debug, "1 entry in sparse bin! Not used... (increase binsize)");
        calibConstant = 0;
        isSelected = false;
      }

      delete tempProj;
    }

    return calibConstant;
  }

  void process(UsedCollisions::iterator const& collision,
               BCsRun3 const& /*bcs*/,
               aod::Zdcs const& /*zdcs*/)
  {
    // for Q-vector calculation
    //  A[0] & C[1]
    std::vector<double> sumZN(2, 0.), sumZN_noEq(2, 0.);
    std::vector<double> xEnZN(2, 0.), xEnZN_noEq(2, 0.);
    std::vector<double> yEnZN(2, 0.), yEnZN_noEq(2, 0.);

    isSelected = true;

    // TODO Implement other ZDC estimators
    auto cent = collision.centFT0C();
    if (cfgFT0Cvariant1)
      cent = collision.centFT0CVariant1();
    if (cfgFT0M)
      cent = collision.centFT0M();
    if (cfgFV0A)
      cent = collision.centFV0A();
    if (cfgNGlobal)
      cent = collision.centNGlobal();

    if (cfgFillCommonRegistry)
      registry.fill(HIST("QA/centrality_before"), cent);

    registry.fill(HIST("hEventCount"), evSel_FilteredEvent);

    if (!eventSelected(collision, cent)) {
      // event not selected
      isSelected = false;
      spTableZDC(runnumber, cent, v[0], v[1], v[2], 0, 0, 0, 0, isSelected, 0, 0);
      counter++;
      return;
    }

    const auto& foundBC = collision.foundBC_as<BCsRun3>();

    if (!foundBC.has_zdc()) {
      isSelected = false;
      spTableZDC(runnumber, cent, v[0], v[1], v[2], 0, 0, 0, 0, isSelected, 0, 0);
      counter++;
      return;
    }

    v[0] = collision.posX();
    v[1] = collision.posY();
    v[2] = collision.posZ();
    centrality = cent;
    runnumber = foundBC.runNumber();

    // load new calibrations for new runs only
    // UPLOAD Energy calibration and vmean in 1 histogram!
    if (runnumber != lastRunNumber) {
      cal.calibfilesLoaded[2] = false;
      cal.calibList[2] = nullptr;
      lastRunNumber = runnumber;
    }

    const auto& zdcCol = foundBC.zdc();

    // Get the raw energies eZN[8] (not the common A,C)
    int nTowers = 8;
    int nTowersPerSide = 4;

    for (int tower = 0; tower < nTowers; tower++) {
      eZN[tower] = (tower < nTowersPerSide) ? zdcCol.energySectorZNA()[tower] : zdcCol.energySectorZNC()[tower % nTowersPerSide];
    }

    // load the calibration histos for iteration 0 step 0 (Energy Calibration)
    loadCalibrations<kEnergyCal>(foundBC.timestamp(), cfgEnergyCal.value);

    if (!cal.calibfilesLoaded[0]) {
      if (counter < 1) {
        LOGF(info, " --> No Energy calibration files found.. -> Only Energy calibration will be done. ");
      }
    }
    // load the calibrations for the mean v
    loadCalibrations<kMeanv>(foundBC.timestamp(), cfgMeanv.value);

    registry.get<TProfile>(HIST("vmean/hvertex_vx"))->Fill(Form("%d", runnumber), v[0]);
    registry.get<TProfile>(HIST("vmean/hvertex_vy"))->Fill(Form("%d", runnumber), v[1]);
    registry.get<TProfile>(HIST("vmean/hvertex_vz"))->Fill(Form("%d", runnumber), v[2]);

    bool isZNAhit = true;
    bool isZNChit = true;

    for (int i = 0; i < nTowers; ++i) {
      if (i < nTowersPerSide && eZN[i] <= 0)
        isZNAhit = false;
      if (i >= nTowersPerSide && eZN[i] <= 0)
        isZNChit = false;
    }

    if (zdcCol.energyCommonZNA() <= 0)
      isZNAhit = false;
    if (zdcCol.energyCommonZNC() <= 0)
      isZNChit = false;

    // Fill to get mean energy per tower in 1% centrality bins
    if (isZNAhit) {
      registry.get<TProfile2D>(HIST("Energy/hZNA_mean_t0_cent"))->Fill(Form("%d", runnumber), cent, zdcCol.energyCommonZNA(), 1);
      registry.get<TProfile2D>(HIST("Energy/hZNA_mean_t1_cent"))->Fill(Form("%d", runnumber), cent, eZN[0], 1);
      registry.get<TProfile2D>(HIST("Energy/hZNA_mean_t2_cent"))->Fill(Form("%d", runnumber), cent, eZN[1], 1);
      registry.get<TProfile2D>(HIST("Energy/hZNA_mean_t3_cent"))->Fill(Form("%d", runnumber), cent, eZN[2], 1);
      registry.get<TProfile2D>(HIST("Energy/hZNA_mean_t4_cent"))->Fill(Form("%d", runnumber), cent, eZN[3], 1);
    }
    if (isZNChit) {
      registry.get<TProfile2D>(HIST("Energy/hZNC_mean_t0_cent"))->Fill(Form("%d", runnumber), cent, zdcCol.energyCommonZNC(), 1);
      registry.get<TProfile2D>(HIST("Energy/hZNC_mean_t1_cent"))->Fill(Form("%d", runnumber), cent, eZN[4], 1);
      registry.get<TProfile2D>(HIST("Energy/hZNC_mean_t2_cent"))->Fill(Form("%d", runnumber), cent, eZN[5], 1);
      registry.get<TProfile2D>(HIST("Energy/hZNC_mean_t3_cent"))->Fill(Form("%d", runnumber), cent, eZN[6], 1);
      registry.get<TProfile2D>(HIST("Energy/hZNC_mean_t4_cent"))->Fill(Form("%d", runnumber), cent, eZN[7], 1);
    }

    // if ZNA or ZNC not hit correctly.. do not use event in q-vector calculation
    if (!isZNAhit || !isZNChit) {
      counter++;
      isSelected = false;
      spTableZDC(runnumber, centrality, v[0], v[1], v[2], 0, 0, 0, 0, isSelected, 0, 0);
      return;
    }

    registry.fill(HIST("hEventCount"), evSel_isSelectedZDC);

    // Do not continue if Energy calibration is not loaded
    if (!cal.calibfilesLoaded[0]) {
      counter++;
      isSelected = false;
      spTableZDC(runnumber, centrality, v[0], v[1], v[2], 0, 0, 0, 0, isSelected, 0, 0);
      return;
    }

    // Now start gain equalisation!
    // Fill the list with calibration constants.
    for (int tower = 0; tower < (nTowers + 2); tower++) {
      meanEZN[tower] = getCorrection<TProfile2D, kEnergyCal>(namesEcal[tower].Data());
    }

    // Use the calibration constants but now only loop over towers 1-4
    int calibtower = 0;
    std::vector<int> towersNocom = {1, 2, 3, 4, 6, 7, 8, 9};

    for (const auto& tower : towersNocom) {
      if (meanEZN[tower] > 0) {
        double ecommon = (tower > nTowersPerSide) ? meanEZN[5] : meanEZN[0];
        e[calibtower] = eZN[calibtower] * (0.25 * ecommon) / meanEZN[tower];
      }
      calibtower++;
    }

    for (int i = 0; i < nTowersPerSide; i++) {
      float bincenter = i + .5;
      registry.fill(HIST("QA/ZNA_Energy"), bincenter, eZN[i]);
      registry.fill(HIST("QA/ZNA_Energy"), bincenter + 4, e[i]);
      registry.fill(HIST("QA/ZNC_Energy"), bincenter, eZN[i + 4]);
      registry.fill(HIST("QA/ZNC_Energy"), bincenter + 4, e[i + 4]);

      registry.get<TProfile>(HIST("QA/before/ZNA_pmC"))->Fill(Form("%d", runnumber), meanEZN[0]);
      registry.get<TProfile>(HIST("QA/before/ZNA_pm1"))->Fill(Form("%d", runnumber), eZN[0]);
      registry.get<TProfile>(HIST("QA/before/ZNA_pm2"))->Fill(Form("%d", runnumber), eZN[1]);
      registry.get<TProfile>(HIST("QA/before/ZNA_pm3"))->Fill(Form("%d", runnumber), eZN[2]);
      registry.get<TProfile>(HIST("QA/before/ZNA_pm4"))->Fill(Form("%d", runnumber), eZN[3]);

      registry.get<TProfile>(HIST("QA/before/ZNC_pmC"))->Fill(Form("%d", runnumber), meanEZN[5]);
      registry.get<TProfile>(HIST("QA/before/ZNC_pm1"))->Fill(Form("%d", runnumber), eZN[4]);
      registry.get<TProfile>(HIST("QA/before/ZNC_pm2"))->Fill(Form("%d", runnumber), eZN[5]);
      registry.get<TProfile>(HIST("QA/before/ZNC_pm3"))->Fill(Form("%d", runnumber), eZN[6]);
      registry.get<TProfile>(HIST("QA/before/ZNC_pm4"))->Fill(Form("%d", runnumber), eZN[7]);

      registry.get<TProfile>(HIST("QA/after/ZNA_pm1"))->Fill(Form("%d", runnumber), e[0]);
      registry.get<TProfile>(HIST("QA/after/ZNA_pm2"))->Fill(Form("%d", runnumber), e[1]);
      registry.get<TProfile>(HIST("QA/after/ZNA_pm3"))->Fill(Form("%d", runnumber), e[2]);
      registry.get<TProfile>(HIST("QA/after/ZNA_pm4"))->Fill(Form("%d", runnumber), e[3]);
      registry.get<TProfile>(HIST("QA/after/ZNC_pm1"))->Fill(Form("%d", runnumber), e[4]);
      registry.get<TProfile>(HIST("QA/after/ZNC_pm2"))->Fill(Form("%d", runnumber), e[5]);
      registry.get<TProfile>(HIST("QA/after/ZNC_pm3"))->Fill(Form("%d", runnumber), e[6]);
      registry.get<TProfile>(HIST("QA/after/ZNC_pm4"))->Fill(Form("%d", runnumber), e[7]);

      double sumZNAbefore = eZN[0] + eZN[1] + eZN[2] + eZN[3];
      double sumZNAafter = e[0] + e[1] + e[2] + e[3];

      double sumZNCbefore = eZN[4] + eZN[5] + eZN[6] + eZN[7];
      double sumZNCafter = e[4] + e[5] + e[6] + e[7];

      registry.fill(HIST("QA/") + HIST("before") + HIST("/ZNA_pmC_vs_Centrality"), centrality, zdcCol.energyCommonZNA());
      registry.fill(HIST("QA/") + HIST("before") + HIST("/ZNA_pmSUM_vs_Centrality"), centrality, sumZNAbefore);
      registry.fill(HIST("QA/") + HIST("before") + HIST("/ZNA_pm1_vs_Centrality"), centrality, eZN[0] / sumZNAbefore);
      registry.fill(HIST("QA/") + HIST("before") + HIST("/ZNA_pm2_vs_Centrality"), centrality, eZN[1] / sumZNAbefore);
      registry.fill(HIST("QA/") + HIST("before") + HIST("/ZNA_pm3_vs_Centrality"), centrality, eZN[2] / sumZNAbefore);
      registry.fill(HIST("QA/") + HIST("before") + HIST("/ZNA_pm4_vs_Centrality"), centrality, eZN[3] / sumZNAbefore);

      registry.fill(HIST("QA/") + HIST("before") + HIST("/ZNC_pmC_vs_Centrality"), centrality, zdcCol.energyCommonZNC());
      registry.fill(HIST("QA/") + HIST("before") + HIST("/ZNC_pmSUM_vs_Centrality"), centrality, sumZNCbefore);
      registry.fill(HIST("QA/") + HIST("before") + HIST("/ZNC_pm1_vs_Centrality"), centrality, eZN[4] / sumZNCbefore);
      registry.fill(HIST("QA/") + HIST("before") + HIST("/ZNC_pm2_vs_Centrality"), centrality, eZN[5] / sumZNCbefore);
      registry.fill(HIST("QA/") + HIST("before") + HIST("/ZNC_pm3_vs_Centrality"), centrality, eZN[6] / sumZNCbefore);
      registry.fill(HIST("QA/") + HIST("before") + HIST("/ZNC_pm4_vs_Centrality"), centrality, eZN[7] / sumZNCbefore);

      registry.fill(HIST("QA/") + HIST("after") + HIST("/ZNA_pmC_vs_Centrality"), centrality, zdcCol.energyCommonZNA());
      registry.fill(HIST("QA/") + HIST("after") + HIST("/ZNA_pmSUM_vs_Centrality"), centrality, sumZNAafter);
      registry.fill(HIST("QA/") + HIST("after") + HIST("/ZNA_pm1_vs_Centrality"), centrality, e[0] / sumZNAafter);
      registry.fill(HIST("QA/") + HIST("after") + HIST("/ZNA_pm2_vs_Centrality"), centrality, e[1] / sumZNAafter);
      registry.fill(HIST("QA/") + HIST("after") + HIST("/ZNA_pm3_vs_Centrality"), centrality, e[2] / sumZNAafter);
      registry.fill(HIST("QA/") + HIST("after") + HIST("/ZNA_pm4_vs_Centrality"), centrality, e[3] / sumZNAafter);

      registry.fill(HIST("QA/") + HIST("after") + HIST("/ZNC_pmC_vs_Centrality"), centrality, zdcCol.energyCommonZNC());
      registry.fill(HIST("QA/") + HIST("after") + HIST("/ZNC_pmSUM_vs_Centrality"), centrality, sumZNCafter);
      registry.fill(HIST("QA/") + HIST("after") + HIST("/ZNC_pm1_vs_Centrality"), centrality, e[4] / sumZNCafter);
      registry.fill(HIST("QA/") + HIST("after") + HIST("/ZNC_pm2_vs_Centrality"), centrality, e[5] / sumZNCafter);
      registry.fill(HIST("QA/") + HIST("after") + HIST("/ZNC_pm3_vs_Centrality"), centrality, e[6] / sumZNCafter);
      registry.fill(HIST("QA/") + HIST("after") + HIST("/ZNC_pm4_vs_Centrality"), centrality, e[7] / sumZNCafter);
    }

    // Now calculate Q-vector
    for (int tower = 0; tower < nTowers; tower++) {
      int side = (tower >= nTowersPerSide) ? 1 : 0;
      int sector = tower % nTowersPerSide;
      double energy = std::pow(e[tower], alphaZDC);
      sumZN[side] += energy;
      xEnZN[side] += (side == 0) ? -1.0 * pxZDC[sector] * energy : pxZDC[sector] * energy;
      yEnZN[side] += pyZDC[sector] * energy;

      // Also calculate the Q-vector for the non-equalized energy
      double energyNoEq = std::pow(eZN[tower], alphaZDC);
      sumZN_noEq[side] += energyNoEq;
      xEnZN_noEq[side] += (side == 0) ? -1.0 * pxZDC[sector] * energyNoEq : pxZDC[sector] * energyNoEq;
      yEnZN_noEq[side] += pyZDC[sector] * energyNoEq;
    }

    // "QXA", "QYA", "QXC", "QYC"
    int sides = 2;
    for (int i = 0; i < sides; ++i) {
      if (sumZN[i] > 0) {
        q[i * 2] = xEnZN[i] / sumZN[i];     // for QXA[0] and QXC[2]
        q[i * 2 + 1] = yEnZN[i] / sumZN[i]; // for QYA[1] and QYC[3]
      }
      if (sumZN_noEq[i] > 0) {
        qNoEq[i * 2] = xEnZN_noEq[i] / sumZN_noEq[i];     // for QXA[0] and QXC[2]
        qNoEq[i * 2 + 1] = yEnZN_noEq[i] / sumZN_noEq[i]; // for QYA[1] and QYC[3]
      }
    }

    if (cal.calibfilesLoaded[1]) {
      v[0] = v[0] - getCorrection<TProfile, kMeanv>(vnames[0].Data());
      v[1] = v[1] - getCorrection<TProfile, kMeanv>(vnames[1].Data());
    } else {
      LOGF(warning, " --> No mean V found.. -> THis wil lead to wrong axis for vx, vy (will be created in vmean/)");
      return;
    }

    loadCalibrations<kRec>(foundBC.timestamp(), cfgRec.value);

    std::vector<double> qRec(q);

    registry.get<TProfile>(HIST("QA/before/ZNA_Qx"))->Fill(Form("%d", runnumber), q[0]);
    registry.get<TProfile>(HIST("QA/before/ZNA_Qy"))->Fill(Form("%d", runnumber), q[1]);
    registry.get<TProfile>(HIST("QA/before/ZNC_Qx"))->Fill(Form("%d", runnumber), q[2]);
    registry.get<TProfile>(HIST("QA/before/ZNC_Qy"))->Fill(Form("%d", runnumber), q[3]);

    registry.get<TProfile>(HIST("QA/before/ZNA_Qx_noEq"))->Fill(Form("%d", runnumber), qNoEq[0]);
    registry.get<TProfile>(HIST("QA/before/ZNA_Qy_noEq"))->Fill(Form("%d", runnumber), qNoEq[1]);
    registry.get<TProfile>(HIST("QA/before/ZNC_Qx_noEq"))->Fill(Form("%d", runnumber), qNoEq[2]);
    registry.get<TProfile>(HIST("QA/before/ZNC_Qy_noEq"))->Fill(Form("%d", runnumber), qNoEq[3]);

    if (cal.atIteration == 0) {
      if (isSelected && cfgFillCommonRegistry)
        fillCommonRegistry<kBefore>(q[0], q[1], q[2], q[3], v, centrality);

      spTableZDC(runnumber, centrality, v[0], v[1], v[2], q[0], q[1], q[2], q[3], isSelected, 0, 0);
      counter++;
      return;
    } else {
      if (cfgFillCommonRegistry)
        fillCommonRegistry<kBefore>(q[0], q[1], q[2], q[3], v, centrality);

      // vector of 4
      std::vector<double> corrQxA;
      std::vector<double> corrQyA;
      std::vector<double> corrQxC;
      std::vector<double> corrQyC;

      int pb = 0;

      int nIterations = 5;
      int nSteps = 5;

      for (int it = 1; it <= nIterations; it++) {
        corrQxA.push_back(getCorrection<THnSparse, kRec>(names[0][0].Data(), it, 1));
        corrQyA.push_back(getCorrection<THnSparse, kRec>(names[0][1].Data(), it, 1));
        corrQxC.push_back(getCorrection<THnSparse, kRec>(names[0][2].Data(), it, 1));
        corrQyC.push_back(getCorrection<THnSparse, kRec>(names[0][3].Data(), it, 1));
        pb++;

        for (int step = 2; step <= nSteps; step++) {
          corrQxA.push_back(getCorrection<TProfile, kRec>(names[step - 1][0].Data(), it, step));
          corrQyA.push_back(getCorrection<TProfile, kRec>(names[step - 1][1].Data(), it, step));
          corrQxC.push_back(getCorrection<TProfile, kRec>(names[step - 1][2].Data(), it, step));
          corrQyC.push_back(getCorrection<TProfile, kRec>(names[step - 1][3].Data(), it, step));
          pb++;
        }
      }

      for (int cor = 0; cor < pb; cor++) {
        qRec[0] -= corrQxA[cor];
        qRec[1] -= corrQyA[cor];
        qRec[2] -= corrQxC[cor];
        qRec[3] -= corrQyC[cor];
      }

      if (isSelected && cfgFillCommonRegistry) {
        fillCommonRegistry<kAfter>(qRec[0], qRec[1], qRec[2], qRec[3], v, centrality);
        registry.fill(HIST("QA/centrality_after"), centrality);
        registry.get<TProfile>(HIST("QA/after/ZNA_Qx"))->Fill(Form("%d", runnumber), qRec[0]);
        registry.get<TProfile>(HIST("QA/after/ZNA_Qy"))->Fill(Form("%d", runnumber), qRec[1]);
        registry.get<TProfile>(HIST("QA/after/ZNC_Qx"))->Fill(Form("%d", runnumber), qRec[2]);
        registry.get<TProfile>(HIST("QA/after/ZNC_Qy"))->Fill(Form("%d", runnumber), qRec[3]);
      }

      spTableZDC(runnumber, centrality, v[0], v[1], v[2], qRec[0], qRec[1], qRec[2], qRec[3], isSelected, cal.atIteration, cal.atStep);

      qRec.clear();

      counter++;
      return;
    }
    LOGF(warning, "We return without saving table... -> THis is a problem");
  } // end of process
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ZdcQVectors>(cfgc)};
}
