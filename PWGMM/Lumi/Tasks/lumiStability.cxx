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
/// \file lumiStability.cxx
/// \brief Analysis over BCs to study the luminosity stability along time.
///
/// \author Josue Martinez Garcia, josuem@cern.ch

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonDataFormat/BunchFilling.h"
#include "DataFormatsFDD/Digit.h"
#include "DataFormatsFT0/Digit.h"
#include "DataFormatsFV0/Digit.h"
#include "DataFormatsParameters/GRPECSObject.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include <map>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;

using BCsWithTimestamps = soa::Join<aod::BCs, aod::Timestamps>;
// using CollisionWithFDD = soa::Join<aod::FDDs, aod::Collision>;

struct LumiStabilityTask {
  // Histogram registry: an object to hold your histograms
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Declare configurables
  Configurable<int> myMaxDeltaBCFDD{"myMaxDeltaBCFDD", 5, {"My BC cut"}};
  Configurable<int> myMaxDeltaBCFT0{"myMaxDeltaBCFT0", 5, {"My BC cut"}};
  Configurable<int> myMaxDeltaBCFV0{"myMaxDeltaBCFV0", 5, {"My BC cut"}};
  Configurable<int> nOrbitsConf{"nOrbitsConf", 972'288'000, "number of orbits"};
  Configurable<int> nOrbitsPerTF{"nOrbitsPerTF", 128, "number of orbits per time frame"};
  Configurable<double> minOrbitConf{"minOrbitConf", 0, "minimum orbit"};
  Configurable<bool> is2022Data{"is2022Data", true, "To 2022 data"};
  Configurable<int> minEmpty{"minEmpty", 5, "number of BCs empty for leading BC"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  parameters::GRPLHCIFData* grplhcif = nullptr;
  int nBCsPerOrbit = 3564;
  int lastRunNumber = -1;
  int64_t currentTFid = -1;
  int nOrbits = nOrbitsConf;
  double minOrbit = minOrbitConf;
  int minTimeFDD = 30;
  int64_t bcSOR = 0; // global bc of the start of the first orbit, setting 0 by default for unanchored MC
  int64_t tsSOR;
  int64_t tsEOR;
  int64_t nBCsPerTF = nOrbitsPerTF * nBCsPerOrbit; // duration of TF in bcs, should be 128*3564 or 32*3564, setting 128 orbits by default sfor unanchored MC
  std::bitset<o2::constants::lhc::LHCMaxBunches> beamPatternA;
  std::bitset<o2::constants::lhc::LHCMaxBunches> beamPatternC;
  std::bitset<o2::constants::lhc::LHCMaxBunches> bcPatternA;
  std::bitset<o2::constants::lhc::LHCMaxBunches> bcPatternC;
  std::bitset<o2::constants::lhc::LHCMaxBunches> bcPatternB;
  std::bitset<o2::constants::lhc::LHCMaxBunches> bcPatternE;

  void init(InitContext const&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    const AxisSpec axisCounts{6, -0.5, 5.5};
    const AxisSpec axisV0Counts{5, -0.5, 4.5};
    const AxisSpec axisTrigger{nBCsPerOrbit, -0.5f, nBCsPerOrbit - 0.5f};
    const AxisSpec axisPos{1000, -1, 1};
    const AxisSpec axisPosZ{1000, -25, 25};
    const AxisSpec axisNumContrib{1001, -0.5, 1000};
    const AxisSpec axisCollisionTime{1000, -50, 50};
    const AxisSpec axisTime{1000, -10, 40};
    const AxisSpec axisTimeFDD{1000, -20, 100};
    const AxisSpec axisCountsTime{2, -0.5, 1.5};
    const AxisSpec axisOrbits{static_cast<int>(nOrbits / nOrbitsPerTF), 0., static_cast<double>(nOrbits), ""};
    const AxisSpec axisTimeRate{static_cast<int>(static_cast<double>(43200) / (nOrbitsPerTF * 89e-6)), 0., 43200, ""}; // t in seconds. Histo for 12 hrs. Each bin contain one time frame (128/32 orbits for 2022/2023).
    const AxisSpec timeAxis{1200, 0., 1200., "#bf{t-t_{SOF} (min)}"};

    histos.add("hBcA", "BC pattern A; BC ; It is present", kTH1F, {axisTrigger});
    histos.add("hBcC", "BC pattern C; BC ; It is present", kTH1F, {axisTrigger});
    histos.add("hBcB", "BC pattern B; BC ; It is present", kTH1F, {axisTrigger});
    histos.add("hBcBL", "BC pattern B - Leading BC; BC ; It is present", kTH1F, {axisTrigger});
    histos.add("hBcE", "BC pattern Empty; BC ; It is present", kTH1F, {axisTrigger});
    histos.add("hvertexX", "Pos X vertex trigger; Pos x; Count ", kTH1F, {axisPos});
    histos.add("hvertexXvsTime", "Pos X vertex vs Collision Time; vertex X (cm) ; time (ns)", {HistType::kTH2F, {{axisPos}, {axisCollisionTime}}});
    histos.add("hvertexY", "Pos Y vertex trigger; Pos y; Count ", kTH1F, {axisPos});
    histos.add("hvertexZ", "Pos Z vertex trigger; Pos z; Count ", kTH1F, {axisPosZ});
    histos.add("hnumContrib", "Num of contributors; Num of contributors; Count ", kTH1I, {axisNumContrib});
    histos.add("hcollisinTime", "Collision Time; ns; Count ", kTH1F, {axisCollisionTime});
    histos.add("hOrbitFDDVertexCoinc", "", kTH1F, {axisOrbits});
    histos.add("hOrbitFDDVertex", "", kTH1F, {axisOrbits});
    histos.add("hOrbitFT0vertex", "", kTH1F, {axisOrbits});
    histos.add("hOrbitFV0Central", "", kTH1F, {axisOrbits});
    histos.add("tsValues", "", kTH1D, {{2, -0.5, 1.5}});
    histos.add("TFsPerMinute", "TFs seen in this minute (to account for failed jobs);#bf{t-t_{SOF} (min)};#bf{#it{N}_{TFs}}", kTH1F, {timeAxis});

    // time 32.766 is dummy time

    // histo about triggers
    histos.add("FDD/hCounts", "0 FDDCount - 1 FDDVertexCount - 2 FDDPPVertexCount - 3 FDDCoincidencesVertexCount - 4 FDDPPCoincidencesVertexCount - 5 FDDPPBotSidesCount; Number; counts", kTH1F, {axisCounts});
    histos.add("FDD/nBCsVsTime", "Time of TVX triggered BCs since the start of fill. FDD;;#bf{#it{N}_{BC}}", kTH1F, {timeAxis});
    histos.add("FDD/nBCsVsTimeLeadingBC", "Time of TVX triggered BCs since the start of fill. FDD;;#bf{#it{N}_{BC}}", kTH1F, {timeAxis});
    histos.add("FDD/bcVertexTriggerCTP", "vertex trigger per BC (FDD);BC in FDD; counts", kTH1F, {axisTrigger});
    histos.add("FDD/bcVertexTrigger", "vertex trigger per BC (FDD);BC in FDD; counts", kTH1F, {axisTrigger});
    histos.add("FDD/bcVertexTriggerPP", "vertex trigger per BC (FDD);BC in FDD; counts", kTH1F, {axisTrigger});
    histos.add("FDD/bcVertexTriggerCoincidence", "vertex trigger per BC (FDD) with coincidences;BC in FDD; counts", kTH1F, {axisTrigger});
    histos.add("FDD/bcVertexTriggerCoincidencePP", "vertex trigger per BC (FDD) with coincidences and Past Protection;BC in FDD; counts", kTH1F, {axisTrigger});
    histos.add("FDD/bcVertexTriggerBothSidesCoincidencePP", "vertex per BC (FDD) with coincidences, at least one side trigger and Past Protection;BC in FDD; counts", kTH1F, {axisTrigger});
    histos.add("FDD/bcSCentralTrigger", "scentral trigger per BC (FDD);BC in FDD; counts", kTH1F, {axisTrigger});
    histos.add("FDD/bcSCentralTriggerCoincidence", "scentral trigger per BC (FDD) with coincidences;BC in FDD; counts", kTH1F, {axisTrigger});
    histos.add("FDD/bcVSCTrigger", "vertex and scentral trigger per BC (FDD);BC in FDD; counts", kTH1F, {axisTrigger});
    histos.add("FDD/bcVSCTriggerCoincidence", "vertex and scentral trigger per BC (FDD) with coincidences;BC in FDD; counts", kTH1F, {axisTrigger});
    histos.add("FDD/bcCentralTrigger", "central trigger per BC (FDD);BC in FDD; counts", kTH1F, {axisTrigger});
    histos.add("FDD/bcCentralTriggerCoincidence", "central trigger per BC (FDD) with coincidences;BC in FDD; counts", kTH1F, {axisTrigger});
    histos.add("FDD/bcVCTrigger", "vertex and central trigger per BC (FDD);BC in FDD; counts", kTH1F, {axisTrigger});
    histos.add("FDD/bcVCTriggerCoincidence", "vertex and central trigger per BC (FDD) with coincidences;BC in FDD; counts", kTH1F, {axisTrigger});
    histos.add("FDD/hBcAVertex", "BC pattern A in FDD; BC in FDD ; It is present", kTH1F, {axisTrigger});
    histos.add("FDD/hBcCVertex", "BC pattern C in FDD; BC in FDD ; It is present", kTH1F, {axisTrigger});
    histos.add("FDD/hBcBVertex", "BC pattern B in FDD; BC in FDD ; It is present", kTH1F, {axisTrigger});
    histos.add("FDD/hBcBVertexL", "BC pattern B in FDD - Leading BC; BC in FDD ; It is present", kTH1F, {axisTrigger});
    histos.add("FDD/hBcEVertex", "BC pattern Empty in FDD; BC in FDD ; It is present", kTH1F, {axisTrigger});
    histos.add("FDD/timeACbcBVertex", "time bcB ; A (ns); C (ns)", {HistType::kTH2F, {{300, -15, 15}, {300, -15, 15}}});
    histos.add("FDD/timeACbcAVertex", "time bcA ; A (ns); C (ns)", {HistType::kTH2F, {{300, -15, 15}, {300, -15, 15}}});
    histos.add("FDD/timeACbcCVertex", "time bcC ; A (ns); C (ns)", {HistType::kTH2F, {{300, -15, 15}, {300, -15, 15}}});
    histos.add("FDD/timeACbcEVertex", "time bcE ; A (ns); C (ns)", {HistType::kTH2F, {{300, -15, 15}, {300, -15, 15}}});
    histos.add("FDD/hBcA", "BC pattern A in FDD; BC in FDD ; It is present", kTH1F, {axisTrigger});
    histos.add("FDD/hBcC", "BC pattern C in FDD; BC in FDD ; It is present", kTH1F, {axisTrigger});
    histos.add("FDD/hBcB", "BC pattern B in FDD; BC in FDD ; It is present", kTH1F, {axisTrigger});
    histos.add("FDD/hBcBL", "BC pattern B in FDD - Leading BC; BC in FDD ; It is present", kTH1F, {axisTrigger});
    histos.add("FDD/hBcE", "BC pattern Empty in FDD; BC in FDD ; It is present", kTH1F, {axisTrigger});
    histos.add("FDD/timeACbcB", "time bcB ; A (ns); C (ns)", {HistType::kTH2F, {{300, -15, 15}, {300, -15, 15}}});
    histos.add("FDD/timeACbcA", "time bcA ; A (ns); C (ns)", {HistType::kTH2F, {{300, -15, 15}, {300, -15, 15}}});
    histos.add("FDD/timeACbcC", "time bcC ; A (ns); C (ns)", {HistType::kTH2F, {{300, -15, 15}, {300, -15, 15}}});
    histos.add("FDD/timeACbcE", "time bcE ; A (ns); C (ns)", {HistType::kTH2F, {{300, -15, 15}, {300, -15, 15}}});
    histos.add("FDD/hTimeAVertex", "FDD time A; ns ; count", kTH1F, {axisTimeFDD});
    histos.add("FDD/hTimeCVertex", "FDD time C; ns ; count", kTH1F, {axisTimeFDD});
    histos.add("FDD/hTimeACoinc", "FDD time A; ns ; count", kTH1F, {axisTimeFDD});
    histos.add("FDD/hTimeCCoinc", "FDD time C; ns ; count", kTH1F, {axisTimeFDD});
    histos.add("FDD/hCountsTimeA2022", "0 Dummy Time - 1 Valid Time ; Kind of Time; counts", kTH1F, {axisCounts});
    histos.add("FDD/hCountsTimeC2022", "0 Dummy Time - 1 Valid Time ; Kind of Time; counts", kTH1F, {axisCounts});
    histos.add("FDD/hCountsTime2022", "0 Dummy Time - 1 Valid Time ; Kind of Time; counts", kTH1F, {axisCounts});
    histos.add("FDD/hValidTimeAvsBC2022", "Valid Time A vs BC id;BC in FT0;valid time counts", kTH1F, {axisTrigger});
    histos.add("FDD/hInvTimeAvsBC2022", "Invalid Time A vs BC id;BC in FT0;invalid time counts", kTH1F, {axisTrigger});
    histos.add("FDD/hValidTimeCvsBC2022", "Valid Time C vs BC id;BC in FT0;valid time counts", kTH1F, {axisTrigger});
    histos.add("FDD/hInvTimeCvsBC2022", "Invalid Time C vs BC id;BC in FT0;invalid time counts", kTH1F, {axisTrigger});
    histos.add("FDD/hValidTimevsBC2022", "Valid Time vs BC id;BC in FT0;valid time counts", kTH1F, {axisTrigger});
    histos.add("FDD/hInvTimevsBC2022", "Invalid Time vs BC id;BC in FT0;invalid time counts", kTH1F, {axisTrigger});
    histos.add("FDD/hCountsTimeA", "0 Dummy Time - 1 Valid Time ; Kind of Time; counts", kTH1F, {axisCounts});
    histos.add("FDD/hCountsTimeC", "0 Dummy Time - 1 Valid Time ; Kind of Time; counts", kTH1F, {axisCounts});
    histos.add("FDD/hCountsTime", "0 Dummy Time - 1 Valid Time ; Kind of Time; counts", kTH1F, {axisCounts});
    histos.add("FDD/hValidTimeAvsBC", "Valid Time A vs BC id;BC in FT0;valid time counts", kTH1F, {axisTrigger});
    histos.add("FDD/hInvTimeAvsBC", "Invalid Time A vs BC id;BC in FT0;invalid time counts", kTH1F, {axisTrigger});
    histos.add("FDD/hValidTimeCvsBC", "Valid Time C vs BC id;BC in FT0;valid time counts", kTH1F, {axisTrigger});
    histos.add("FDD/hInvTimeCvsBC", "Invalid Time C vs BC id;BC in FT0;invalid time counts", kTH1F, {axisTrigger});
    histos.add("FDD/hValidTimevsBC", "Valid Time vs BC id;BC in FT0;valid time counts", kTH1F, {axisTrigger});
    histos.add("FDD/hInvTimevsBC", "Invalid Time vs BC id;BC in FT0;invalid time counts", kTH1F, {axisTrigger});
    histos.add("FDD/hTimeForRate", "Counts by time in FDD;t (in seconds) in FDD; counts", kTH1F, {axisTimeRate});
    histos.add("FDD/hTimeForRateCTP", "Counts by time in FDD;t (in seconds) in FDD; counts", kTH1F, {axisTimeRate});
    histos.add("FDD/hTimeForRateLeadingBC", "Counts by time in FDD;t (in seconds) in FDD; counts", kTH1F, {axisTimeRate});
    histos.add("FDD/hTimeForRateLeadingBCCTP", "Counts by time in FDD;t (in seconds) in FDD; counts", kTH1F, {axisTimeRate});

    histos.add("FT0/hCounts", "0 FT0Count - 1 FT0VertexCount - 2 FT0PPVertexCount - 3 FT0PPBothSidesCount; Number; counts", kTH1F, {axisCounts});
    histos.add("FT0/nBCsVsTime", "Time of TVX triggered BCs since the start of fill. FT0;;#bf{#it{N}_{BC}}", kTH1F, {timeAxis});
    histos.add("FT0/nBCsVsTimeLeadingBC", "Time of TVX triggered BCs since the start of fill. FT0;;#bf{#it{N}_{BC}}", kTH1F, {timeAxis});
    histos.add("FT0/bcVertexTriggerCTP", "vertex trigger per BC (FT0);BC in FT0; counts", kTH1F, {axisTrigger});
    histos.add("FT0/bcVertexTrigger", "vertex trigger per BC (FT0);BC in FT0; counts", kTH1F, {axisTrigger});
    histos.add("FT0/bcVertexTriggerPP", "vertex trigger per BC (FT0) with Past Protection;BC in FT0; counts", kTH1F, {axisTrigger});
    histos.add("FT0/bcVertexTriggerBothSidesPP", "vertex per BC (FDD) with coincidences, at least one side trigger and Past Protection;BC in FDD; counts", kTH1F, {axisTrigger});
    histos.add("FT0/bcSCentralTrigger", "Scentral trigger per BC (FT0);BC in FT0; counts", kTH1F, {axisTrigger});
    histos.add("FT0/bcVSCTrigger", "vertex and Scentral trigger per BC (FT0);BC in FT0; counts", kTH1F, {axisTrigger});
    histos.add("FT0/bcCentralTrigger", "central trigger per BC (FT0);BC in FT0; counts", kTH1F, {axisTrigger});
    histos.add("FT0/bcVCTrigger", "vertex and central trigger per BC (FT0);BC in FT0; counts", kTH1F, {axisTrigger});
    histos.add("FT0/hBcA", "BC pattern A in FT0; BC in FT0 ; It is present", kTH1F, {axisTrigger});
    histos.add("FT0/hBcC", "BC pattern C in FT0; BC in FT0 ; It is present", kTH1F, {axisTrigger});
    histos.add("FT0/hBcB", "BC pattern B in FT0; BC in FT0 ; It is present", kTH1F, {axisTrigger});
    histos.add("FT0/hBcBL", "BC pattern B in FT0 - Leading BC; BC in FT0 ; It is present", kTH1F, {axisTrigger});
    histos.add("FT0/hBcE", "BC pattern Empty in FT0; BC in FT0 ; It is present", kTH1F, {axisTrigger});
    histos.add("FT0/timeACbcB", "time bcB ; A (ns); C (ns)", {HistType::kTH2F, {{300, -15, 15}, {300, -15, 15}}});
    histos.add("FT0/timeACbcA", "time bcA ; A (ns); C (ns)", {HistType::kTH2F, {{300, -15, 15}, {300, -15, 15}}});
    histos.add("FT0/timeACbcC", "time bcC ; A (ns); C (ns)", {HistType::kTH2F, {{300, -15, 15}, {300, -15, 15}}});
    histos.add("FT0/timeACbcE", "time bcE ; A (ns); C (ns)", {HistType::kTH2F, {{300, -15, 15}, {300, -15, 15}}});
    histos.add("FT0/hTimeA", "FT0 time A; ns ; count", kTH1F, {axisTime});
    histos.add("FT0/hTimeC", "FT0 time C; ns ; count", kTH1F, {axisTime});
    histos.add("FT0/hCountsTimeA", "0 Dummy Time - 1 Valid Time ; Kind of Time; counts", kTH1F, {axisCounts});
    histos.add("FT0/hCountsTimeC", "0 Dummy Time - 1 Valid Time ; Kind of Time; counts", kTH1F, {axisCounts});
    histos.add("FT0/hCountsTime", "0 Dummy Time - 1 Valid Time ; Kind of Time; counts", kTH1F, {axisCounts});
    histos.add("FT0/hValidTimeAvsBC", "Valid Time A vs BC id;BC in FT0;valid time counts", kTH1F, {axisTrigger});
    histos.add("FT0/hInvTimeAvsBC", "Invalid Time A vs BC id;BC in FT0;invalid time counts", kTH1F, {axisTrigger});
    histos.add("FT0/hValidTimeCvsBC", "Valid Time C vs BC id;BC in FT0;valid time counts", kTH1F, {axisTrigger});
    histos.add("FT0/hInvTimeCvsBC", "Invalid Time C vs BC id;BC in FT0;invalid time counts", kTH1F, {axisTrigger});
    histos.add("FT0/hValidTimevsBC", "Valid Time vs BC id;BC in FT0;valid time counts", kTH1F, {axisTrigger});
    histos.add("FT0/hInvTimevsBC", "Invalid Time vs BC id;BC in FT0;invalid time counts", kTH1F, {axisTrigger});
    histos.add("FT0/hTimeForRate", "Counts by time in FT0;t (in seconds) in FT0; counts", kTH1F, {axisTimeRate});
    histos.add("FT0/hTimeForRateCTP", "Counts by time in FT0;t (in seconds) in FT0; counts", kTH1F, {axisTimeRate});
    histos.add("FT0/hTimeForRateLeadingBC", "Counts by time in FT0;t (in seconds) in FT0; counts", kTH1F, {axisTimeRate});
    histos.add("FT0/hTimeForRateLeadingBCCTP", "Counts by time in FT0;t (in seconds) in FT0; counts", kTH1F, {axisTimeRate});

    histos.add("FV0/hCounts", "0 CountCentralFV0 - 1 CountPFPCentralFV0 - 2 CountPFPOutInFV0 - 3 CountPPCentralFV0 - 4 CountPPOutInFV0; Number; counts", kTH1F, {axisV0Counts});
    histos.add("FV0/bcChargeTriggerCTP", "Out trigger per BC (FV0);BC in V0; counts", kTH1F, {axisTrigger});
    histos.add("FV0/bcOutTrigger", "Out trigger per BC (FV0);BC in V0; counts", kTH1F, {axisTrigger});
    histos.add("FV0/bcInTrigger", "In trigger per BC (FV0);BC in V0; counts", kTH1F, {axisTrigger});
    histos.add("FV0/bcSCenTrigger", "SCen trigger per BC (FV0);BC in V0; counts", kTH1F, {axisTrigger});
    histos.add("FV0/bcCenTrigger", "Central trigger per BC (FV0);BC in V0; counts", kTH1F, {axisTrigger});
    histos.add("FV0/bcCenTriggerPFPCentral", "Central trigger per BC (FV0) with PFP in central trigger;BC in V0; counts", kTH1F, {axisTrigger});
    histos.add("FV0/bcCenTriggerPPCentral", "Central trigger per BC (FV0) with PP in central trigger;BC in V0; counts", kTH1F, {axisTrigger});
    histos.add("FV0/bcCenTriggerPFPOutIn", "Central trigger per BC (FV0) with PFP in Out and In trigger;BC in V0; counts", kTH1F, {axisTrigger});
    histos.add("FV0/bcCenTriggerPPOutIn", "Central trigger per BC (FV0) with PP in Out and In trigger;BC in V0; counts", kTH1F, {axisTrigger});
    histos.add("FV0/hBcA", "BC pattern A in FV0; BC in FV0 ; It is present", kTH1F, {axisTrigger});
    histos.add("FV0/hBcC", "BC pattern C in FV0; BC in FV0 ; It is present", kTH1F, {axisTrigger});
    histos.add("FV0/hBcB", "BC pattern B in FV0; BC in FV0 ; It is present", kTH1F, {axisTrigger});
    histos.add("FV0/hBcE", "BC pattern Empty in FV0; BC in FV0 ; It is present", kTH1F, {axisTrigger});
    histos.add("FV0/timeAbcB", "time bcB ; A (ns)", kTH1F, {{300, -15, 15}});
    histos.add("FV0/timeAbcA", "time bcA ; A (ns)", kTH1F, {{300, -15, 15}});
    histos.add("FV0/timeAbcC", "time bcC ; A (ns)", kTH1F, {{300, -15, 15}});
    histos.add("FV0/timeAbcE", "time bcE ; A (ns)", kTH1F, {{300, -15, 15}});
    histos.add("FV0/hTimeForRateCTP", "Counts by time in FV0;t (in seconds) in FV0; counts", kTH1F, {axisTimeRate});
    histos.add("FV0/hTimeForRateLeadingBCCTP", "Counts by time in FV0;t (in seconds) in FV0; counts", kTH1F, {axisTimeRate});
  }

  bool checkAnyCoincidence(const std::vector<int>& channels)
  {
    constexpr std::pair<int, int> kPair0 = {0, 4};
    constexpr std::pair<int, int> kPair1 = {1, 5};
    constexpr std::pair<int, int> kPair2 = {2, 6};
    constexpr std::pair<int, int> kPair3 = {3, 7};
    constexpr std::array<std::pair<int, int>, 4> kChannelPairs = {kPair0, kPair1, kPair2, kPair3};
    // std::map<int, int> kChannelPairs = {{0, 4}, {1, 5}, {2, 6}, {3, 7}};
    for (const auto& pair : kChannelPairs) {
      if (std::find(channels.begin(), channels.end(), pair.first) != channels.end() &&
          std::find(channels.begin(), channels.end(), pair.second) != channels.end()) {
        return true;
      }
    }
    return false;
  }

  float getTimeSinceSOF(const auto& bc)
  {
    return (bc.timestamp() - grplhcif->getFillNumberTime()) / 1e3 / 60; // Convert to minutes
  }

  void processMain(aod::FDDs const& fdds, aod::FT0s const& ft0s, aod::FV0As const& fv0s, aod::BCsWithTimestamps const& bcs)
  {
    int executionCounter = 0;
    int nbin = o2::constants::lhc::LHCMaxBunches;
    uint32_t nOrbitsPerTF = 0; // 128 in 2022, 32 in 2023
    if (is2022Data) {
      nOrbitsPerTF = 128; // 128 in 2022, 32 in 2023
    } else {
      nOrbitsPerTF = 32; // 128 in 2022, 32 in 2023
    }
    int runNumber = bcs.iteratorAt(0).runNumber();
    // std::string histName = "hOrbitFDDVertexCoinc_" + std::to_string(runNumber);
    if (runNumber != lastRunNumber && executionCounter < 1) {
      tsSOR = 0;
      tsEOR = 1;
      lastRunNumber = runNumber; // do it only once
      executionCounter++;

      // access CCDB for data or anchored MC only
      int64_t ts = bcs.iteratorAt(0).timestamp();

      // access colliding and beam-gas bc patterns
      grplhcif = ccdb->getForTimeStamp<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", ts);
      beamPatternA = grplhcif->getBunchFilling().getBeamPattern(0);
      beamPatternC = grplhcif->getBunchFilling().getBeamPattern(1);
      bcPatternA = beamPatternA & ~beamPatternC;
      bcPatternC = ~beamPatternA & beamPatternC;
      bcPatternB = beamPatternA & beamPatternC;
      bcPatternE = ~beamPatternA & ~beamPatternC;

      for (int i = 0; i < nBCsPerOrbit; i++) {
        if (bcPatternA[i]) {
          histos.fill(HIST("hBcA"), i);
        }
        if (bcPatternC[i]) {
          histos.fill(HIST("hBcC"), i);
        }
        if (bcPatternB[i]) {
          histos.fill(HIST("hBcB"), i);
          bool isLeadBC = true;
          for (int jbit = i - minEmpty; jbit < i; jbit++) {
            int kbit = jbit;
            if (kbit < 0)
              kbit += nbin;
            if (bcPatternB[kbit]) {
              isLeadBC = false;
              break;
            }
          }
          if (isLeadBC)
            histos.fill(HIST("hBcBL"), i);
        }
        if (bcPatternE[i]) {
          histos.fill(HIST("hBcE"), i);
        }
      }

      EventSelectionParams* par = ccdb->getForTimeStamp<EventSelectionParams>("EventSelection/EventSelectionParams", ts);
      // access orbit-reset timestamp
      auto ctpx = ccdb->getForTimeStamp<std::vector<int64_t>>("CTP/Calib/OrbitReset", ts);
      int64_t tsOrbitReset = (*ctpx)[0]; // us
      // access TF duration, start-of-run and end-of-run timestamps from ECS GRP
      std::map<std::string, std::string> metadata;
      metadata["runNumber"] = Form("%d", runNumber);
      auto grpecs = ccdb->getSpecific<o2::parameters::GRPECSObject>("GLO/Config/GRPECS", ts, metadata);
      nOrbitsPerTF = grpecs->getNHBFPerTF(); // assuming 1 orbit = 1 HBF;  nOrbitsPerTF=128 in 2022, 32 in 2023
      tsSOR = grpecs->getTimeStart();        // ms
      tsEOR = grpecs->getTimeEnd();          // ms
      // calculate SOR and EOR orbits
      int64_t orbitSOR = (tsSOR * 1000 - tsOrbitReset) / o2::constants::lhc::LHCOrbitMUS;
      int64_t orbitEOR = (tsEOR * 1000 - tsOrbitReset) / o2::constants::lhc::LHCOrbitMUS;
      // adjust to the nearest TF edge
      orbitSOR = orbitSOR / nOrbitsPerTF * nOrbitsPerTF + par->fTimeFrameOrbitShift;
      orbitEOR = orbitEOR / nOrbitsPerTF * nOrbitsPerTF + par->fTimeFrameOrbitShift;
      // set nOrbits and minOrbit used for orbit-axis binning
      nOrbits = orbitEOR - orbitSOR;
      minOrbit = orbitSOR;
      // first bc of the first orbit (should coincide with TF start)
      bcSOR = orbitSOR * o2::constants::lhc::LHCMaxBunches;
      // duration of TF in bcs
      nBCsPerTF = nOrbitsPerTF * o2::constants::lhc::LHCMaxBunches;
      LOGP(info, "tsOrbitReset={} us, SOR = {} ms, EOR = {} ms, orbitSOR = {}, nBCsPerTF = {}", tsOrbitReset, tsSOR, tsEOR, orbitSOR, nBCsPerTF);

      auto hTsValues = histos.get<TH1>(HIST("tsValues"));
      hTsValues->GetXaxis()->SetBinLabel(1, "tsSOR");
      hTsValues->GetXaxis()->SetBinLabel(2, "tsEOR");
      hTsValues->SetBinContent(1, tsSOR / 1000); // seconds
      hTsValues->SetBinContent(2, tsEOR / 1000); // seconds

      // create orbit-axis histograms on the fly with binning based on info from GRP if GRP is available
      // otherwise default minOrbit and nOrbits will be used
      // const AxisSpec axisOrbits{static_cast<int>(nOrbits / nOrbitsPerTF), 0., static_cast<double>(nOrbits), ""};
      // histos.add("hOrbitFDDVertexCoinc", "", kTH1F, {axisOrbits});
      // histos.add("hOrbitFDDVertex", "", kTH1F, {axisOrbits});
      // histos.add("hOrbitFT0vertex", "", kTH1F, {axisOrbits});
      // histos.add("hOrbitFV0Central", "", kTH1F, {axisOrbits});
    }

    for (auto const& bc : bcs) {
      if (bc.timestamp() == 0) {
        continue;
      }

      std::bitset<64> ctpInputMask(bc.inputMask());
      bool trgFDD = ctpInputMask[15];
      bool trgFT0 = ctpInputMask[2];
      bool trgFV0 = ctpInputMask[9];

      int64_t globalBC = bc.globalBC();
      int localBC = globalBC % nBCsPerOrbit;

      float timeSinceSOF = getTimeSinceSOF(bc);

      int64_t thisTFid = (globalBC - bcSOR) / nBCsPerTF;

      if (thisTFid != currentTFid) {
        currentTFid = thisTFid;
        histos.fill(HIST("TFsPerMinute"), timeSinceSOF);
      }

      if (trgFDD) {
        histos.fill(HIST("FDD/bcVertexTriggerCTP"), localBC + 7);
        if (bcPatternB[localBC]) {
          histos.fill(HIST("FDD/nBCsVsTime"), timeSinceSOF);
          histos.fill(HIST("FDD/hTimeForRateCTP"), (bc.timestamp() - tsSOR) * 1.e-3); // Converting ms into seconds
        }
      }
      if (trgFT0) {
        histos.fill(HIST("FT0/bcVertexTriggerCTP"), localBC);
        if (bcPatternB[localBC]) {
          histos.fill(HIST("FT0/nBCsVsTime"), timeSinceSOF);
          histos.fill(HIST("FT0/hTimeForRateCTP"), (bc.timestamp() - tsSOR) * 1.e-3); // Converting ms into seconds
        }
      }
      if (trgFV0) {
        histos.fill(HIST("FV0/bcChargeTriggerCTP"), localBC);
        if (bcPatternB[localBC]) {
          histos.fill(HIST("FV0/hTimeForRateCTP"), (bc.timestamp() - tsSOR) * 1.e-3); // Converting ms into seconds
        }
      }
      bool isLeadBC = true;
      for (int jbit = localBC - minEmpty; jbit < localBC; jbit++) {
        int kbit = jbit;
        if (kbit < 0)
          kbit += nbin;
        if (bcPatternB[kbit]) {
          isLeadBC = false;
          break;
        }
      }
      if (isLeadBC) {
        if (trgFDD) {
          histos.fill(HIST("FDD/nBCsVsTimeLeadingBC"), timeSinceSOF);
          histos.fill(HIST("FDD/hTimeForRateLeadingBCCTP"), (bc.timestamp() - tsSOR) * 1.e-3);
        }
        if (trgFT0) {
          histos.fill(HIST("FT0/nBCsVsTimeLeadingBC"), timeSinceSOF);
          histos.fill(HIST("FT0/hTimeForRateLeadingBCCTP"), (bc.timestamp() - tsSOR) * 1.e-3);
        }
        if (trgFV0) {
          histos.fill(HIST("FV0/hTimeForRateLeadingBCCTP"), (bc.timestamp() - tsSOR) * 1.e-3);
        }
      }
      // }
    } // loop over bcs

    for (auto const& fdd : fdds) {
      auto bc = fdd.bc_as<BCsWithTimestamps>();
      if (bc.timestamp() == 0) {
        continue;
      }

      int64_t globalBC = bc.globalBC();
      int localBC = globalBC % nBCsPerOrbit;
      uint64_t orbit = globalBC / nBCsPerOrbit;

      std::bitset<8> fddTriggers = fdd.triggerMask();
      bool vertex = fddTriggers[o2::fdd::Triggers::bitVertex];
      bool scentral = fddTriggers[o2::fdd::Triggers::bitSCen];
      bool central = fddTriggers[o2::fdd::Triggers::bitCen];

      auto sideA = fdd.chargeA();
      auto sideC = fdd.chargeC();
      std::vector<int> channelA;
      std::vector<int> channelC;
      int minLimit = 0;
      int maxNChanels = 8;
      for (auto i = 0; i < maxNChanels; i++) {
        if (sideA[i] > minLimit) {
          channelA.push_back(i);
        }
        if (sideC[i] > minLimit) {
          channelC.push_back(i);
        }
      }

      bool isCoinA = checkAnyCoincidence(channelA);
      bool isCoinC = checkAnyCoincidence(channelC);

      histos.fill(HIST("FDD/hCounts"), 0);
      if (vertex) {
        histos.fill(HIST("FDD/bcVertexTrigger"), localBC);
        histos.fill(HIST("FDD/hCounts"), 1);
        histos.fill(HIST("hOrbitFDDVertex"), orbit - minOrbit);

        if (bcPatternB[localBC]) {
          histos.fill(HIST("FDD/hTimeForRate"), (bc.timestamp() - tsSOR) * 1.e-3); // Converting ms into seconds
          bool isLeadBC = true;
          for (int jbit = localBC - minEmpty; jbit < localBC; jbit++) {
            int kbit = jbit;
            if (kbit < 0)
              kbit += nbin;
            if (bcPatternB[kbit]) {
              isLeadBC = false;
              break;
            }
          }
          if (isLeadBC)
            histos.fill(HIST("FDD/hTimeForRateLeadingBC"), (bc.timestamp() - tsSOR) * 1.e-3);
        }

        int deltaIndex = 0; // backward move counts
        int deltaBC = 0;    // current difference wrt globalBC
        bool pastActivityFDDVertex = false;
        while (deltaBC < myMaxDeltaBCFDD) {
          deltaIndex++;
          if (fdd.globalIndex() - deltaIndex < 0) {
            break;
          }
          const auto& fddPast = fdds.iteratorAt(fdd.globalIndex() - deltaIndex);
          auto bcPast = fddPast.bc_as<BCsWithTimestamps>();
          deltaBC = globalBC - bcPast.globalBC();

          if (deltaBC < myMaxDeltaBCFDD) {
            std::bitset<8> fddTriggersPast = fddPast.triggerMask();
            bool vertexPast = fddTriggersPast[o2::fdd::Triggers::bitVertex];
            pastActivityFDDVertex |= (vertexPast);
          }
        }
        deltaIndex = 0;
        deltaBC = 0;

        if (pastActivityFDDVertex == false) {
          histos.fill(HIST("FDD/hCounts"), 2);
          histos.fill(HIST("FDD/bcVertexTriggerPP"), localBC);
          if (bcPatternA[localBC]) {
            histos.fill(HIST("FDD/timeACbcAVertex"), fdd.timeA(), fdd.timeC());
            histos.fill(HIST("FDD/hBcAVertex"), localBC);
          }
          if (bcPatternC[localBC]) {
            histos.fill(HIST("FDD/timeACbcCVertex"), fdd.timeA(), fdd.timeC());
            histos.fill(HIST("FDD/hBcCVertex"), localBC);
          }
          if (bcPatternB[localBC]) {
            histos.fill(HIST("FDD/timeACbcBVertex"), fdd.timeA(), fdd.timeC());
            histos.fill(HIST("FDD/hBcBVertex"), localBC);
            bool isLeadBC = true;
            for (int jbit = localBC - minEmpty; jbit < localBC; jbit++) {
              int kbit = jbit;
              if (kbit < 0)
                kbit += nbin;
              if (bcPatternB[kbit]) {
                isLeadBC = false;
                break;
              }
            }
            if (isLeadBC)
              histos.fill(HIST("FDD/hBcBVertexL"), localBC);
            histos.fill(HIST("FDD/hTimeAVertex"), fdd.timeA());
            histos.fill(HIST("FDD/hTimeCVertex"), fdd.timeC());
            if (is2022Data) {
              if (fdd.timeA() > minTimeFDD) {
                histos.fill(HIST("FDD/hCountsTimeA2022"), 0);
                histos.fill(HIST("FDD/hInvTimeAvsBC2022"), localBC);
              } else {
                histos.fill(HIST("FDD/hCountsTimeA2022"), 1);
                histos.fill(HIST("FDD/hValidTimeAvsBC2022"), localBC);
              }

              if (fdd.timeC() > minTimeFDD) {
                histos.fill(HIST("FDD/hCountsTimeC2022"), 0);
                histos.fill(HIST("FDD/hInvTimeCvsBC2022"), localBC);
              } else {
                histos.fill(HIST("FDD/hCountsTimeC2022"), 1);
                histos.fill(HIST("FDD/hValidTimeCvsBC2022"), localBC);
              }

              if (fdd.timeA() > minTimeFDD || fdd.timeC() > minTimeFDD) {
                histos.fill(HIST("FDD/hCountsTime2022"), 0);
                histos.fill(HIST("FDD/hInvTimevsBC2022"), localBC);
              }
              if (fdd.timeA() < minTimeFDD && fdd.timeC() < minTimeFDD) {
                histos.fill(HIST("FDD/hCountsTime2022"), 1);
                histos.fill(HIST("FDD/hValidTimevsBC2022"), localBC);
              }
            }
          }
          if (bcPatternE[localBC]) {
            histos.fill(HIST("FDD/timeACbcEVertex"), fdd.timeA(), fdd.timeC());
            histos.fill(HIST("FDD/hBcEVertex"), localBC);
          }
        }

          if (deltaBC < myMaxDeltaBCFDD) {
            std::bitset<8> fddTriggersPast = fdd_past.triggerMask();
            bool vertexPast = fddTriggersPast[o2::fdd::Triggers::bitVertex];
            pastActivityFDDVertex |= (vertexPast);
          }
        }
        deltaIndex = 0;
        deltaBC = 0;

        if (pastActivityFDDVertex == false) {
          histos.fill(HIST("FDD/hCounts"), 2);
          histos.fill(HIST("FDD/bcVertexTriggerPP"), localBC);
          if (bcPatternA[localBC]) {
            histos.fill(HIST("FDD/timeACbcAVertex"), fdd.timeA(), fdd.timeC());
            histos.fill(HIST("FDD/hBcAVertex"), localBC);
          }
          if (bcPatternC[localBC]) {
            histos.fill(HIST("FDD/timeACbcC"), fdd.timeA(), fdd.timeC());
            histos.fill(HIST("FDD/hBcC"), localBC);
          }
          if (bcPatternB[localBC]) {
            histos.fill(HIST("FDD/timeACbcB"), fdd.timeA(), fdd.timeC());
            histos.fill(HIST("FDD/hBcB"), localBC);
            bool isLeadBC = true;
            for (int jbit = localBC - minEmpty; jbit < localBC; jbit++) {
              int kbit = jbit;
              if (kbit < 0)
                kbit += nbin;
              if (bcPatternB[kbit]) {
                isLeadBC = false;
                break;
              }
            }
            if (isLeadBC)
              histos.fill(HIST("FDD/hBcBL"), localBC);
            histos.fill(HIST("FDD/hTimeACoinc"), fdd.timeA());
            histos.fill(HIST("FDD/hTimeCCoinc"), fdd.timeC());
            if (!is2022Data) {
              if (fdd.timeA() > minTimeFDD) {
                histos.fill(HIST("FDD/hCountsTimeA"), 0);
                histos.fill(HIST("FDD/hInvTimeAvsBC"), localBC);
              } else {
                histos.fill(HIST("FDD/hCountsTimeA"), 1);
                histos.fill(HIST("FDD/hValidTimeAvsBC"), localBC);
              }

              if (fdd.timeC() > minTimeFDD) {
                histos.fill(HIST("FDD/hCountsTimeC"), 0);
                histos.fill(HIST("FDD/hInvTimeCvsBC"), localBC);
              } else {
                histos.fill(HIST("FDD/hCountsTimeC"), 1);
                histos.fill(HIST("FDD/hValidTimeCvsBC"), localBC);
              }

              if (fdd.timeA() > minTimeFDD || fdd.timeC() > minTimeFDD) {
                histos.fill(HIST("FDD/hCountsTime"), 0);
                histos.fill(HIST("FDD/hInvTimevsBC"), localBC);
              }
              if (fdd.timeA() < minTimeFDD && fdd.timeC() < minTimeFDD) {
                histos.fill(HIST("FDD/hCountsTime"), 1);
                histos.fill(HIST("FDD/hValidTimevsBC"), localBC);
              }
            }
          }
        }


        if (isCoinA && isCoinC) {
          histos.fill(HIST("FDD/bcVertexTriggerCoincidence"), localBC);
          histos.fill(HIST("FDD/hCounts"), 3);
          histos.fill(HIST("hOrbitFDDVertexCoinc"), orbit - minOrbit);

          if (bcPatternA[localBC]) {
            histos.fill(HIST("FDD/timeACbcA"), fdd.timeA(), fdd.timeC());
            histos.fill(HIST("FDD/hBcA"), localBC);
          }
          if (bcPatternC[localBC]) {
            histos.fill(HIST("FDD/timeACbcC"), fdd.timeA(), fdd.timeC());
            histos.fill(HIST("FDD/hBcC"), localBC);
          }
          if (bcPatternB[localBC]) {
            histos.fill(HIST("FDD/timeACbcB"), fdd.timeA(), fdd.timeC());
            histos.fill(HIST("FDD/hBcB"), localBC);
            histos.fill(HIST("FDD/hTimeACoinc"), fdd.timeA());
            histos.fill(HIST("FDD/hTimeCCoinc"), fdd.timeC());
          }
          if (bcPatternE[localBC]) {
            histos.fill(HIST("FDD/timeACbcE"), fdd.timeA(), fdd.timeC());
            histos.fill(HIST("FDD/hBcE"), localBC);
          }

          int deltaIndex = 0; // backward move counts
          int deltaBC = 0;    // current difference wrt globalBC
          bool pastActivityFDDVertexCoincidences = false;
          bool pastActivityFDDTriggerACoincidenceA = false;
          bool pastActivityFDDTriggerCCoincidenceC = false;
          while (deltaBC < myMaxDeltaBCFDD) {
            deltaIndex++;
            if (fdd.globalIndex() - deltaIndex < 0) {
              break;
            }
            const auto& fddPast = fdds.iteratorAt(fdd.globalIndex() - deltaIndex);
            auto bcPast = fddPast.bc_as<BCsWithTimestamps>();
            deltaBC = globalBC - bcPast.globalBC();

            if (deltaBC < myMaxDeltaBCFDD) {
              std::bitset<8> fddTriggersPast = fddPast.triggerMask();
              bool vertexPast = fddTriggersPast[o2::fdd::Triggers::bitVertex];
              bool triggerAPast = fddTriggersPast[o2::fdd::Triggers::bitA];
              bool triggerCPast = fddTriggersPast[o2::fdd::Triggers::bitC];
              auto sideAPast = fddPast.chargeA();
              auto sideCPast = fddPast.chargeC();
              std::vector<int> channelAPast;
              std::vector<int> channelCPast;
              int maxNChanels = 8;
              for (auto i = 0; i < maxNChanels; i++) {
                if (sideAPast[i] > 0) {
                  channelAPast.push_back(i);
                }
                if (sideCPast[i] > 0) {
                  channelCPast.push_back(i);
                }
              }

              bool isCoinAPast = checkAnyCoincidence(channelAPast);
              bool isCoinCPast = checkAnyCoincidence(channelCPast);
              pastActivityFDDVertexCoincidences |= (vertexPast & isCoinAPast & isCoinCPast);
              pastActivityFDDTriggerACoincidenceA |= (triggerAPast & isCoinAPast);
              pastActivityFDDTriggerCCoincidenceC |= (triggerCPast & isCoinCPast);
            }
          }
          deltaIndex = 0;
          deltaBC = 0;

          if (pastActivityFDDVertexCoincidences == false) {
            histos.fill(HIST("FDD/hCounts"), 4);
            histos.fill(HIST("FDD/bcVertexTriggerCoincidencePP"), localBC);
          }
          if (pastActivityFDDTriggerACoincidenceA == false || pastActivityFDDTriggerCCoincidenceC == false) {
            histos.fill(HIST("FDD/hCounts"), 5);
            histos.fill(HIST("FDD/bcVertexTriggerBothSidesCoincidencePP"), localBC);
          }
        } // coincidences
      } // vertex true

      if (scentral) {
        histos.fill(HIST("FDD/bcSCentralTrigger"), localBC);
        if (isCoinA && isCoinC) {
          histos.fill(HIST("FDD/bcSCentralTriggerCoincidence"), localBC);
        }
      } // central true

      if (vertex && scentral) {
        histos.fill(HIST("FDD/bcVSCTrigger"), localBC);
        if (isCoinA && isCoinC) {
          histos.fill(HIST("FDD/bcVSCTriggerCoincidence"), localBC);
        }
      } // vertex and scentral true

      if (central) {
        histos.fill(HIST("FDD/bcCentralTrigger"), localBC);
        if (isCoinA && isCoinC) {
          histos.fill(HIST("FDD/bcCentralTriggerCoincidence"), localBC);
        }
      }

      if (vertex && central) {
        histos.fill(HIST("FDD/bcVCTrigger"), localBC);
        if (isCoinA && isCoinC) {
          histos.fill(HIST("FDD/bcVCTriggerCoincidence"), localBC);
        }
      } // vertex and scentral true
    } // loop over FDD events

    for (auto const& ft0 : ft0s) {
      auto bc = ft0.bc_as<BCsWithTimestamps>();
      if (bc.timestamp() == 0) {
        continue;
      }

      int64_t globalBC = bc.globalBC();
      int localBC = globalBC % nBCsPerOrbit;
      uint64_t orbit = globalBC / nBCsPerOrbit;

      std::bitset<8> fT0Triggers = ft0.triggerMask();
      bool vertex = fT0Triggers[o2::ft0::Triggers::bitVertex];
      bool sCentral = fT0Triggers[o2::ft0::Triggers::bitSCen];
      bool central = fT0Triggers[o2::ft0::Triggers::bitCen];

      histos.fill(HIST("FT0/hCounts"), 0);
      if (vertex) {
        histos.fill(HIST("FT0/bcVertexTrigger"), localBC);
        histos.fill(HIST("hOrbitFT0vertex"), orbit - minOrbit);

        if (bcPatternA[localBC]) {
          histos.fill(HIST("FT0/timeACbcA"), ft0.timeA(), ft0.timeC());
          histos.fill(HIST("FT0/hBcA"), localBC);
        }
        if (bcPatternC[localBC]) {
          histos.fill(HIST("FT0/timeACbcC"), ft0.timeA(), ft0.timeC());
          histos.fill(HIST("FT0/hBcC"), localBC);
        }
        if (bcPatternB[localBC]) {
          histos.fill(HIST("FT0/timeACbcB"), ft0.timeA(), ft0.timeC());
          histos.fill(HIST("FT0/hBcB"), localBC);
          histos.fill(HIST("FT0/hTimeForRate"), (bc.timestamp() - tsSOR) * 1.e-3); // Converting ms into seconds
          bool isLeadBC = true;
          for (int jbit = localBC - minEmpty; jbit < localBC; jbit++) {
            int kbit = jbit;
            if (kbit < 0)
              kbit += nbin;
            if (bcPatternB[kbit]) {
              isLeadBC = false;
              break;
            }
          }
          if (isLeadBC) {
            histos.fill(HIST("FT0/hTimeForRateLeadingBC"), (bc.timestamp() - tsSOR) * 1.e-3); // Converting ms into seconds
            histos.fill(HIST("FT0/hBcBL"), localBC);
          }
          histos.fill(HIST("FT0/hTimeA"), ft0.timeA());
          histos.fill(HIST("FT0/hTimeC"), ft0.timeC());

          if (ft0.timeA() > minTimeFDD) {
            histos.fill(HIST("FT0/hCountsTimeA"), 0);
            histos.fill(HIST("FT0/hInvTimeAvsBC"), localBC);
          } else {
            histos.fill(HIST("FT0/hCountsTimeA"), 1);
            histos.fill(HIST("FT0/hValidTimeAvsBC"), localBC);
          }

          if (ft0.timeC() > minTimeFDD) {
            histos.fill(HIST("FT0/hCountsTimeC"), 0);
            histos.fill(HIST("FT0/hInvTimeCvsBC"), localBC);
          } else {
            histos.fill(HIST("FT0/hCountsTimeC"), 1);
            histos.fill(HIST("FT0/hValidTimeCvsBC"), localBC);
          }

          if (ft0.timeA() > minTimeFDD || ft0.timeC() > minTimeFDD) {
            histos.fill(HIST("FT0/hCountsTime"), 0);
            histos.fill(HIST("FT0/hInvTimevsBC"), localBC);
          }
          if (ft0.timeA() < minTimeFDD && ft0.timeC() < minTimeFDD) {
            histos.fill(HIST("FT0/hCountsTime"), 1);
            histos.fill(HIST("FT0/hValidTimevsBC"), localBC);
          }
        }
        if (bcPatternE[localBC]) {
          histos.fill(HIST("FT0/timeACbcE"), ft0.timeA(), ft0.timeC());
          histos.fill(HIST("FT0/hBcE"), localBC);
        }

        int deltaIndex = 0; // backward move counts
        int deltaBC = 0;    // current difference wrt globalBC
        bool pastActivityFT0Vertex = false;
        bool pastActivityFT0TriggerA = false;
        bool pastActivityFT0TriggerC = false;
        while (deltaBC < myMaxDeltaBCFT0) {
          deltaIndex++;
          if (ft0.globalIndex() - deltaIndex < 0) {
            break;
          }
          const auto& ft0Past = ft0s.iteratorAt(ft0.globalIndex() - deltaIndex);
          auto bcPast = ft0Past.bc_as<BCsWithTimestamps>();
          deltaBC = globalBC - bcPast.globalBC();

          if (deltaBC < myMaxDeltaBCFT0) {
            std::bitset<8> fT0TriggersPast = ft0Past.triggerMask();
            bool vertexPast = fT0TriggersPast[o2::ft0::Triggers::bitVertex];
            bool triggerAPast = fT0TriggersPast[o2::ft0::Triggers::bitA];
            bool triggerCPast = fT0TriggersPast[o2::ft0::Triggers::bitC];

            pastActivityFT0Vertex |= vertexPast;
            pastActivityFT0TriggerA |= triggerAPast;
            pastActivityFT0TriggerC |= triggerCPast;
          }
        }
        deltaIndex = 0;
        deltaBC = 0;

        histos.fill(HIST("FT0/hCounts"), 1);
        if (pastActivityFT0Vertex == false) {
          histos.fill(HIST("FT0/hCounts"), 2);
          histos.fill(HIST("FT0/bcVertexTriggerPP"), localBC);
        }
        if (pastActivityFT0TriggerA == false || pastActivityFT0TriggerC == false) {
          histos.fill(HIST("FT0/hCounts"), 3);
          histos.fill(HIST("FT0/bcVertexTriggerBothSidesPP"), localBC);
        }
        if (pastActivityFT0Vertex == true) {
          histos.fill(HIST("FT0/hCounts"), 3);
        } else {
          histos.fill(HIST("FT0/bcVertexTriggerPP"), localBC);
        }
        if (pastActivityFT0TriggerA == false || pastActivityFT0TriggerC == false) {
          histos.fill(HIST("FT0/hCounts"), 3);
          histos.fill(HIST("FT0/bcVertexTriggerBothSidesPP"), localBC);
        }
      } // vertex true

      if (sCentral) {
        histos.fill(HIST("FT0/bcSCentralTrigger"), localBC);
        if (vertex) {
          histos.fill(HIST("FT0/bcVSCTrigger"), localBC);
        }
      } // scentral true

      if (central) {
        histos.fill(HIST("FT0/bcCentralTrigger"), localBC);
        if (sCentral) {
          histos.fill(HIST("FT0/bcSCentralCentralTrigger"), localBC);
        }
        if (vertex) {
          histos.fill(HIST("FT0/bcVCTrigger"), localBC);
        }
      }
    } // loop over FT0 events

    for (auto const& fv0 : fv0s) {
      auto bc = fv0.bc_as<BCsWithTimestamps>();
      if (bc.timestamp() == 0) {
        continue;
      }

      int64_t globalBC = bc.globalBC();
      int localBC = globalBC % nBCsPerOrbit;
      uint64_t orbit = globalBC / nBCsPerOrbit;

      std::bitset<8> fv0Triggers = fv0.triggerMask();
      bool aOut = fv0Triggers[o2::fv0::Triggers::bitAOut];
      bool aIn = fv0Triggers[o2::fv0::Triggers::bitAIn];
      bool aSCen = fv0Triggers[o2::fv0::Triggers::bitTrgNchan];
      bool aCen = fv0Triggers[o2::fv0::Triggers::bitTrgCharge];

      if (aOut) {
        histos.fill(HIST("FV0/bcOutTrigger"), localBC);
      }

      if (aIn) {
        histos.fill(HIST("FV0/bcInTrigger"), localBC);
      }

      if (aSCen) {
        histos.fill(HIST("FV0/bcSCenTrigger"), localBC);
      }

      if (aCen) {
        histos.fill(HIST("hOrbitFV0Central"), orbit - minOrbit);
        histos.fill(HIST("FV0/bcCenTrigger"), localBC);

        if (bcPatternA[localBC]) {
          histos.fill(HIST("FV0/timeAbcA"), fv0.time());
          histos.fill(HIST("FV0/hBcA"), localBC);
        }
        if (bcPatternC[localBC]) {
          histos.fill(HIST("FV0/timeAbcC"), fv0.time());
          histos.fill(HIST("FV0/hBcC"), localBC);
        }
        if (bcPatternB[localBC]) {
          histos.fill(HIST("FV0/timeAbcB"), fv0.time());
          histos.fill(HIST("FV0/hBcB"), localBC);
        }
        if (bcPatternE[localBC]) {
          histos.fill(HIST("FV0/timeAbcE"), fv0.time());
          histos.fill(HIST("FV0/hBcE"), localBC);
        }

        int deltaIndex = 0; // backward move counts
        int deltaBC = 0;    // current difference wrt globalBC
        bool pastActivityFV0Cen = false;
        bool pastActivityFV0TriggerOut = false;
        bool pastActivityFV0TriggerIn = false;
        while (deltaBC < myMaxDeltaBCFV0) {
          deltaIndex++;
          if (fv0.globalIndex() - deltaIndex < 0) {
            break;
          }
          const auto& fv0Past = fv0s.iteratorAt(fv0.globalIndex() - deltaIndex);
          auto bcPast = fv0Past.bc_as<BCsWithTimestamps>();
          deltaBC = globalBC - bcPast.globalBC();

          if (deltaBC < myMaxDeltaBCFV0) {
            std::bitset<8> fv0Triggers = fv0Past.triggerMask();
            bool centralPast = fv0Triggers[o2::fv0::Triggers::bitTrgCharge];
            bool triggerOutPast = fv0Triggers[o2::fv0::Triggers::bitAOut];
            bool triggerInPast = fv0Triggers[o2::fv0::Triggers::bitAIn];

            pastActivityFV0Cen |= centralPast;
            pastActivityFV0TriggerOut |= triggerOutPast;
            pastActivityFV0TriggerIn |= triggerInPast;
          }
        }
        deltaIndex = 0;
        deltaBC = 0;

        bool futureActivityFV0Cen = false;
        bool futureActivityFV0TriggerOut = false;
        bool futureActivityFV0TriggerIn = false;
        while (deltaBC < myMaxDeltaBCFV0) {
          deltaIndex++;
          if (fv0.globalIndex() + deltaIndex >= fv0s.size()) {
            break;
          }
          const auto& fv0Future = fv0s.iteratorAt(fv0.globalIndex() + deltaIndex);
          deltaBC = fv0Future.bcId() - fv0.bcId();

          if (deltaBC < myMaxDeltaBCFV0) {
            std::bitset<8> fv0Triggers = fv0Future.triggerMask();
            bool centralFuture = fv0Triggers[o2::fv0::Triggers::bitTrgCharge];
            bool triggerOutFuture = fv0Triggers[o2::fv0::Triggers::bitAOut];
            bool triggerInFuture = fv0Triggers[o2::fv0::Triggers::bitAIn];

            futureActivityFV0Cen |= centralFuture;
            futureActivityFV0TriggerOut |= triggerOutFuture;
            futureActivityFV0TriggerIn |= triggerInFuture;
          }
        }

        histos.fill(HIST("FV0/hCounts"), 0);
        if ((pastActivityFV0TriggerOut || futureActivityFV0TriggerOut) == true || (pastActivityFV0TriggerIn || futureActivityFV0TriggerIn) == true) {
          histos.fill(HIST("FV0/hCounts"), 2);
        } else {
          histos.fill(HIST("FV0/bcCenTriggerPFPOutIn"), localBC);
        }
        if (pastActivityFV0TriggerOut == true || pastActivityFV0TriggerIn == true) {
          histos.fill(HIST("FV0/hCounts"), 4);
        } else {
          histos.fill(HIST("FV0/bcCenTriggerPPOutIn"), localBC);
        }
        if (pastActivityFV0Cen == true || futureActivityFV0Cen == true) {
          histos.fill(HIST("FV0/hCounts"), 1);
        } else {
          histos.fill(HIST("FV0/bcCenTriggerPFPCentral"), localBC);
        }
        if (pastActivityFV0Cen == true) {
          histos.fill(HIST("FV0/hCounts"), 3);
        } else {
          histos.fill(HIST("FV0/bcCenTriggerPPCentral"), localBC);
        }
      }
    } // loop over V0 events
  } // end processMain

  PROCESS_SWITCH(LumiStabilityTask, processMain, "Process FDD and FT0 to lumi stability analysis", true);

  void processCollisions(aod::Collision const& collision)
  {
    histos.fill(HIST("hvertexX"), collision.posX());
    histos.fill(HIST("hvertexY"), collision.posY());
    histos.fill(HIST("hvertexZ"), collision.posZ());
    histos.fill(HIST("hnumContrib"), collision.numContrib());
    histos.fill(HIST("hcollisinTime"), collision.collisionTime());
    histos.fill(HIST("hvertexXvsTime"), collision.posX(), collision.collisionTime());
  }

  PROCESS_SWITCH(LumiStabilityTask, processCollisions, "Process collision to get vertex position", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<LumiStabilityTask>(cfgc)};
}
