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
///
/// \brief This task is an empty skeleton that fills a simple eta histogram.
///        it is meant to be a blank page for further developments.
/// \author everyone

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/ASoAHelpers.h"
#include "DataFormatsFDD/Digit.h"
#include "DataFormatsFT0/Digit.h"
#include "DataFormatsFV0/Digit.h"
#include "Framework/ASoA.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "CCDB/BasicCCDBManager.h"
#include "CommonDataFormat/BunchFilling.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "DataFormatsParameters/GRPECSObject.h"

using namespace o2;
using namespace o2::framework;

using BCsWithTimestamps = soa::Join<aod::BCs, aod::Timestamps>;

struct lumiStabilityTask {
  // Histogram registry: an object to hold your histograms
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Declare configurables
  Configurable<int> myMaxDeltaBCFDD{"myMaxDeltaBCFDD", 5, {"My BC cut"}};
  Configurable<int> myMaxDeltaBCFT0{"myMaxDeltaBCFT0", 5, {"My BC cut"}};
  Configurable<int> myMaxDeltaBCFV0{"myMaxDeltaBCFV0", 5, {"My BC cut"}};
  Configurable<int> nOrbitsConf{"nOrbits", 10000, "number of orbits"};
  Configurable<double> minOrbitConf{"minOrbit", 0, "minimum orbit"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int nBCsPerOrbit = 3564;
  int lastRunNumber = -1;
  int nOrbits = nOrbitsConf;
  double minOrbit = minOrbitConf;
  int64_t bcSOR = 0;                      // global bc of the start of the first orbit, setting 0 by default for unanchored MC
  int64_t nBCsPerTF = 128 * nBCsPerOrbit; // duration of TF in bcs, should be 128*3564 or 32*3564, setting 128 orbits by default sfor unanchored MC
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
    const AxisSpec axisTriggger{nBCsPerOrbit, -0.5f, nBCsPerOrbit - 0.5f};

    histos.add("hBcA", "BC pattern A; BC ; It is present", kTH1F, {axisTriggger});
    histos.add("hBcC", "BC pattern C; BC ; It is present", kTH1F, {axisTriggger});
    histos.add("hBcB", "BC pattern B; BC ; It is present", kTH1F, {axisTriggger});
    histos.add("hBcE", "BC pattern Empty; BC ; It is present", kTH1F, {axisTriggger});

    // histo about triggers
    histos.add("FDD/hCounts", "0 FDDCount - 1 FDDVertexCount - 2 FDDPPVertexCount - 3 FDDCoincidencesVertexCount - 4 FDDPPCoincidencesVertexCount - 5 FDDPPBotSidesCount; Number; counts", kTH1F, {axisCounts});
    histos.add("FDD/bcVertexTrigger", "vertex trigger per BC (FDD);BC in FDD; counts", kTH1F, {axisTriggger});
    histos.add("FDD/bcVertexTriggerPP", "vertex trigger per BC (FDD);BC in FDD; counts", kTH1F, {axisTriggger});
    histos.add("FDD/bcVertexTriggerCoincidence", "vertex trigger per BC (FDD) with coincidences;BC in FDD; counts", kTH1F, {axisTriggger});
    histos.add("FDD/bcVertexTriggerCoincidencePP", "vertex trigger per BC (FDD) with coincidences and Past Protection;BC in FDD; counts", kTH1F, {axisTriggger});
    histos.add("FDD/bcVertexTriggerBothSidesCoincidencePP", "vertex per BC (FDD) with coincidences, at least one side trigger and Past Protection;BC in FDD; counts", kTH1F, {axisTriggger});
    histos.add("FDD/bcSCentralTrigger", "scentral trigger per BC (FDD);BC in FDD; counts", kTH1F, {axisTriggger});
    histos.add("FDD/bcSCentralTriggerCoincidence", "scentral trigger per BC (FDD) with coincidences;BC in FDD; counts", kTH1F, {axisTriggger});
    histos.add("FDD/bcVSCTrigger", "vertex and scentral trigger per BC (FDD);BC in FDD; counts", kTH1F, {axisTriggger});
    histos.add("FDD/bcVSCTriggerCoincidence", "vertex and scentral trigger per BC (FDD) with coincidences;BC in FDD; counts", kTH1F, {axisTriggger});
    histos.add("FDD/bcCentralTrigger", "central trigger per BC (FDD);BC in FDD; counts", kTH1F, {axisTriggger});
    histos.add("FDD/bcCentralTriggerCoincidence", "central trigger per BC (FDD) with coincidences;BC in FDD; counts", kTH1F, {axisTriggger});
    histos.add("FDD/bcVCTrigger", "vertex and central trigger per BC (FDD);BC in FDD; counts", kTH1F, {axisTriggger});
    histos.add("FDD/bcVCTriggerCoincidence", "vertex and central trigger per BC (FDD) with coincidences;BC in FDD; counts", kTH1F, {axisTriggger});
    histos.add("FDD/hBcAVertex", "BC pattern A in FDD; BC in FDD ; It is present", kTH1F, {axisTriggger});
    histos.add("FDD/hBcCVertex", "BC pattern C in FDD; BC in FDD ; It is present", kTH1F, {axisTriggger});
    histos.add("FDD/hBcBVertex", "BC pattern B in FDD; BC in FDD ; It is present", kTH1F, {axisTriggger});
    histos.add("FDD/hBcEVertex", "BC pattern Empty in FDD; BC in FDD ; It is present", kTH1F, {axisTriggger});
    histos.add("FDD/timeACbcBVertex", "time bcB ; A (ns); C (ns)", {HistType::kTH2F, {{300, -15, 15}, {300, -15, 15}}});
    histos.add("FDD/timeACbcAVertex", "time bcA ; A (ns); C (ns)", {HistType::kTH2F, {{300, -15, 15}, {300, -15, 15}}});
    histos.add("FDD/timeACbcCVertex", "time bcC ; A (ns); C (ns)", {HistType::kTH2F, {{300, -15, 15}, {300, -15, 15}}});
    histos.add("FDD/timeACbcEVertex", "time bcE ; A (ns); C (ns)", {HistType::kTH2F, {{300, -15, 15}, {300, -15, 15}}});
    histos.add("FDD/hBcA", "BC pattern A in FDD; BC in FDD ; It is present", kTH1F, {axisTriggger});
    histos.add("FDD/hBcC", "BC pattern C in FDD; BC in FDD ; It is present", kTH1F, {axisTriggger});
    histos.add("FDD/hBcB", "BC pattern B in FDD; BC in FDD ; It is present", kTH1F, {axisTriggger});
    histos.add("FDD/hBcE", "BC pattern Empty in FDD; BC in FDD ; It is present", kTH1F, {axisTriggger});
    histos.add("FDD/timeACbcB", "time bcB ; A (ns); C (ns)", {HistType::kTH2F, {{300, -15, 15}, {300, -15, 15}}});
    histos.add("FDD/timeACbcA", "time bcA ; A (ns); C (ns)", {HistType::kTH2F, {{300, -15, 15}, {300, -15, 15}}});
    histos.add("FDD/timeACbcC", "time bcC ; A (ns); C (ns)", {HistType::kTH2F, {{300, -15, 15}, {300, -15, 15}}});
    histos.add("FDD/timeACbcE", "time bcE ; A (ns); C (ns)", {HistType::kTH2F, {{300, -15, 15}, {300, -15, 15}}});

    histos.add("FT0/hCounts", "0 FT0Count - 1 FT0VertexCount - 2 FT0PPVertexCount - 3 FT0PPBothSidesCount; Number; counts", kTH1F, {axisCounts});
    histos.add("FT0/bcVertexTrigger", "vertex trigger per BC (FT0);BC in FT0; counts", kTH1F, {axisTriggger});
    histos.add("FT0/bcVertexTriggerPP", "vertex trigger per BC (FT0) with Past Protection;BC in FT0; counts", kTH1F, {axisTriggger});
    histos.add("FT0/bcVertexTriggerBothSidesPP", "vertex per BC (FDD) with coincidences, at least one side trigger and Past Protection;BC in FDD; counts", kTH1F, {axisTriggger});
    histos.add("FT0/bcSCentralTrigger", "Scentral trigger per BC (FT0);BC in FT0; counts", kTH1F, {axisTriggger});
    histos.add("FT0/bcVSCTrigger", "vertex and Scentral trigger per BC (FT0);BC in FT0; counts", kTH1F, {axisTriggger});
    histos.add("FT0/bcCentralTrigger", "central trigger per BC (FT0);BC in FT0; counts", kTH1F, {axisTriggger});
    histos.add("FT0/bcVCTrigger", "vertex and central trigger per BC (FT0);BC in FT0; counts", kTH1F, {axisTriggger});
    histos.add("FT0/hBcA", "BC pattern A in FT0; BC in FT0 ; It is present", kTH1F, {axisTriggger});
    histos.add("FT0/hBcC", "BC pattern C in FT0; BC in FT0 ; It is present", kTH1F, {axisTriggger});
    histos.add("FT0/hBcB", "BC pattern B in FT0; BC in FT0 ; It is present", kTH1F, {axisTriggger});
    histos.add("FT0/hBcE", "BC pattern Empty in FT0; BC in FT0 ; It is present", kTH1F, {axisTriggger});
    histos.add("FT0/timeACbcB", "time bcB ; A (ns); C (ns)", {HistType::kTH2F, {{300, -15, 15}, {300, -15, 15}}});
    histos.add("FT0/timeACbcA", "time bcA ; A (ns); C (ns)", {HistType::kTH2F, {{300, -15, 15}, {300, -15, 15}}});
    histos.add("FT0/timeACbcC", "time bcC ; A (ns); C (ns)", {HistType::kTH2F, {{300, -15, 15}, {300, -15, 15}}});
    histos.add("FT0/timeACbcE", "time bcE ; A (ns); C (ns)", {HistType::kTH2F, {{300, -15, 15}, {300, -15, 15}}});

    histos.add("FV0/hCounts", "0 CountCentralFV0 - 1 CountPFPCentralFV0 - 2 CountPFPOutInFV0 - 3 CountPPCentralFV0 - 4 CountPPOutInFV0; Number; counts", kTH1F, {axisV0Counts});
    histos.add("FV0/bcOutTrigger", "Out trigger per BC (FV0);BC in V0; counts", kTH1F, {axisTriggger});
    histos.add("FV0/bcInTrigger", "In trigger per BC (FV0);BC in V0; counts", kTH1F, {axisTriggger});
    histos.add("FV0/bcSCenTrigger", "SCen trigger per BC (FV0);BC in V0; counts", kTH1F, {axisTriggger});
    histos.add("FV0/bcCenTrigger", "Central trigger per BC (FV0);BC in V0; counts", kTH1F, {axisTriggger});
    histos.add("FV0/bcCenTriggerPFPCentral", "Central trigger per BC (FV0) with PFP in central trigger;BC in V0; counts", kTH1F, {axisTriggger});
    histos.add("FV0/bcCenTriggerPPCentral", "Central trigger per BC (FV0) with PP in central trigger;BC in V0; counts", kTH1F, {axisTriggger});
    histos.add("FV0/bcCenTriggerPFPOutIn", "Central trigger per BC (FV0) with PFP in Out and In trigger;BC in V0; counts", kTH1F, {axisTriggger});
    histos.add("FV0/bcCenTriggerPPOutIn", "Central trigger per BC (FV0) with PP in Out and In trigger;BC in V0; counts", kTH1F, {axisTriggger});
  }

  bool checkAnyCoincidence(const std::vector<int>& channels)
  {
    constexpr std::pair<int, int> pair0 = {0, 4};
    constexpr std::pair<int, int> pair1 = {1, 5};
    constexpr std::pair<int, int> pair2 = {2, 6};
    constexpr std::pair<int, int> pair3 = {3, 7};
    constexpr std::array<std::pair<int, int>, 4> channelPairs = {pair0, pair1, pair2, pair3};
    // std::map<int, int> channelPairs = {{0, 4}, {1, 5}, {2, 6}, {3, 7}};
    for (const auto& pair : channelPairs) {
      if (std::find(channels.begin(), channels.end(), pair.first) != channels.end() &&
          std::find(channels.begin(), channels.end(), pair.second) != channels.end()) {
        return true;
      }
    }
    return false;
  }

  void processMain(aod::FDDs const& fdds, aod::FT0s const& ft0s, aod::FV0As const& fv0s, aod::BCsWithTimestamps const& bcs)
  {
    uint32_t nOrbitsPerTF = 128; // 128 in 2022, 32 in 2023
    int runNumber = bcs.iteratorAt(0).runNumber();
    if (runNumber != lastRunNumber) {
      lastRunNumber = runNumber; // do it only once
      int64_t tsSOR = 0;
      int64_t tsEOR = 1;

      if (runNumber >= 500000) { // access CCDB for data or anchored MC only
        int64_t ts = bcs.iteratorAt(0).timestamp();

        // access colliding and beam-gas bc patterns
        auto grplhcif = ccdb->getForTimeStamp<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", ts);
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
          }
          if (bcPatternE[i]) {
            histos.fill(HIST("hBcE"), i);
          }
        }

        EventSelectionParams* par = ccdb->getForTimeStamp<EventSelectionParams>("EventSelection/EventSelectionParams", ts);
        // access orbit-reset timestamp
        auto ctpx = ccdb->getForTimeStamp<std::vector<Long64_t>>("CTP/Calib/OrbitReset", ts);
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
      }

      // create orbit-axis histograms on the fly with binning based on info from GRP if GRP is available
      // otherwise default minOrbit and nOrbits will be used
      const AxisSpec axisOrbits{static_cast<int>(nOrbits / nOrbitsPerTF), 0., static_cast<double>(nOrbits), ""};
      histos.add("hOrbitFDD", "FDD Orbits; Orbit; Entries", kTH1F, {axisOrbits});
      histos.add("hOrbitFDDVertex", "FDD Orbits; Orbit; Entries", kTH1F, {axisOrbits});
      histos.add("hOrbitFT0", "FT0 Orbits; Orbit; Entries", kTH1F, {axisOrbits});
      histos.add("hOrbitFV0", "FV0 Orbits; Orbit; Entries", kTH1F, {axisOrbits});
    }

    for (auto const& fdd : fdds) {
      auto bc = fdd.bc_as<BCsWithTimestamps>();
      if (bc.timestamp() == 0) {
        continue;
      }

      Long64_t globalBC = bc.globalBC();
      int localBC = globalBC % nBCsPerOrbit;
      uint64_t orbit = globalBC / nBCsPerOrbit;

      std::bitset<8> fddTriggers = fdd.triggerMask();
      bool vertex = fddTriggers[o2::fdd::Triggers::bitVertex];
      bool scentral = fddTriggers[o2::fdd::Triggers::bitSCen];
      bool central = fddTriggers[o2::fdd::Triggers::bitCen];

      auto SideA = fdd.chargeA();
      auto SideC = fdd.chargeC();
      std::vector<int> channelA;
      std::vector<int> channelC;
      for (auto i = 0; i < 8; i++) {
        if (SideA[i] > 0) {
          channelA.push_back(i);
        }
        if (SideC[i] > 0) {
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
        }
        if (bcPatternE[localBC]) {
          histos.fill(HIST("FDD/timeACbcEVertex"), fdd.timeA(), fdd.timeC());
          histos.fill(HIST("FDD/hBcEVertex"), localBC);
        }

        int deltaIndex = 0; // backward move counts
        int deltaBC = 0;    // current difference wrt globalBC
        bool pastActivityFDDVertex = false;
        while (deltaBC < myMaxDeltaBCFDD) {
          deltaIndex++;
          if (fdd.globalIndex() - deltaIndex < 0) {
            break;
          }
          const auto& fdd_past = fdds.iteratorAt(fdd.globalIndex() - deltaIndex);
          auto bc_past = fdd_past.bc_as<BCsWithTimestamps>();
          deltaBC = globalBC - bc_past.globalBC();

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
        }

        if (isCoinA && isCoinC) {
          histos.fill(HIST("FDD/bcVertexTriggerCoincidence"), localBC);
          histos.fill(HIST("FDD/hCounts"), 1);
          histos.fill(HIST("hOrbitFDD"), orbit - minOrbit);

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
            const auto& fdd_past = fdds.iteratorAt(fdd.globalIndex() - deltaIndex);
            auto bc_past = fdd_past.bc_as<BCsWithTimestamps>();
            deltaBC = globalBC - bc_past.globalBC();

            if (deltaBC < myMaxDeltaBCFDD) {
              std::bitset<8> fddTriggersPast = fdd_past.triggerMask();
              bool vertexPast = fddTriggersPast[o2::fdd::Triggers::bitVertex];
              bool triggerAPast = fddTriggersPast[o2::fdd::Triggers::bitA];
              bool triggerCPast = fddTriggersPast[o2::fdd::Triggers::bitC];
              auto SideAPast = fdd_past.chargeA();
              auto SideCPast = fdd_past.chargeC();
              std::vector<int> channelAPast;
              std::vector<int> channelCPast;
              for (auto i = 0; i < 8; i++) {
                if (SideAPast[i] > 0) {
                  channelAPast.push_back(i);
                }
                if (SideCPast[i] > 0) {
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
            histos.fill(HIST("FDD/hCounts"), 2);
            histos.fill(HIST("FDD/bcVertexTriggerCoincidencePP"), localBC);
          }
          if (pastActivityFDDTriggerACoincidenceA == false || pastActivityFDDTriggerCCoincidenceC == false) {
            histos.fill(HIST("FDD/hCounts"), 3);
            histos.fill(HIST("FDD/bcVertexTriggerBothSidesCoincidencePP"), localBC);
          }
        } // coincidences
      }   // vertex true

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
    }   // loop over FDD events

    for (auto const& ft0 : ft0s) {
      auto bc = ft0.bc_as<BCsWithTimestamps>();
      if (bc.timestamp() == 0) {
        continue;
      }

      Long64_t globalBC = bc.globalBC();
      int localBC = globalBC % nBCsPerOrbit;
      uint64_t orbit = globalBC / nBCsPerOrbit;

      std::bitset<8> fT0Triggers = ft0.triggerMask();
      bool vertex = fT0Triggers[o2::ft0::Triggers::bitVertex];
      bool sCentral = fT0Triggers[o2::ft0::Triggers::bitSCen];
      bool central = fT0Triggers[o2::ft0::Triggers::bitCen];

      histos.fill(HIST("FT0/hCounts"), 0);
      if (vertex) {
        histos.fill(HIST("FT0/bcVertexTrigger"), localBC);
        histos.fill(HIST("hOrbitFT0"), orbit - minOrbit);

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
          const auto& ft0_past = ft0s.iteratorAt(ft0.globalIndex() - deltaIndex);
          auto bc_past = ft0_past.bc_as<BCsWithTimestamps>();
          deltaBC = globalBC - bc_past.globalBC();

          if (deltaBC < myMaxDeltaBCFT0) {
            std::bitset<8> fT0TriggersPast = ft0_past.triggerMask();
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

      Long64_t globalBC = bc.globalBC();
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
        histos.fill(HIST("hOrbitFV0"), orbit - minOrbit);
        histos.fill(HIST("FV0/bcCenTrigger"), localBC);

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
          const auto& fv0_past = fv0s.iteratorAt(fv0.globalIndex() - deltaIndex);
          auto bc_past = fv0_past.bc_as<BCsWithTimestamps>();
          deltaBC = globalBC - bc_past.globalBC();

          if (deltaBC < myMaxDeltaBCFV0) {
            std::bitset<8> fv0Triggers = fv0_past.triggerMask();
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
          const auto& fv0_future = fv0s.iteratorAt(fv0.globalIndex() + deltaIndex);
          deltaBC = fv0_future.bcId() - fv0.bcId();

          if (deltaBC < myMaxDeltaBCFV0) {
            std::bitset<8> fv0Triggers = fv0_future.triggerMask();
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
  }   // end processMain

  PROCESS_SWITCH(lumiStabilityTask, processMain, "Process FDD and FT0 to lumi stability analysis", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lumiStabilityTask>(cfgc)};
}
