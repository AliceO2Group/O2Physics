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

using namespace o2;
using namespace o2::framework;

using BCsWithTimestamps = soa::Join<aod::BCs, aod::Timestamps>;
using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;

int nBCsPerOrbit = 3564;

struct lumiStabilityTask {
  // Histogram registry: an object to hold your histograms
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Declare configurables
  Configurable<int> myMaxDeltaBCFDD{"myMaxDeltaBCFDD", 5, {"My BC cut"}};
  Configurable<int> myMaxDeltaBCFT0{"myMaxDeltaBCFT0", 5, {"My BC cut"}};
  Configurable<int> myMaxDeltaBCFV0{"myMaxDeltaBCFV0", 5, {"My BC cut"}};

  void init(InitContext const&)
  {
    const AxisSpec axisCounts{3, -0.5, 2.5};
    const AxisSpec axisTriggger{nBCsPerOrbit, -0.5f, nBCsPerOrbit - 0.5f};

    // histo about triggers
    histos.add("FDD/hCounts", "0 CountVertexFDD - 1 CountPFPVertexCoincidencesFDD - 2 CountPFPTriggerCoincidencesFDD; Number of Count; counts", kTH1F, {axisCounts});
    histos.add("FDD/bcVertexTrigger", "vertex trigger per BC (FDD);BC in FDD; counts", kTH1F, {axisTriggger});
    histos.add("FDD/bcVertexTriggerCoincidence", "vertex trigger per BC (FDD) with coincidences;BC in FDD; counts", kTH1F, {axisTriggger});
    histos.add("FDD/bcVertexTriggerCoincidencePFP", "vertex trigger per BC (FDD) with coincidences and Past Future Protection;BC in FDD; counts", kTH1F, {axisTriggger});
    histos.add("FDD/bcVertexTriggerBothSidesCoincidencePFP", "vertex per BC (FDD) with coincidences, at least one side trigger and Past Future Protection;BC in FDD; counts", kTH1F, {axisTriggger});
    histos.add("FDD/bcSCentralTrigger", "scentral trigger per BC (FDD);BC in FDD; counts", kTH1F, {axisTriggger});
    histos.add("FDD/bcSCentralTriggerCoincidence", "scentral trigger per BC (FDD) with coincidences;BC in FDD; counts", kTH1F, {axisTriggger});
    histos.add("FDD/bcVSCTrigger", "vertex and scentral trigger per BC (FDD);BC in FDD; counts", kTH1F, {axisTriggger});
    histos.add("FDD/bcVSCTriggerCoincidence", "vertex and scentral trigger per BC (FDD) with coincidences;BC in FDD; counts", kTH1F, {axisTriggger});
    histos.add("FDD/bcCentralTrigger", "central trigger per BC (FDD);BC in FDD; counts", kTH1F, {axisTriggger});
    histos.add("FDD/bcCentralTriggerCoincidence", "central trigger per BC (FDD) with coincidences;BC in FDD; counts", kTH1F, {axisTriggger});
    histos.add("FDD/bcVCTrigger", "vertex and central trigger per BC (FDD);BC in FDD; counts", kTH1F, {axisTriggger});
    histos.add("FDD/bcVCTriggerCoincidence", "vertex and central trigger per BC (FDD) with coincidences;BC in FDD; counts", kTH1F, {axisTriggger});

    histos.add("FT0/hCounts", "0 CountVertexFT0 - 1 CountPFPVertexCoincidencesFT0 - 2 CountPFPTriggerCoincidencesFT0; Number of Count; counts", kTH1F, {axisCounts});
    histos.add("FT0/bcVertexTrigger", "vertex trigger per BC (FT0);BC in FT0; counts", kTH1F, {axisTriggger});
    histos.add("FT0/bcVertexTriggerPFP", "vertex trigger per BC (FT0) with Past Future Protection;BC in FT0; counts", kTH1F, {axisTriggger});
    histos.add("FT0/bcVertexTriggerBothSidesPFP", "vertex per BC (FDD) with coincidences, at least one side trigger and Past Future Protection;BC in FDD; counts", kTH1F, {axisTriggger});
    histos.add("FT0/bcSCentralTrigger", "Scentral trigger per BC (FT0);BC in FT0; counts", kTH1F, {axisTriggger});
    histos.add("FT0/bcVSCTrigger", "vertex and Scentral trigger per BC (FT0);BC in FT0; counts", kTH1F, {axisTriggger});
    histos.add("FT0/bcCentralTrigger", "central trigger per BC (FT0);BC in FT0; counts", kTH1F, {axisTriggger});
    histos.add("FT0/bcVCTrigger", "vertex and central trigger per BC (FT0);BC in FT0; counts", kTH1F, {axisTriggger});

    histos.add("FV0/hCounts", "0 CountCentralFV0 - 1 CountPFPCentralFV0 - 2 CountPFPOutInFV0; Number of Count; counts", kTH1F, {axisCounts});
    histos.add("FV0/bcOutTrigger", "Out trigger per BC (FV0);BC in V0; counts", kTH1F, {axisTriggger});
    histos.add("FV0/bcInTrigger", "In trigger per BC (FV0);BC in V0; counts", kTH1F, {axisTriggger});
    histos.add("FV0/bcSCenTrigger", "SCen trigger per BC (FV0);BC in V0; counts", kTH1F, {axisTriggger});
    histos.add("FV0/bcCenTrigger", "Central trigger per BC (FV0);BC in V0; counts", kTH1F, {axisTriggger});
    histos.add("FV0/bcCenTriggerPFPCentral", "Central trigger per BC (FV0) with PFP in central trigger;BC in V0; counts", kTH1F, {axisTriggger});
    histos.add("FV0/bcCenTriggerPFPOutIn", "Central trigger per BC (FV0) with PFP in Out and In trigger;BC in V0; counts", kTH1F, {axisTriggger});
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

  void processMain(aod::FDDs const& fdds, aod::FT0s const& ft0s, aod::FV0As const& fv0s, aod::BCsWithTimestamps const&)
  {
    for (auto const& fdd : fdds) {
      auto bc = fdd.bc_as<BCsWithTimestamps>();
      if (bc.timestamp() == 0) {
        continue;
      }

      Long64_t globalBC = bc.globalBC();
      int localBC = globalBC % nBCsPerOrbit;

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

      if (vertex) {
        histos.fill(HIST("FDD/bcVertexTrigger"), localBC);
        if (isCoinA && isCoinC) {
          histos.fill(HIST("FDD/bcVertexTriggerCoincidence"), localBC);

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
            deltaBC = fdd.bcId() - fdd_past.bcId();

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

          bool futureActivityFDDVertexCoincidences = false;
          bool futureActivityFDDTriggerACoincidenceA = false;
          bool futureActivityFDDTriggerCCoincidenceC = false;
          while (deltaBC < myMaxDeltaBCFDD) {
            deltaIndex++;
            if (fdd.globalIndex() + deltaIndex >= fdds.size()) {
              break;
            }
            const auto& fdd_future = fdds.iteratorAt(fdd.globalIndex() + deltaIndex);
            deltaBC = fdd_future.bcId() - fdd.bcId();

            if (deltaBC < myMaxDeltaBCFDD) {
              std::bitset<8> fddTriggersFuture = fdd_future.triggerMask();
              bool vertexFuture = fddTriggersFuture[o2::fdd::Triggers::bitVertex];
              bool triggerAFuture = fddTriggersFuture[o2::fdd::Triggers::bitA];
              bool triggerCFuture = fddTriggersFuture[o2::fdd::Triggers::bitC];
              auto SideAFuture = fdd_future.chargeA();
              auto SideCFuture = fdd_future.chargeC();
              std::vector<int> channelAFuture;
              std::vector<int> channelCFuture;
              for (auto i = 0; i < 8; i++) {
                if (SideAFuture[i] > 0) {
                  channelAFuture.push_back(i);
                }
                if (SideCFuture[i] > 0) {
                  channelCFuture.push_back(i);
                }
              }

              bool isCoinAFuture = checkAnyCoincidence(channelAFuture);
              bool isCoinCFuture = checkAnyCoincidence(channelCFuture);
              futureActivityFDDVertexCoincidences |= (vertexFuture & isCoinAFuture & isCoinCFuture);
              futureActivityFDDTriggerACoincidenceA |= (triggerAFuture & isCoinAFuture);
              futureActivityFDDTriggerCCoincidenceC |= (triggerCFuture & isCoinCFuture);
            }
          }

          histos.fill(HIST("FDD/hCounts"), 0);
          if ((pastActivityFDDTriggerACoincidenceA || futureActivityFDDTriggerACoincidenceA) == true || (pastActivityFDDTriggerCCoincidenceC || futureActivityFDDTriggerCCoincidenceC) == true) {
            histos.fill(HIST("FDD/hCounts"), 2);
          } else {
            histos.fill(HIST("FDD/bcVertexTriggerBothSidesCoincidencePFP"), localBC);
          }
          if (pastActivityFDDVertexCoincidences == true || futureActivityFDDVertexCoincidences == true) {
            histos.fill(HIST("FDD/hCounts"), 1);
          } else {
            histos.fill(HIST("FDD/bcVertexTriggerCoincidencePFP"), localBC);
          }
        }
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
    }   // loop over FDD events

    for (auto const& ft0 : ft0s) {
      auto bc = ft0.bc_as<BCsWithTimestamps>();
      if (bc.timestamp() == 0) {
        continue;
      }

      Long64_t globalBC = bc.globalBC();
      int localBC = globalBC % nBCsPerOrbit;

      std::bitset<8> fT0Triggers = ft0.triggerMask();
      bool vertex = fT0Triggers[o2::ft0::Triggers::bitVertex];
      bool sCentral = fT0Triggers[o2::ft0::Triggers::bitSCen];
      bool central = fT0Triggers[o2::ft0::Triggers::bitCen];

      if (vertex) {
        histos.fill(HIST("FT0/bcVertexTrigger"), localBC);

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
          deltaBC = ft0.bcId() - ft0_past.bcId();

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

        bool futureActivityFT0Vertex = false;
        bool futureActivityFT0TriggerA = false;
        bool futureActivityFT0TriggerC = false;
        while (deltaBC < myMaxDeltaBCFT0) {
          deltaIndex++;
          if (ft0.globalIndex() + deltaIndex >= ft0s.size()) {
            break;
          }
          const auto& ft0_future = ft0s.iteratorAt(ft0.globalIndex() + deltaIndex);
          deltaBC = ft0_future.bcId() - ft0.bcId();

          if (deltaBC < myMaxDeltaBCFT0) {
            std::bitset<8> fT0TriggersFuture = ft0_future.triggerMask();
            bool vertexFuture = fT0TriggersFuture[o2::ft0::Triggers::bitVertex];
            bool triggerAFuture = fT0TriggersFuture[o2::ft0::Triggers::bitA];
            bool triggerCFuture = fT0TriggersFuture[o2::ft0::Triggers::bitC];

            futureActivityFT0Vertex |= vertexFuture;
            futureActivityFT0TriggerA |= triggerAFuture;
            futureActivityFT0TriggerC |= triggerCFuture;
          }
        }

        histos.fill(HIST("FT0/hCounts"), 0);
        if ((pastActivityFT0TriggerA || futureActivityFT0TriggerA) == true || (pastActivityFT0TriggerC || futureActivityFT0TriggerC) == true) {
          histos.fill(HIST("FT0/hCounts"), 2);
        } else {
          histos.fill(HIST("FT0/bcVertexTriggerBothSidesPFP"), localBC);
        }
        if (pastActivityFT0Vertex == true || futureActivityFT0Vertex == true) {
          histos.fill(HIST("FT0/hCounts"), 1);
        } else {
          histos.fill(HIST("FT0/bcVertexTriggerPFP"), localBC);
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
          deltaBC = fv0.bcId() - fv0_past.bcId();

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
        if (pastActivityFV0Cen == true || futureActivityFV0Cen == true) {
          histos.fill(HIST("FV0/hCounts"), 1);
        } else {
          histos.fill(HIST("FV0/bcCenTriggerPFPCentral"), localBC);
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
