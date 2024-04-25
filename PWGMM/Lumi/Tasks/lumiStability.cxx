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
#include "Framework/ASoA.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::evsel;

using BCsWithTimestamps = soa::Join<aod::BCs, aod::Timestamps>;
using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;

int nBCsPerOrbit = 3564;

struct lumiStabilityTask {
  // Histogram registry: an object to hold your histograms
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    const AxisSpec axisTriggger{nBCsPerOrbit, -0.5f, nBCsPerOrbit - 0.5f};

    // histo about triggers
    histos.add("FDD/bcVertexTrigger", "vertex trigger per BC (FDD);BC in FDD; counts", kTH1F, {axisTriggger});
    histos.add("FDD/bcVertexTriggerCoincidence", "vertex trigger per BC (FDD) with coincidences;BC in FDD; counts", kTH1F, {axisTriggger});
    histos.add("FDD/bcSCentralTrigger", "scentral trigger per BC (FDD);BC in FDD; counts", kTH1F, {axisTriggger});
    histos.add("FDD/bcSCentralTriggerCoincidence", "scentral trigger per BC (FDD) with coincidences;BC in FDD; counts", kTH1F, {axisTriggger});
    histos.add("FDD/bcVSCTrigger", "vertex and scentral trigger per BC (FDD);BC in FDD; counts", kTH1F, {axisTriggger});
    histos.add("FDD/bcVSCTriggerCoincidence", "vertex and scentral trigger per BC (FDD) with coincidences;BC in FDD; counts", kTH1F, {axisTriggger});
    histos.add("FDD/bcCentralTrigger", "central trigger per BC (FDD);BC in FDD; counts", kTH1F, {axisTriggger});
    histos.add("FDD/bcCentralTriggerCoincidence", "central trigger per BC (FDD) with coincidences;BC in FDD; counts", kTH1F, {axisTriggger});
    histos.add("FDD/bcVCTrigger", "vertex and central trigger per BC (FDD);BC in FDD; counts", kTH1F, {axisTriggger});
    histos.add("FDD/bcVCTriggerCoincidence", "vertex and central trigger per BC (FDD) with coincidences;BC in FDD; counts", kTH1F, {axisTriggger});

    histos.add("FT0/bcVertexTrigger", "vertex trigger per BC (FT0);BC in FT0; counts", kTH1F, {axisTriggger});
    histos.add("FT0/bcSCentralTrigger", "Scentral trigger per BC (FT0);BC in FT0; counts", kTH1F, {axisTriggger});
    histos.add("FT0/bcVSCTrigger", "vertex and Scentral trigger per BC (FT0);BC in FT0; counts", kTH1F, {axisTriggger});
    histos.add("FT0/bcCentralTrigger", "central trigger per BC (FT0);BC in FT0; counts", kTH1F, {axisTriggger});
    histos.add("FT0/bcVCTrigger", "vertex and central trigger per BC (FT0);BC in FT0; counts", kTH1F, {axisTriggger});
    histos.add("FT0/bcSCentralCentralTrigger", "Scentral and central trigger per BC (FT0);BC in FT0; counts", kTH1F, {axisTriggger});

    histos.add("FV0/bcOutTrigger", "Out trigger per BC (FV0);BC in V0; counts", kTH1F, {axisTriggger});
    histos.add("FV0/bcInTrigger", "In trigger per BC (FV0);BC in V0; counts", kTH1F, {axisTriggger});
    histos.add("FV0/bcSCenTrigger", "SCen trigger per BC (FV0);BC in V0; counts", kTH1F, {axisTriggger});
    histos.add("FV0/bcCenTrigger", "Out trigger per BC (FV0);BC in V0; counts", kTH1F, {axisTriggger});
    histos.add("FV0/bcSCenCenTrigger", "SCen and Cen trigger per BC (FV0);BC in V0; counts", kTH1F, {axisTriggger});
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

  void processMain(aod::FT0s const& ft0s, aod::FDDs const& fdds, aod::FV0As const& fv0s, BCsRun3 const& bcs)
  {
    int CountNormal(0), CountPastProtec(0);
    for (auto const& fdd : fdds) {
      auto bc = fdd.bc_as<BCsRun3>();
      if (bc.timestamp() == 0) {
        continue;
      }

      Long64_t globalBC = bc.globalBC();
      int localBC = globalBC % nBCsPerOrbit;

      int deltaIndex = 0; // backward move counts
      int deltaBC = 0;    // current difference wrt globalBC
      int maxDeltaBC = 5; // maximum difference
      bool pastActivityFDD = false;
      while (deltaBC < maxDeltaBC) {
        deltaIndex++;
        if (bc.globalIndex() - deltaIndex < 0) {
          break;
        }
        const auto& bc_past = bcs.iteratorAt(bc.globalIndex() - deltaIndex);
        deltaBC = globalBC - bc_past.globalBC();
        if (deltaBC < maxDeltaBC) {
          pastActivityFDD |= bc_past.has_fdd();
        }
      }
      deltaIndex = 0;
      deltaBC = 0;

      bool futureActivityFDD = false;
      while (deltaBC < maxDeltaBC) {
        deltaIndex++;
        if (bc.globalIndex() + deltaIndex >= bcs.size()) {
          break;
        }
        const auto& bc_future = bcs.iteratorAt(bc.globalIndex() + deltaIndex);
        deltaBC = bc_future.globalBC() - globalBC;
        if (deltaBC < maxDeltaBC) {
          futureActivityFDD |= bc_future.has_fdd();
        }
      }

      CountNormal++;
      if (pastActivityFDD == true || futureActivityFDD == true) {
        CountPastProtec++;
        continue;
      }
      /*if (pastActivityFDD == true) {
        CountPastProtec++;
        continue;
      }*/

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
      auto bc = ft0.bc_as<BCsRun3>();
      if (bc.timestamp() == 0) {
        continue;
      }

      Long64_t globalBC = bc.globalBC();
      int localBC = globalBC % nBCsPerOrbit;

      int deltaIndex = 0; // backward move counts
      int deltaBC = 0;    // current difference wrt globalBC
      int maxDeltaBC = 5; // maximum difference
      bool pastActivityFT0 = false;
      while (deltaBC < maxDeltaBC) {
        deltaIndex++;
        if (bc.globalIndex() - deltaIndex < 0) {
          break;
        }
        const auto& bc_past = bcs.iteratorAt(bc.globalIndex() - deltaIndex);
        deltaBC = globalBC - bc_past.globalBC();
        if (deltaBC < maxDeltaBC) {
          pastActivityFT0 |= bc_past.has_ft0();
        }
      }
      deltaIndex = 0;
      deltaBC = 0;

      bool futureActivityFT0 = false;
      while (deltaBC < maxDeltaBC) {
        deltaIndex++;
        if (bc.globalIndex() + deltaIndex >= bcs.size()) {
          break;
        }
        const auto& bc_future = bcs.iteratorAt(bc.globalIndex() + deltaIndex);
        deltaBC = bc_future.globalBC() - globalBC;
        if (deltaBC < maxDeltaBC) {
          futureActivityFT0 |= bc_future.has_ft0();
        }
      }

      if (pastActivityFT0 == true || futureActivityFT0 == true) {
        continue;
      }
      /*if (pastActivityFT0 == true) {
        continue;
      }*/

      std::bitset<8> fT0Triggers = ft0.triggerMask();
      bool vertex = fT0Triggers[o2::fdd::Triggers::bitVertex];

      if (vertex) {
        histos.fill(HIST("FT0/bcVertexTrigger"), localBC);
      } // vertex true

      bool sCentral = fT0Triggers[o2::fdd::Triggers::bitSCen];
      bool central = fT0Triggers[o2::fdd::Triggers::bitCen];

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
      auto bc = fv0.bc_as<BCsRun3>();
      if (bc.timestamp() == 0) {
        continue;
      }

      Long64_t globalBC = bc.globalBC();
      int localBC = globalBC % nBCsPerOrbit;

      int deltaIndex = 0; // backward move counts
      int deltaBC = 0;    // current difference wrt globalBC
      int maxDeltaBC = 5; // maximum difference
      bool pastActivityV0A = false;
      while (deltaBC < maxDeltaBC) {
        deltaIndex++;
        if (bc.globalIndex() - deltaIndex < 0) {
          break;
        }
        const auto& bc_past = bcs.iteratorAt(bc.globalIndex() - deltaIndex);
        deltaBC = globalBC - bc_past.globalBC();
        if (deltaBC < maxDeltaBC) {
          pastActivityV0A |= bc_past.has_fv0a();
        }
      }
      deltaIndex = 0;
      deltaBC = 0;

      bool futureActivityV0A = false;
      while (deltaBC < maxDeltaBC) {
        deltaIndex++;
        if (bc.globalIndex() + deltaIndex >= bcs.size()) {
          break;
        }
        const auto& bc_future = bcs.iteratorAt(bc.globalIndex() + deltaIndex);
        deltaBC = bc_future.globalBC() - globalBC;
        if (deltaBC < maxDeltaBC) {
          futureActivityV0A |= bc_future.has_fv0a();
        }
      }

      if (pastActivityV0A == true || futureActivityV0A == true) {
        continue;
      }

      /*if (pastActivityV0A == true) {
        continue;
      }*/

      std::bitset<8> fv0Triggers = fv0.triggerMask();
      bool aOut = fv0Triggers[o2::fdd::Triggers::bitAOut];
      bool aIn = fv0Triggers[o2::fdd::Triggers::bitAIn];
      bool aSCen = fv0Triggers[o2::fdd::Triggers::bitTrgNchan];
      bool aCen = fv0Triggers[o2::fdd::Triggers::bitTrgCharge];

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
        if (aSCen) {
          histos.fill(HIST("FV0/bcSCenCenTrigger"), localBC);
        }
      }
    } // loop over V0 events
    std::cout << "************ >>>>>>>>>>>>>> "
              << "Whithout Past Protection: " << CountNormal << " "
              << "Avoided Whith Past Protection: " << CountPastProtec << "<<<<<<<<<<<< ********************" << std::endl;
  } // end processMain

  PROCESS_SWITCH(lumiStabilityTask, processMain, "Process FDD and FT0 to lumi stability analysis", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lumiStabilityTask>(cfgc)};
}
