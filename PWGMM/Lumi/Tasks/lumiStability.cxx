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
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/ASoAHelpers.h"
#include "DataFormatsFDD/Digit.h"
#include "Framework/ASoA.h"

using namespace o2;
using namespace o2::framework;
using BCsWithTimestamps = soa::Join<aod::BCs, aod::Timestamps>;

int nBCsPerOrbit = 3564;

struct lumiStabilityTask {
  // Histogram registry: an object to hold your histograms
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    const AxisSpec axisFDDTriggger{nBCsPerOrbit, -0.5f, nBCsPerOrbit - 0.5f};
    const AxisSpec axisFT0Triggger{nBCsPerOrbit, -0.5f, nBCsPerOrbit - 0.5f};
    const AxisSpec axisFV0Triggger{nBCsPerOrbit, -0.5f, nBCsPerOrbit - 0.5f};

    // histo about triggers
    histos.add("FDD/bcVertexTrigger", "vertex trigger per BC (FDD);BC in FDD; counts", kTH1F, {axisFDDTriggger});
    histos.add("FDD/bcVertexTriggerCoincidence", "vertex trigger per BC (FDD) with coincidences;BC in FDD; counts", kTH1F, {axisFDDTriggger});
    histos.add("FDD/bcSCentralTrigger", "scentral trigger per BC (FDD);BC in FDD; counts", kTH1F, {axisFDDTriggger});
    histos.add("FDD/bcSCentralTriggerCoincidence", "scentral trigger per BC (FDD) with coincidences;BC in FDD; counts", kTH1F, {axisFDDTriggger});
    histos.add("FDD/bcVSCTrigger", "vertex and scentral trigger per BC (FDD);BC in FDD; counts", kTH1F, {axisFDDTriggger});
    histos.add("FDD/bcVSCTriggerCoincidence", "vertex and scentral trigger per BC (FDD) with coincidences;BC in FDD; counts", kTH1F, {axisFDDTriggger});
    histos.add("FDD/bcCentralTrigger", "central trigger per BC (FDD);BC in FDD; counts", kTH1F, {axisFDDTriggger});
    histos.add("FDD/bcCentralTriggerCoincidence", "central trigger per BC (FDD) with coincidences;BC in FDD; counts", kTH1F, {axisFDDTriggger});
    histos.add("FDD/bcVCTrigger", "vertex and central trigger per BC (FDD);BC in FDD; counts", kTH1F, {axisFDDTriggger});
    histos.add("FDD/bcVCTriggerCoincidence", "vertex and central trigger per BC (FDD) with coincidences;BC in FDD; counts", kTH1F, {axisFDDTriggger});

    histos.add("FT0/bcVertexTrigger", "vertex trigger per BC (FT0);BC in FT0; counts", kTH1F, {axisFT0Triggger});
    histos.add("FT0/bcSCentralTrigger", "Scentral trigger per BC (FT0);BC in FT0; counts", kTH1F, {axisFT0Triggger});
    histos.add("FT0/bcVSCTrigger", "vertex and Scentral trigger per BC (FT0);BC in FT0; counts", kTH1F, {axisFT0Triggger});
    histos.add("FT0/bcCentralTrigger", "central trigger per BC (FT0);BC in FT0; counts", kTH1F, {axisFT0Triggger});
    histos.add("FT0/bcVCTrigger", "vertex and central trigger per BC (FT0);BC in FT0; counts", kTH1F, {axisFT0Triggger});
    histos.add("FT0/bcSCentralCentralTrigger", "Scentral and central trigger per BC (FT0);BC in FT0; counts", kTH1F, {axisFT0Triggger});

    histos.add("FV0/bcOutTrigger", "Out trigger per BC (FV0);BC in V0; counts", kTH1F, {axisFV0Triggger});
    histos.add("FV0/bcInTrigger", "In trigger per BC (FV0);BC in V0; counts", kTH1F, {axisFV0Triggger});
    histos.add("FV0/bcSCenTrigger", "SCen trigger per BC (FV0);BC in V0; counts", kTH1F, {axisFV0Triggger});
    histos.add("FV0/bcCenTrigger", "Out trigger per BC (FV0);BC in V0; counts", kTH1F, {axisFV0Triggger});
    histos.add("FV0/bcSCenCenTrigger", "SCen and Cen trigger per BC (FV0);BC in V0; counts", kTH1F, {axisFV0Triggger});
  }

  bool checkAnyCoincidence(const std::vector<int>& channels)
  {
    std::map<int, int> channelPairs = {{0, 4}, {1, 5}, {2, 6}, {4, 7}};
    for (const auto& pair : channelPairs) {
      if (std::find(channels.begin(), channels.end(), pair.first) != channels.end() &&
          std::find(channels.begin(), channels.end(), pair.second) != channels.end()) {
        return true;
      }
    }
    return false;
  }

  void processFDDFT0(aod::FT0s const& ft0s, aod::FDDs const& fdds, aod::BCsWithTimestamps const&)
  {
    for (auto const& fdd : fdds) {
      auto bc = fdd.bc_as<BCsWithTimestamps>();
      if (bc.timestamp() == false) {
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
      if (bc.timestamp() == false) {
        continue;
      }

      Long64_t globalBC = bc.globalBC();
      int localBC = globalBC % nBCsPerOrbit;

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
  }   // end processFDDFT0

  PROCESS_SWITCH(lumiStabilityTask, processFDDFT0, "Process FDD and FT0 to lumi stability analysis", true);

  void processV0(aod::FV0As const& fv0s, aod::BCsWithTimestamps const&)
  {
    for (auto const& fv0 : fv0s) {
      auto bc = fv0.bc_as<BCsWithTimestamps>();
      if (bc.timestamp() == false) {
        continue;
      }

      Long64_t globalBC = bc.globalBC();
      int localBC = globalBC % nBCsPerOrbit;

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
  }   // end processV0

  PROCESS_SWITCH(lumiStabilityTask, processV0, "Process V0 to lumi stability analysis", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lumiStabilityTask>(cfgc)};
}
