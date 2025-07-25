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

// author: arvind.khuntia@cern.ch

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/EventSelection.h"
#include "DataFormatsFDD/Digit.h"
#include "DataFormatsFIT/Triggers.h"
#include "Common/DataModel/FT0Corrected.h"

#include "CCDB/CcdbApi.h"
#include "CommonDataFormat/BunchFilling.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "TH1F.h"
#include "TH2F.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using BCPattern = std::bitset<o2::constants::lhc::LHCMaxBunches>;
const int nBCsPerOrbit = o2::constants::lhc::LHCMaxBunches;

struct VdMAO2D {
  Configurable<std::vector<int>> collBCArray{"collBCArray", {-999, -999}, "Colliding BC for the VdM"};
  Configurable<int> fillNumber{"fillNumber", -1, "fill nuber to select the SOR and EOR ts"};

  HistogramRegistry registry;
  BCPattern CollidingBunch;
  std::vector<int> collBCArrayFromCCDB;
  int newRunNumber = -999;
  int oldRunNumber = -999;
  int nTF = 0;
  int nBins, nCollBins;
  int relTS;
  uint64_t startTimeInS, endTimeInS;

  void init(InitContext&)
  {
    if (fillNumber == 8379) { /// VdM Nov22, pp 13.6 TeV
      startTimeInS = 1668079200;
      endTimeInS = 1668098700;
    } else if (fillNumber == 9126) { /// VdM 23, pp 13.6 TeV
      startTimeInS = 1694118600;
      endTimeInS = 1694138220;
    } else if (fillNumber == 9140) { /// VdM 23, Pb-Pb 5.36 ATeV
      startTimeInS = 1696962082;
      endTimeInS = 1696973263;
    }
    if (!fillNumber) {
      LOG(info) << " No valid fill number";
      return;
    }

    nBins = static_cast<int>(endTimeInS - startTimeInS);
    nCollBins = collBCArray->size();
    // Hist for FT0
    registry.add("FT0/VtxTrig", "vertex trigger;ts (s); Counts", {HistType::kTH1F, {{nBins, 0., static_cast<double>(endTimeInS - startTimeInS)}}});
    registry.add("FT0/VtxTrigPerBC", "vertex trigger per BC;ts (s); Counts", {HistType::kTH2F, {{nBins, 0., static_cast<double>(endTimeInS - startTimeInS)}, {nCollBins - 1, 0., static_cast<double>(nCollBins)}}});
    registry.add("FT0/TF", "TF ;ts (s); Counts (nTF)", {HistType::kTH1F, {{nBins, 0., static_cast<double>(endTimeInS - startTimeInS)}}});
    registry.add("FT0/CollBCMap", "BC map for the colliding bcsr;BC entry; entries", {HistType::kTH1F, {{nCollBins - 1, 0., static_cast<double>(nCollBins)}}});
    registry.add("FT0/bcVertex", "BC distribution (FT0 vertex);ts (s); Counts", {HistType::kTH1F, {{nBCsPerOrbit, 0., nBCsPerOrbit}}});
    registry.add("FT0/bcVertexCollBC", "BC distribution (FT0 vertex-CollBC);ts (s); Counts", {HistType::kTH1F, {{nBCsPerOrbit, 0., nBCsPerOrbit}}});
    // Hist for FDD
    registry.add("FDD/VtxTrig", "vertex trigger;ts (s); Counts", {HistType::kTH1F, {{nBins, 0., static_cast<double>(endTimeInS - startTimeInS)}}});
    registry.add("FDD/TrigCoincidence", "Coincidence trigger;ts (s); Counts", {HistType::kTH1F, {{nBins, 0., static_cast<double>(endTimeInS - startTimeInS)}}});
    registry.add("FDD/VtxTrigPerBC", "vertex trigger per BC;ts (s); Counts", {HistType::kTH2F, {{nBins, 0., static_cast<double>(endTimeInS - startTimeInS)}, {nCollBins - 1, 0., static_cast<double>(nCollBins)}}});
    registry.add("FDD/CoincidenceTrigPerBC", "coincidence trigger per BC;ts (s); Counts", {HistType::kTH2F, {{nBins, 0., static_cast<double>(endTimeInS - startTimeInS)}, {nCollBins - 1, 0., static_cast<double>(nCollBins)}}});
    registry.add("FDD/bcVertex", "BC distribution (FDD vertex);ts (s); Counts", {HistType::kTH1F, {{nBCsPerOrbit, 0., nBCsPerOrbit}}});
    registry.add("FDD/bcVertexCollBC", "BC distribution (FDD vertex-CollBC);ts (s); Counts", {HistType::kTH1F, {{nBCsPerOrbit, 0., nBCsPerOrbit}}});
  }

  using BCsWithTimestamps = soa::Join<aod::BCs, aod::Timestamps>;
  void processFT0(aod::FT0s const& ft0s, aod::BCsWithTimestamps const&)
  {
    if (nTF < 1) {
      for (int iBin = 0; iBin < nCollBins; iBin++) {
        registry.get<TH1>(HIST("FT0/CollBCMap"))->SetBinContent(iBin + 1, collBCArray->at(iBin));
        LOG(debug) << "bin " << iBin << " value " << collBCArray->at(iBin);
      }
    }

    for (auto const& ft0 : ft0s) {
      auto bc = ft0.bc_as<BCsWithTimestamps>();
      if (!bc.timestamp()) {
        continue;
      }
      std::bitset<8> fT0Triggers = ft0.triggerMask();
      auto localBC = bc.globalBC() % nBCsPerOrbit;
      auto pos = std::find(collBCArray->begin(), collBCArray->end(), localBC);
      bool vertex = fT0Triggers[o2::fdd::Triggers::bitVertex];
      auto tsInSecond = ((bc.timestamp() * 1.e-3) - startTimeInS); // convert ts from ms to second
      if (vertex) {
        registry.get<TH1>(HIST("FT0/VtxTrig"))->Fill(tsInSecond);
        if (tsInSecond > 1600 && tsInSecond < 2400)
          registry.get<TH1>(HIST("FT0/bcVertex"))->Fill(localBC);
        if (pos != collBCArray->end()) {
          registry.get<TH2>(HIST("FT0/VtxTrigPerBC"))->Fill(tsInSecond, 1 + (std::distance(collBCArray->begin(), pos)));
          registry.get<TH1>(HIST("FT0/bcVertexCollBC"))->Fill(localBC);
        }
      } // vertex
    }   // ft0
    nTF++;
  } // process

  PROCESS_SWITCH(VdMAO2D, processFT0, "Process FT0 trigger rates for VdM", true);

  bool checkAnyCoincidence(const std::vector<int>& channels)
  {
    std::map<int, int> channelPairs = {{0, 4}, {1, 5}, {2, 6}, {3, 7}};
    for (const auto& pair : channelPairs) {
      if (std::find(channels.begin(), channels.end(), pair.first) != channels.end() &&
          std::find(channels.begin(), channels.end(), pair.second) != channels.end()) {
        return true;
      }
    }
    return false;
  }

  void processFDD(aod::FDDs const& fdds, aod::BCsWithTimestamps const&)
  {
    for (auto const& fdd : fdds) {
      auto bc = fdd.bc_as<BCsWithTimestamps>();
      if (!bc.timestamp()) {
        continue;
      }
      std::bitset<8> fddTriggers = fdd.triggerMask();
      auto localBC = bc.globalBC() % nBCsPerOrbit;
      auto pos = std::find(collBCArray->begin(), collBCArray->end(), localBC);
      bool vertex = fddTriggers[o2::fdd::Triggers::bitVertex];
      auto tsInSecond = ((bc.timestamp() * 1.e-3) - startTimeInS); // convert ts from ms to second
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
        registry.get<TH1>(HIST("FDD/VtxTrig"))->Fill(tsInSecond);
        if (isCoinA && isCoinC) {
          registry.get<TH1>(HIST("FDD/TrigCoincidence"))->Fill(tsInSecond);
        }
        if (tsInSecond > 1600 && tsInSecond < 2400)
          registry.get<TH1>(HIST("FDD/bcVertex"))->Fill(localBC);
        if (pos != collBCArray->end()) {
          registry.get<TH2>(HIST("FDD/VtxTrigPerBC"))->Fill(tsInSecond, 1 + (std::distance(collBCArray->begin(), pos)));
          if (isCoinA && isCoinC) {
            registry.get<TH2>(HIST("FDD/CoincidenceTrigPerBC"))->Fill(tsInSecond, 1 + (std::distance(collBCArray->begin(), pos)));
          }
          registry.get<TH1>(HIST("FDD/bcVertexCollBC"))->Fill(localBC);
        }

      } // vertex
    }   // fdd
  }
  PROCESS_SWITCH(VdMAO2D, processFDD, "Process FDD trigger rates for VdM", true);

  void processBCTable(aod::BCsWithTimestamps const& bcWTimeStamps)
  {
    for (auto const& bc : bcWTimeStamps) {
      if (!bc.timestamp()) {
        continue;
      } else {
        auto timeFirstInTFinS = bc.timestamp() * 1.e-3;
        registry.get<TH1>(HIST("FT0/TF"))->Fill(timeFirstInTFinS - startTimeInS);
      }
    }
  }

  PROCESS_SWITCH(VdMAO2D, processBCTable, "Process BC table for VdM", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<VdMAO2D>(cfgc, TaskName{"VdMAO2D"})};
}
