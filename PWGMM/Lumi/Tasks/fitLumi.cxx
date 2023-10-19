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
#include <iostream>

// author: arvind.khuntia@cern.ch

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using BCPattern = std::bitset<o2::constants::lhc::LHCMaxBunches>;
const int nBCsPerOrbit = o2::constants::lhc::LHCMaxBunches;

#define FIND_INDEX(vector, element) (std::find((vector)->begin(), (vector)->end(), (element)) != (vector)->end() ? std::distance((vector)->begin(), std::find((vector)->begin(), (vector)->end(), (element))) : -1)

struct fitLumi {
  int nTF = 0;
  HistogramRegistry registry;
  int nBins, nCollBins;
  int relTS;
  Configurable<ULong64_t> startTimeInS{"startTime", 1668079200, "unix start time"};
  Configurable<ULong64_t> endTimeInS{"endTime", 1668098700, "unix end time"};
  Configurable<std::vector<int>> collBCArray{"collBCArray", {1022, 1062, 1102, 1142, 2015, 2055, 2161, 2201, 2241, 2281, 2321, 2361, 2401, 2441, 2481, 2521, 2561, 2601, 2641, 2681}, "Colliding BC for the VdM"};

  void init(InitContext&)
  {
    nBins = ULong64_t(endTimeInS - startTimeInS) / 2.;
    nCollBins = collBCArray->size();
    // Hist for FT0
    registry.add("FT0/VtxTrig", "vertex trigger;ts (s); Counts", {HistType::kTH1F, {{nBins, 0., static_cast<double>(endTimeInS - startTimeInS)}}});
    registry.add("FT0/VtxTrigPerBC", "vertex trigger per BC;ts (s); Counts", {HistType::kTH2F, {{nBins, 0., static_cast<double>(endTimeInS - startTimeInS)}, {nCollBins - 1, 0., static_cast<double>(nCollBins)}}});
    registry.add("FT0/TF", "TF ;ts (s); Counts (nTF)", {HistType::kTH1F, {{nBins, 0., static_cast<double>(endTimeInS - startTimeInS)}}});
    registry.add("FT0/CollBCMap", "BC map for the colliding bcsr;BC entry; entries", {HistType::kTH1F, {{nCollBins - 1, 0., static_cast<double>(nCollBins)}}});
    // Hist for FDD
    registry.add("FDD/VtxTrig", "vertex trigger;ts (s); Counts", {HistType::kTH1F, {{nBins, 0., static_cast<double>(endTimeInS - startTimeInS)}}});
    registry.add("FDD/VtxTrigPerBC", "vertex trigger per BC;ts (s); Counts", {HistType::kTH2F, {{nBins, 0., static_cast<double>(endTimeInS - startTimeInS)}, {nCollBins - 1, 0., static_cast<double>(nCollBins)}}});
  }

  using BCsWithTimestamps = soa::Join<aod::BCs, aod::Timestamps>;
  void process(aod::FT0s const& ft0s, aod::FDDs const& fdds, aod::BCsWithTimestamps const&)
  {
    if (nTF < 1) {
      for (int iBin = 0; iBin < nCollBins; iBin++) {
        registry.get<TH1>(HIST("FT0/CollBCMap"))->SetBinContent(iBin + 1, collBCArray->at(iBin));
        std::cout << "bin " << iBin << " value " << collBCArray->at(iBin) << std::endl;
      }
    }

    for (auto const& ft0 : ft0s) {
      auto bc = ft0.bc_as<BCsWithTimestamps>();
      if (!bc.timestamp())
        continue;
      std::bitset<8> fT0Triggers = ft0.triggerMask();
      [[maybe_unused]] auto localBC = bc.globalBC() % nBCsPerOrbit;
      int index = FIND_INDEX(collBCArray, localBC);
      bool vertex = fT0Triggers[o2::fdd::Triggers::bitVertex];
      auto tsInSecond = ((bc.timestamp() * 1.e-3) - startTimeInS); // covert ts from ms to second
      if (vertex) {
        registry.get<TH1>(HIST("FT0/VtxTrig"))->Fill(tsInSecond);
        if (index) {
          registry.get<TH2>(HIST("FT0/VtxTrigPerBC"))->Fill(tsInSecond, index);
        }
      } // vertex
    }   // ft0

    for (auto const& fdd : fdds) {
      auto bc = fdd.bc_as<BCsWithTimestamps>();
      if (!bc.timestamp())
        continue;
      std::bitset<8> fddTriggers = fdd.triggerMask();
      [[maybe_unused]] auto localBC = bc.globalBC() % nBCsPerOrbit;
      int index = FIND_INDEX(collBCArray, localBC);
      bool vertex = fddTriggers[o2::fdd::Triggers::bitVertex];
      auto tsInSecond = ((bc.timestamp() * 1.e-3) - startTimeInS); // covert ts from ms to second
      if (vertex) {
        registry.get<TH1>(HIST("FDD/VtxTrig"))->Fill(tsInSecond);
        if (index) {
          registry.get<TH2>(HIST("FDD/VtxTrigPerBC"))->Fill(tsInSecond, index);
        }
      } // vertex
    }   // fdd

    auto timeFirstInTFinS = ft0s.iteratorAt(0).bc_as<BCsWithTimestamps>().timestamp() * 1.e-3;
    registry.get<TH1>(HIST("FT0/TF"))->Fill(timeFirstInTFinS - startTimeInS);
    nTF++;
  } // process
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<fitLumi>(cfgc, TaskName{"ft0qavdm"})};
}
