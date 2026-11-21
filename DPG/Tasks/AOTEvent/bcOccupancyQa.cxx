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

/// \file bcOccupancyQa.cxx
/// \brief QA of the occupancy wrt bunch crossings: aod::BcOccsAll (all BCs of the TFs), produced by the event selection service
///
/// \author Igor Altsybeev <Igor.Altsybeev@cern.ch>

#include "Common/DataModel/EventSelection.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>

#include <array>
#include <cstdint>

using namespace o2;
using namespace o2::framework;

struct BcOccupancyQaTask {
  Configurable<int> nBinsOccupancyByTracks{"nBinsOccupancyByTracks", 500, "number of bins for the occupancy by ITS tracks"};
  Configurable<float> maxOccupancyByTracks{"maxOccupancyByTracks", 25000., "upper edge of the axis for the occupancy by ITS tracks"};
  Configurable<int> nBinsOccupancyByFT0C{"nBinsOccupancyByFT0C", 500, "number of bins for the occupancy by FT0C amplitudes"};
  Configurable<float> maxOccupancyByFT0C{"maxOccupancyByFT0C", 250000., "upper edge of the axis for the occupancy by FT0C amplitudes"};
  Configurable<int> nBinsNumCollsInRof{"nBinsNumCollsInRof", 50, "number of bins for the number of collisions in the ITS ROF (also the upper edge of the axis)"};
  Configurable<int> nBinsBcInTF{"nBinsBcInTF", 1188, "number of bins for the BC position in the TF"};
  Configurable<int> nOrbitsPerTF{"nOrbitsPerTF", 32, "number of orbits per TF, for the histograms vs BC position in the TF (128 in 2022)"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  static constexpr int NFlags = 7;
  static constexpr std::array<int, NFlags> FlagBits = {aod::bcocc::kNoCollInTimeRangeNarrow, aod::bcocc::kNoCollInTimeRangeStrict, aod::bcocc::kNoCollInTimeRangeStandard,
                                                       aod::bcocc::kNoCollInRofStrict, aod::bcocc::kNoCollInRofStandard, aod::bcocc::kNoHighMultCollInPrevRof, aod::bcocc::kUpcMode};
  static constexpr std::array<const char*, NFlags> FlagNames = {"NoCollInTimeRangeNarrow", "NoCollInTimeRangeStrict", "NoCollInTimeRangeStandard",
                                                                "NoCollInRofStrict", "NoCollInRofStandard", "NoHighMultCollInPrevRof", "UpcMode"};
  static constexpr int NBCsPerOrbit = 3564;
  int64_t nBCsPerTF = 0;

  void init(InitContext&)
  {
    nBCsPerTF = static_cast<int64_t>(nOrbitsPerTF.value) * NBCsPerOrbit;
    const AxisSpec axisOccupancyByTracks{nBinsOccupancyByTracks.value, 0., maxOccupancyByTracks.value, "occupancy by ITS tracks"};
    const AxisSpec axisOccupancyByFT0C{nBinsOccupancyByFT0C.value, 0., maxOccupancyByFT0C.value, "occupancy by FT0C amplitudes"};
    const AxisSpec axisMedianTime{280, -40., 100., "median time of the occupancy (#mus)"};
    const AxisSpec axisNumCollsInBC{11, -0.5, 10.5, "number of collisions found in the BC"};
    const AxisSpec axisNumCollsInRof{nBinsNumCollsInRof.value + 1, -1.5, nBinsNumCollsInRof.value - 0.5, "number of collisions found in the ITS ROF of the BC"};
    const AxisSpec axisBcInTF{nBinsBcInTF.value, 0., static_cast<double>(nBCsPerTF), "BC in TF"};
    const AxisSpec axisBcInOrbit{NBCsPerOrbit, -0.5, NBCsPerOrbit - 0.5, "BC in orbit"};
    const AxisSpec axisFlags{NFlags + 1, -0.5, NFlags + 0.5, ""};
    const AxisSpec axisBcCounter{3, -0.5, 2.5, ""};
    const AxisSpec axisIndexCheck{6, -0.5, 5.5, ""};

    // all BCs
    histos.add("hBcCounter", "number of BCs", kTH1D, {axisBcCounter});
    histos.add("hNumCollsInBC", "all BCs", kTH1D, {axisNumCollsInBC});
    histos.add("hNumCollsInRof", "all BCs", kTH1D, {axisNumCollsInRof});
    histos.add("hNumCollsInRofVsInBC", "all BCs", kTH2D, {axisNumCollsInBC, axisNumCollsInRof});
    histos.add("hNumCollsInRofUpcMode", "all BCs in ITS ROFs in UPC mode (kUpcMode)", kTH1D, {axisNumCollsInRof});
    histos.add("hBcInTF", "all BCs", kTH1D, {axisBcInTF});
    // BCs with calculated occupancy
    histos.add("hBcInTFOccupancyDefined", "BCs with calculated occupancy", kTH1D, {axisBcInTF});
    histos.add("hOccupancyByTracks", "BCs with calculated occupancy", kTH1D, {axisOccupancyByTracks});
    histos.add("hOccupancyByFT0C", "BCs with calculated occupancy", kTH1D, {axisOccupancyByFT0C});
    histos.add("hOccupancyByTracksVsFT0C", "BCs with calculated occupancy", kTH2D, {axisOccupancyByFT0C, axisOccupancyByTracks});
    histos.add("hMedianTime", "BCs with calculated occupancy", kTH1D, {axisMedianTime});
    histos.add("pOccupancyByTracksVsBcInTF", "BCs with calculated occupancy;BC in TF;mean occupancy by ITS tracks", kTProfile, {axisBcInTF});
    histos.add("pOccupancyByTracksVsBcInOrbit", "BCs with calculated occupancy;BC in orbit;mean occupancy by ITS tracks", kTProfile, {axisBcInOrbit});
    // occupancy flags: first bin - all BCs with calculated occupancy, then the BCs with a given flag set
    histos.add("hOccupancyFlags", "BCs with calculated occupancy", kTH1D, {axisFlags});
    histos.add("hOccupancyFlagsVsOccupancyByTracks", "BCs with calculated occupancy", kTH2D, {axisFlags, axisOccupancyByTracks});
    // cross-check of the lookup by global BC (aod::bcocc::BcOccAllIndex)
    histos.add("hBcIndexCheck", "lookup in BcOccsAll by global BC", kTH1D, {axisIndexCheck});

    auto hBcCounter = histos.get<TH1>(HIST("hBcCounter"));
    hBcCounter->GetXaxis()->SetBinLabel(1, "all");
    hBcCounter->GetXaxis()->SetBinLabel(2, "occupancy calculated");
    hBcCounter->GetXaxis()->SetBinLabel(3, "occupancy not available");
    auto hFlags = histos.get<TH1>(HIST("hOccupancyFlags"));
    auto hFlagsVsOcc = histos.get<TH2>(HIST("hOccupancyFlagsVsOccupancyByTracks"));
    hFlags->GetXaxis()->SetBinLabel(1, "all");
    hFlagsVsOcc->GetXaxis()->SetBinLabel(1, "all");
    for (int i = 0; i < NFlags; i++) {
      hFlags->GetXaxis()->SetBinLabel(i + 2, FlagNames[i]);
      hFlagsVsOcc->GetXaxis()->SetBinLabel(i + 2, FlagNames[i]);
    }
    auto hBcIndexCheck = histos.get<TH1>(HIST("hBcIndexCheck"));
    hBcIndexCheck->GetXaxis()->SetBinLabel(1, "BcOccsAll rows");
    hBcIndexCheck->GetXaxis()->SetBinLabel(2, "row found by its global BC");
    hBcIndexCheck->GetXaxis()->SetBinLabel(3, "BCs in aod::BCs");
    hBcIndexCheck->GetXaxis()->SetBinLabel(4, "BC found in BcOccsAll");
    hBcIndexCheck->GetXaxis()->SetBinLabel(5, "BC found, same global BC");
    hBcIndexCheck->GetXaxis()->SetBinLabel(6, "BC not found");

    // aod::BcOccs (joinable with aod::BCs)
    histos.add("BcOccs/hBcCounter", "number of BCs", kTH1D, {axisBcCounter});
    histos.add("BcOccs/hNumCollsInBC", "all BCs", kTH1D, {axisNumCollsInBC});
    histos.add("BcOccs/hNumCollsInRof", "all BCs", kTH1D, {axisNumCollsInRof});
    histos.add("BcOccs/hOccupancyByTracks", "BCs with calculated occupancy", kTH1D, {axisOccupancyByTracks});
    histos.add("BcOccs/hOccupancyByFT0C", "BCs with calculated occupancy", kTH1D, {axisOccupancyByFT0C});
    histos.add("BcOccs/hOccupancyByTracksVsFT0C", "BCs with calculated occupancy", kTH2D, {axisOccupancyByFT0C, axisOccupancyByTracks});
    histos.add("BcOccs/hOccupancyFlags", "BCs with calculated occupancy", kTH1D, {axisFlags});
    auto hBcCounterBcOccs = histos.get<TH1>(HIST("BcOccs/hBcCounter"));
    hBcCounterBcOccs->GetXaxis()->SetBinLabel(1, "all");
    hBcCounterBcOccs->GetXaxis()->SetBinLabel(2, "occupancy calculated");
    hBcCounterBcOccs->GetXaxis()->SetBinLabel(3, "occupancy not available");
    auto hFlagsBcOccs = histos.get<TH1>(HIST("BcOccs/hOccupancyFlags"));
    hFlagsBcOccs->GetXaxis()->SetBinLabel(1, "all");
    for (int i = 0; i < NFlags; i++) {
      hFlagsBcOccs->GetXaxis()->SetBinLabel(i + 2, FlagNames[i]);
    }

    // comparison of aod::BcOccs with aod::BcOccsAll for the BCs of aod::BCs
    histos.add("BcOccsVsBcOccsAll/hCheck", "BcOccs vs BcOccsAll", kTH1D, {{7, -0.5, 6.5, ""}});
    auto hCheck = histos.get<TH1>(HIST("BcOccsVsBcOccsAll/hCheck"));
    hCheck->GetXaxis()->SetBinLabel(1, "BCs in BcOccs");
    hCheck->GetXaxis()->SetBinLabel(2, "found in BcOccsAll");
    hCheck->GetXaxis()->SetBinLabel(3, "same numCollsInBC, numCollsInRof");
    hCheck->GetXaxis()->SetBinLabel(4, "same kUpcMode");
    hCheck->GetXaxis()->SetBinLabel(5, "occupancy calculated in both");
    hCheck->GetXaxis()->SetBinLabel(6, "same occupancy and flags");
    hCheck->GetXaxis()->SetBinLabel(7, "not found in BcOccsAll");
  }

  // QA of aod::BcOccsAll (all BCs of the TFs), with the cross-check of the lookup by global BC
  void processBcOccsAll(aod::BcOccsAll const& bcOccs, aod::BCs const& bcs)
  {
    // cross-check of the lookup by global BC: every BC of the BC table should be found, with the same global BC
    // (not found: BCs before the start of run, or no BcOccsAll entries for this DF)
    aod::bcocc::BcOccAllIndex bcOccAllIndexed(bcOccs);
    for (const auto& bc : bcs) {
      histos.fill(HIST("hBcIndexCheck"), 2);
      if (!bcOccAllIndexed.contains(bc.globalBC())) {
        histos.fill(HIST("hBcIndexCheck"), 5);
        continue;
      }
      histos.fill(HIST("hBcIndexCheck"), 3);
      if (bcOccAllIndexed[bc.globalBC()].globalBC() == bc.globalBC()) {
        histos.fill(HIST("hBcIndexCheck"), 4);
      }
    }

    if (bcOccs.size() == 0) {
      return;
    }
    // BC position in the TF: the table contains complete TFs and TF starts differ by multiples of nBCsPerTF,
    // so bcInTF = (globalBC - globalBC of the first row) % nBCsPerTF; check that nOrbitsPerTF is consistent with the table
    const uint64_t firstGlobalBC = bcOccs.iteratorAt(0).globalBC();
    const uint64_t lastGlobalBC = bcOccs.iteratorAt(bcOccs.size() - 1).globalBC();
    const bool isBcInTFValid = (bcOccs.size() % nBCsPerTF == 0) && (static_cast<int64_t>((lastGlobalBC - firstGlobalBC) % nBCsPerTF) == nBCsPerTF - 1);
    if (!isBcInTFValid) {
      LOGP(warning, "BC occupancy QA: table size {} is not consistent with {} orbits per TF, no histograms vs BC in TF for this DF", bcOccs.size(), nOrbitsPerTF.value);
    }

    for (const auto& bcOcc : bcOccs) {
      const uint64_t globalBC = bcOcc.globalBC();
      const int bcInTF = isBcInTFValid ? static_cast<int>((globalBC - firstGlobalBC) % nBCsPerTF) : -1;
      const int bcInOrbit = static_cast<int>(globalBC % NBCsPerOrbit);
      histos.fill(HIST("hBcCounter"), 0);
      histos.fill(HIST("hNumCollsInBC"), bcOcc.numCollsInBC());
      histos.fill(HIST("hNumCollsInRof"), bcOcc.numCollsInRof());
      histos.fill(HIST("hNumCollsInRofVsInBC"), bcOcc.numCollsInBC(), bcOcc.numCollsInRof());
      if (bcOcc.occupancyFlags_bit(aod::bcocc::kUpcMode)) { // set for all BCs, also without calculated occupancy
        histos.fill(HIST("hNumCollsInRofUpcMode"), bcOcc.numCollsInRof());
      }
      histos.fill(HIST("hBcIndexCheck"), 0);
      if (bcOccAllIndexed.contains(globalBC) && bcOccAllIndexed[globalBC].globalBC() == globalBC) { // global BCs are unique: same row
        histos.fill(HIST("hBcIndexCheck"), 1);
      }
      if (isBcInTFValid) {
        histos.fill(HIST("hBcInTF"), bcInTF);
      }

      // occupancy is not available (BC is close to TF borders, or not selected for the calculation)
      const int occByTracks = bcOcc.trackOccupancyInTimeRange();
      const float occByFT0C = bcOcc.ft0cOccupancyInTimeRange();
      if (occByTracks < 0 || occByFT0C < 0) {
        histos.fill(HIST("hBcCounter"), 2);
        continue;
      }
      histos.fill(HIST("hBcCounter"), 1);
      histos.fill(HIST("hOccupancyByTracks"), occByTracks);
      histos.fill(HIST("hOccupancyByFT0C"), occByFT0C);
      histos.fill(HIST("hOccupancyByTracksVsFT0C"), occByFT0C, occByTracks);
      histos.fill(HIST("hMedianTime"), bcOcc.occupancyMedianTime());
      histos.fill(HIST("pOccupancyByTracksVsBcInOrbit"), bcInOrbit, occByTracks);
      if (isBcInTFValid) {
        histos.fill(HIST("hBcInTFOccupancyDefined"), bcInTF);
        histos.fill(HIST("pOccupancyByTracksVsBcInTF"), bcInTF, occByTracks);
      }

      histos.fill(HIST("hOccupancyFlags"), 0);
      histos.fill(HIST("hOccupancyFlagsVsOccupancyByTracks"), 0, occByTracks);
      for (int i = 0; i < NFlags; i++) {
        if (bcOcc.occupancyFlags_bit(FlagBits[i])) {
          histos.fill(HIST("hOccupancyFlags"), i + 1);
          histos.fill(HIST("hOccupancyFlagsVsOccupancyByTracks"), i + 1, occByTracks);
        }
      }
    }
  }
  PROCESS_SWITCH(BcOccupancyQaTask, processBcOccsAll, "QA of aod::BcOccsAll (all BCs of the TFs)", true);

  // QA of aod::BcOccs (joinable with aod::BCs)
  void processBcOccs(soa::Join<aod::BCs, aod::BcOccs> const& bcs)
  {
    for (const auto& bc : bcs) {
      histos.fill(HIST("BcOccs/hBcCounter"), 0);
      histos.fill(HIST("BcOccs/hNumCollsInBC"), bc.numCollsInBC());
      histos.fill(HIST("BcOccs/hNumCollsInRof"), bc.numCollsInRof());
      const int occByTracks = bc.trackOccupancyInTimeRange();
      const float occByFT0C = bc.ft0cOccupancyInTimeRange();
      if (occByTracks < 0 || occByFT0C < 0) {
        histos.fill(HIST("BcOccs/hBcCounter"), 2);
        continue;
      }
      histos.fill(HIST("BcOccs/hBcCounter"), 1);
      histos.fill(HIST("BcOccs/hOccupancyByTracks"), occByTracks);
      histos.fill(HIST("BcOccs/hOccupancyByFT0C"), occByFT0C);
      histos.fill(HIST("BcOccs/hOccupancyByTracksVsFT0C"), occByFT0C, occByTracks);
      histos.fill(HIST("BcOccs/hOccupancyFlags"), 0);
      for (int i = 0; i < NFlags; i++) {
        if (bc.occupancyFlags_bit(FlagBits[i])) {
          histos.fill(HIST("BcOccs/hOccupancyFlags"), i + 1);
        }
      }
    }
  }
  PROCESS_SWITCH(BcOccupancyQaTask, processBcOccs, "QA of aod::BcOccs (joinable with aod::BCs)", true);

  // comparison of aod::BcOccs with aod::BcOccsAll for the BCs of aod::BCs: the same calculation, so the values should be identical
  // (occupancy and flags are compared only if calculated in both, since bcsForOccupancy can differ for the two tables)
  void processBcOccsVsBcOccsAll(soa::Join<aod::BCs, aod::BcOccs> const& bcs, aod::BcOccsAll const& bcOccsAll)
  {
    aod::bcocc::BcOccAllIndex bcOccAllIndexed(bcOccsAll);
    for (const auto& bc : bcs) {
      histos.fill(HIST("BcOccsVsBcOccsAll/hCheck"), 0);
      if (!bcOccAllIndexed.contains(bc.globalBC())) {
        histos.fill(HIST("BcOccsVsBcOccsAll/hCheck"), 6);
        continue;
      }
      histos.fill(HIST("BcOccsVsBcOccsAll/hCheck"), 1);
      const auto occAll = bcOccAllIndexed[bc.globalBC()];
      if (bc.numCollsInBC() == occAll.numCollsInBC() && bc.numCollsInRof() == occAll.numCollsInRof()) {
        histos.fill(HIST("BcOccsVsBcOccsAll/hCheck"), 2);
      }
      if (bc.occupancyFlags_bit(aod::bcocc::kUpcMode) == occAll.occupancyFlags_bit(aod::bcocc::kUpcMode)) {
        histos.fill(HIST("BcOccsVsBcOccsAll/hCheck"), 3);
      }
      if (bc.trackOccupancyInTimeRange() >= 0 && occAll.trackOccupancyInTimeRange() >= 0) {
        histos.fill(HIST("BcOccsVsBcOccsAll/hCheck"), 4);
        if (bc.trackOccupancyInTimeRange() == occAll.trackOccupancyInTimeRange() && bc.ft0cOccupancyInTimeRange() == occAll.ft0cOccupancyInTimeRange() &&
            bc.occupancyMedianTime() == occAll.occupancyMedianTime() && bc.occupancyFlags_raw() == occAll.occupancyFlags_raw()) {
          histos.fill(HIST("BcOccsVsBcOccsAll/hCheck"), 5);
        }
      }
    }
  }
  PROCESS_SWITCH(BcOccupancyQaTask, processBcOccsVsBcOccsAll, "comparison of aod::BcOccs with aod::BcOccsAll for the BCs of aod::BCs (requires both tables)", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<BcOccupancyQaTask>(cfgc)};
}
