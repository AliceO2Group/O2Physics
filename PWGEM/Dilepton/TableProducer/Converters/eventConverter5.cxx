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
// This code runs loop over ULS ee pars for virtual photon QC.
//    Please write to: daiki.sekihata@cern.ch

#include "PWGEM/Dilepton/DataModel/dileptonTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct eventConverter5 {
  Produces<aod::EMEvents_005> event_005;

  // template <typename TBC>
  // uint32_t reduceSelectionBit(TBC const& bc)
  // {
  //   //input should be aod::BcSels or aod::EvSels.
  //   uint32_t bitMap = 0;
  //   if (bc.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
  //     SETBIT(bitMap, o2::aod::emevsel::kIsTriggerTVX);
  //   }
  //   if (bc.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
  //     SETBIT(bitMap, o2::aod::emevsel::kNoTimeFrameBorder);
  //   }
  //   if (bc.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
  //     SETBIT(bitMap, o2::aod::emevsel::kNoITSROFrameBorder);
  //   }
  //   if (bc.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
  //     SETBIT(bitMap, o2::aod::emevsel::kNoSameBunchPileup);
  //   }
  //   if (bc.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
  //     SETBIT(bitMap, o2::aod::emevsel::kIsGoodZvtxFT0vsPV);
  //   }
  //   if (bc.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
  //     SETBIT(bitMap, o2::aod::emevsel::kIsVertexITSTPC);
  //   }
  //   if (bc.selection_bit(o2::aod::evsel::kIsVertexTRDmatched)) {
  //     SETBIT(bitMap, o2::aod::emevsel::kIsVertexTRDmatched);
  //   }
  //   if (bc.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
  //     SETBIT(bitMap, o2::aod::emevsel::kIsVertexTOFmatched);
  //   }
  //   if (bc.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
  //     SETBIT(bitMap, o2::aod::emevsel::kNoCollInTimeRangeStandard);
  //   }
  //   if (bc.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict)) {
  //     SETBIT(bitMap, o2::aod::emevsel::kNoCollInTimeRangeStandard);
  //   }
  //   if (bc.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow)) {
  //     SETBIT(bitMap, o2::aod::emevsel::kNoCollInTimeRangeNarrow);
  //   }
  //   if (bc.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
  //     SETBIT(bitMap, o2::aod::emevsel::kNoCollInRofStandard);
  //   }
  //   if (bc.selection_bit(o2::aod::evsel::kNoCollInRofStrict)) {
  //     SETBIT(bitMap, o2::aod::emevsel::kNoCollInRofStrict);
  //   }
  //   if (bc.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
  //     SETBIT(bitMap, o2::aod::emevsel::kNoHighMultCollInPrevRof);
  //   }
  //   if (bc.selection_bit(o2::aod::evsel::kIsGoodITSLayer3)) {
  //     SETBIT(bitMap, o2::aod::emevsel::kIsGoodITSLayer3);
  //   }
  //   if (bc.selection_bit(o2::aod::evsel::kIsGoodITSLayer0123)) {
  //     SETBIT(bitMap, o2::aod::emevsel::kIsGoodITSLayer0123);
  //   }
  //   if (bc.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
  //     SETBIT(bitMap, o2::aod::emevsel::kIsGoodITSLayersAll);
  //   }
  //   return bitMap;
  // }

  void process004to005(aod::EMEvents_004 const& collisions)
  {
    for (const auto& collision : collisions) {
      event_005(
        collision.globalIndex(),
        collision.runNumber(),
        collision.globalBC(),
        o2::aod::emevsel::reduceSelectionBit(collision),
        collision.rct_raw(),
        collision.timestamp(),
        collision.posZ(),
        collision.numContrib(),
        collision.trackOccupancyInTimeRange(),
        collision.ft0cOccupancyInTimeRange());
    } // end of collision loop
  }
  PROCESS_SWITCH(eventConverter5, process004to005, "convert from 004 into 005", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<eventConverter5>(cfgc, TaskName{"event-converter5"})};
}
