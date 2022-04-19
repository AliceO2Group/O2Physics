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
// O2 includes

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "CommonConstants/LHCConstants.h"
#include "CommonDataFormat/InteractionRecord.h"
#include "CommonDataFormat/IRFrame.h"
#include "ReconstructionDataFormats/BCRange.h"

#include "filterTables.h"
#include "PWGUD/diffractionSelectors.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// .............................................................................
// Run 3
struct BCRangeSelector {

  Configurable<double> fillFac{"fillFactor", 0.0, "Factor of MB events to add"};

  using CCs = soa::Join<aod::Collisions, aod::EvSels>;
  using CC = CCs::iterator;
  using BBCs = soa::Join<aod::BCs, aod::BcSels, aod::Run3MatchedToBCSparse>;
  using BBC = BBCs::iterator;

  uint64_t nBCs, nCompBCs, nNotCompBCs;
  uint64_t clast, cnew;
  o2::dataformats::bcRanges cbcrs = o2::dataformats::bcRanges("Initial list"); // ranges of compatible BCs

  // buffer for task output
  std::vector<o2::dataformats::IRFrame> res;

  void init(o2::framework::InitContext&)
  {
    cbcrs.reset();
    res.clear();
  }

  void run(ProcessingContext& pc)
  {
    // get and prepare the input tables
    auto t1 = pc.inputs().get<TableConsumer>("BCs")->asArrowTable();
    auto t2 = pc.inputs().get<TableConsumer>("BcSels")->asArrowTable();
    auto t3 = pc.inputs().get<TableConsumer>("Run3MatchedToBCSparse")->asArrowTable();
    auto t4 = pc.inputs().get<TableConsumer>("Collisions")->asArrowTable();
    auto t5 = pc.inputs().get<TableConsumer>("EvSels")->asArrowTable();

    // join tables
    auto bcs = BBCs({t1, t2, t3});
    auto cols = CCs({t4, t5});
    cols.bindExternalIndices(&bcs);

    // 1. loop over collisions
    for (auto collision : cols) {
      auto bcRange = compatibleBCs(collision, 4, bcs);

      // update list of ranges
      auto bcfirst = bcRange.rawIteratorAt(0);
      auto bclast = bcRange.rawIteratorAt(bcRange.size());
      cbcrs.add(bcfirst.globalIndex(), bclast.globalIndex());
    }

    // 2. sort, merge, and extend ranges of compatible BCs
    //LOGP(info, "fillFac {}", fillFac);
    cbcrs.compact(bcs, fillFac);

    // fill res
    InteractionRecord IR1, IR2;
    LOGF(debug, "Merged and extended sorted ranges");
    for (auto limit : cbcrs.list()) {
      LOGF(debug, "  %i - %i", limit.first, limit.second);
      IR1.setFromLong(limit.first);
      IR2.setFromLong(limit.second);
      res.emplace_back(IR1, IR2);
    }
    // make res an output
    pc.outputs().snapshot({"PPF", "BCR", 0, Lifetime::Timeframe}, res);
    
    // clean up
    cbcrs.reset();
    res.clear();
  }

  // need a trivial process method
  // the parameters determine the tables available in the input
  void process(CCs const& collisions, BBCs const& bcs)
  {
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  DataProcessorSpec spec{adaptAnalysisTask<BCRangeSelector>(cfgc, TaskName{"bc-ranges-selector-task"})};

  // add output
  spec.outputs.emplace_back("PPF", "BCR", 0, Lifetime::Timeframe);
  LOGF(debug, "Output %i", spec.outputs.size());

  return WorkflowSpec{spec};
}
