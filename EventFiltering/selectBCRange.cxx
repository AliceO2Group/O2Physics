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

#include "Common/DataModel/EventSelection.h"
#include "CommonConstants/LHCConstants.h"
#include "CommonDataFormat/InteractionRecord.h"
#include "CommonDataFormat/IRFrame.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Logger.h"
#include "Framework/runDataProcessing.h"

#include "filterTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// .............................................................................
// Run 3
struct BCRangeSelector {

  Configurable<int> nTimeRes{"nTimeRes", 4, "Range to consider for search of compatible BCs in units of vertex-time-resolution."};
  Configurable<int> nMinBCs{"nMinBCs", 7, "Minimum width of time window to consider for search of compatible BCs in units of 2*BunchSpacing."};
  Configurable<double> fillFac{"fillFactor", 0.0, "Factor of MB events to add"};

  using FDs = aod::CefpDecisions;
  using CCs = soa::Join<aod::Collisions, aod::EvSels>;
  using CC = CCs::iterator;
  using BCs = soa::Join<aod::BCs, aod::BcSels, aod::Run3MatchedToBCSparse>;
  using BC = BCs::iterator;

  uint64_t nBCs, nCompBCs, nNotCompBCs;
  uint64_t clast, cnew;

  // buffer for task output
  std::vector<o2::dataformats::IRFrame> res;
  Produces<aod::BCRanges> tags;

  void init(o2::framework::InitContext&)
  {
    res.clear();
  }

  void run(ProcessingContext& pc)
  {
    // get and prepare the input tables
    auto TabConsumer1 = pc.inputs().get<TableConsumer>(aod::MetadataTrait<std::decay_t<aod::BCs>>::metadata::tableLabel());
    auto t1{TabConsumer1->asArrowTable()};
    auto TabConsumer2 = pc.inputs().get<TableConsumer>(aod::MetadataTrait<std::decay_t<aod::BcSels>>::metadata::tableLabel());
    auto t2{TabConsumer2->asArrowTable()};
    auto TabConsumer3 = pc.inputs().get<TableConsumer>(aod::MetadataTrait<std::decay_t<aod::Run3MatchedToBCSparse>>::metadata::tableLabel());
    auto t3{TabConsumer3->asArrowTable()};
    auto TabConsumer4 = pc.inputs().get<TableConsumer>(aod::MetadataTrait<std::decay_t<aod::Collisions>>::metadata::tableLabel());
    auto t4{TabConsumer4->asArrowTable()};
    auto TabConsumer5 = pc.inputs().get<TableConsumer>(aod::MetadataTrait<std::decay_t<aod::EvSels>>::metadata::tableLabel());
    auto t5{TabConsumer5->asArrowTable()};
    auto TabConsumer6 = pc.inputs().get<TableConsumer>(aod::MetadataTrait<std::decay_t<aod::CefpDecisions>>::metadata::tableLabel());
    auto t6{TabConsumer6->asArrowTable()};

    // join tables
    auto bcs = BCs({t1, t2, t3});
    auto cols = CCs({t4, t5});
    cols.bindExternalIndices(&bcs);
    FDs fdecs{{t6}};
    if (cols.size() != fdecs.size()) {
      throw std::runtime_error("Collision table and CefpDecision do not have the same number of rows! ");
    }

    // 1. loop over collisions
    auto filt = fdecs.begin();
    std::vector<std::pair<uint64_t, uint64_t>> bcRanges;
    for (auto collision : cols) {
      if (filt.hasCefpSelected()) {

        LOGF(debug, "Collision time / resolution [ns]: %f / %f", collision.collisionTime(), collision.collisionTimeRes());

        // return if collisions has no associated BC
        if (!collision.has_foundBC()) {
          LOGF(warning, "No compatible BCs found for collision that the framework assigned to BC %i", filt.hasGlobalBCId());
          filt++;
          continue;
        }

        // get associated BC
        auto bcIter = collision.foundBC_as<BCs>();

        // due to the filling scheme the most probable BC may not be the one estimated from the collision time
        InteractionRecord mostProbableBC;
        mostProbableBC.setFromLong(bcIter.globalBC());
        InteractionRecord meanBC = mostProbableBC + std::lround(collision.collisionTime() / o2::constants::lhc::LHCBunchSpacingNS);

        // enforce minimum number for deltaBC
        int deltaBC = std::ceil(collision.collisionTimeRes() * nTimeRes / o2::constants::lhc::LHCBunchSpacingNS);
        if (deltaBC < nMinBCs) {
          deltaBC = nMinBCs;
        }
        LOGF(debug, "BC %d,  deltaBC %d", bcIter.globalIndex(), deltaBC);

        auto minBC = meanBC - deltaBC;
        auto maxBC = meanBC + deltaBC;

        uint64_t minBCId = bcIter.globalIndex();
        uint64_t maxBCId = bcIter.globalIndex();

        auto localIter = bcIter;
        while (localIter.globalIndex() > 0) {
          --localIter;
          if (localIter.globalBC() >= minBC.toLong()) {
            minBCId = localIter.globalIndex();
          } else {
            break;
          }
        }
        localIter = bcIter;
        while (localIter.globalIndex() < bcs.size()) {
          ++localIter;
          if (localIter.globalBC() <= maxBC.toLong()) {
            maxBCId = localIter.globalIndex();
          } else {
            break;
          }
        }

        bcRanges.push_back(std::make_pair(minBCId, maxBCId));
      }
      filt++;
    }

    /// We cannot merge the ranges in the previous loop because while collisions are sorted by time, the corresponding minBCs can be unsorted as the collision time resolution is not constant
    std::sort(bcRanges.begin(), bcRanges.end(), [](const std::pair<uint64_t, uint64_t>& a, const std::pair<uint64_t, uint64_t>& b) {
      return a.first < b.first;
    });
    std::vector<std::pair<uint64_t, uint64_t>> bcRangesMerged(1, bcRanges[0]);
    for (uint64_t iR{1}; iR < bcRanges.size(); ++iR) {
      if (bcRanges[iR - 1].second >= bcRanges[iR].first) {
        bcRangesMerged.back().second = bcRanges[iR].second;
      } else {
        bcRangesMerged.push_back(bcRanges[iR]);
      }
    }
    bcRanges.swap(bcRangesMerged);

    // 2. extend ranges
    int nBCselected{0};
    for (auto& range : bcRanges) {
      nBCselected += range.second - range.first + 1;
    }
    int nToBeAdded = std::ceil((bcs.size() - nBCselected) * fillFac);
    int nToBeAddedPerRange = std::ceil(float(nToBeAdded) / bcRanges.size());
    LOGF(debug, "Extending ranges by %d BCs (%d per selected range)", nToBeAdded, nToBeAddedPerRange);

    InteractionRecord IR1, IR2;
    uint64_t first{bcs.iteratorAt(bcRanges[0].first).globalBC()};
    IR1.setFromLong(first);
    while (nToBeAdded > 0 && bcRanges[0].first > 0) { /// TODO: decide if we want to extend the ranges to the beginning of the dataframe
      first = bcs.iteratorAt(bcRanges[0].first - 1).globalBC();
      IR2.setFromLong(first);
      if (IR1.differenceInBC(IR2) > o2::constants::lhc::LHCMaxBunches) { // protection against change of orbit in the DataFrame
        LOGF(debug, "Jump by more than one orbit detected!");
        break;
      } else {
        IR1.setFromLong(first);
      }
      bcRanges[0].first--;
      nToBeAdded--;
    }
    LOGF(debug, "Extending ranges by %d BCs (%d per selected range)", nToBeAdded, nToBeAddedPerRange);

    for (uint64_t iR{0}; iR < bcRanges.size() && nToBeAdded > 0; ++iR) {
      uint64_t second{bcs.iteratorAt(bcRanges[iR].second).globalBC()};
      IR2.setFromLong(second);
      for (int i{0}; i < nToBeAddedPerRange && nToBeAdded > 0; ++i) {
        if (bcRanges[iR].second < bcs.size() - 1) {
          second = bcs.iteratorAt(bcRanges[iR].second + 1).globalBC();
          IR1.setFromLong(second);
          if (IR1.differenceInBC(IR2) > o2::constants::lhc::LHCMaxBunches) { // protection against change of orbit in the DataFrame
            LOGF(debug, "Jump by more than one orbit detected!");
            break;
          } else {
            IR2.setFromLong(second);
          }
          bcRanges[iR].second++;
          nToBeAdded--;
        }
      }
    }
    LOGF(debug, "End extension, remaining to be added %d BCs", nToBeAdded);

    // fill res
    LOGF(debug, "Merged and extended sorted ranges");
    for (auto& range : bcRanges) {
      LOGF(debug, "  %i - %i", range.first, range.second);
      uint64_t first{bcs.iteratorAt(range.first).globalBC()}, second{bcs.iteratorAt(range.second).globalBC()};
      IR1.setFromLong(first);
      IR2.setFromLong(second);
      res.emplace_back(IR1, IR2);
      tags(first, second);
    }
    // make res an output
    pc.outputs().snapshot({"PPF", "IFRAMES", 0, Lifetime::Timeframe}, res);

    // clean up
    res.clear();
  }

  // need a trivial process method
  // the parameters determine the tables available in the input
  void process(CCs const& collisions, BCs const& bcs, FDs const& fdecisions)
  {
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  DataProcessorSpec spec{adaptAnalysisTask<BCRangeSelector>(cfgc, TaskName{"bc-ranges-selector-task"})};

  // add output
  spec.outputs.emplace_back("PPF", "IFRAMES", 0, Lifetime::Timeframe);
  LOGF(debug, "Output %i", spec.outputs.size());

  return WorkflowSpec{spec};
}
