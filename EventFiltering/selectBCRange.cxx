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
using o2::dataformats::IRFrame;

struct BCRangeSelector {

  Configurable<int> nTimeRes{"nTimeRes", 4, "Range to consider for search of compatible BCs in units of vertex-time-resolution."};
  Configurable<int> nMinBCs{"nMinBCs", 7, "Minimum width of time window to consider for search of compatible BCs in units of 2*BunchSpacing."};
  Configurable<double> fillFac{"fillFactor", 0.0, "Factor of MB events to add"};

  using CCs = soa::Join<aod::Collisions, aod::EvSels>;

  // buffer for task output
  Produces<aod::BCRanges> tags;

  template <typename T>
  IRFrame getIRFrame(T& collision)
  {
    auto collBC = collision.bc().globalBC();
    auto evSelBC = collision.has_foundBC() ? collision.foundBC().globalBC() : collBC;
    int deltaBC = std::ceil(collision.collisionTimeRes() * nTimeRes / constants::lhc::LHCBunchSpacingNS);
    deltaBC = std::max(deltaBC, nMinBCs.value);
    IRFrame bcRange{InteractionRecord::long2IR(std::min(collBC, evSelBC)), InteractionRecord::long2IR(std::max(collBC, evSelBC))};
    bcRange.getMax() += deltaBC;
    bcRange.getMin() -= deltaBC;
    return bcRange;
  }

  void run(ProcessingContext& pc)
  {
    auto bcConsumer = pc.inputs().get<TableConsumer>(o2::soa::getTableLabel<aod::BCs>());
    auto bcTable{bcConsumer->asArrowTable()};
    auto collConsumer = pc.inputs().get<TableConsumer>(o2::soa::getTableLabel<aod::Collisions>());
    auto collTable{collConsumer->asArrowTable()};
    auto evSelConsumer = pc.inputs().get<TableConsumer>(o2::soa::getTableLabel<aod::EvSels>());
    auto evSelTable{evSelConsumer->asArrowTable()};
    auto cefpConsumer = pc.inputs().get<TableConsumer>(o2::soa::getTableLabel<aod::CefpDecisions>());
    auto cefpTable{cefpConsumer->asArrowTable()};

    auto bcs = aod::BCs({bcTable});
    auto cols = CCs({collTable, evSelTable});
    cols.bindExternalIndices(&bcs);
    aod::CefpDecisions decisions{{cefpTable}};

    if (cols.size() != decisions.size()) {
      throw std::runtime_error("Collision table and CefpDecision do not have the same number of rows! ");
    }
    if (cols.size() == 0) {
      LOGF(warning, "No collisions found!");
      return;
    }

    auto filt = decisions.begin();
    int firstSelectedCollision{-1};
    std::vector<IRFrame> bcRanges;
    int nColl{0}, nSelected{0};
    for (auto collision : cols) {
      if (filt.cefpSelected0() || filt.cefpSelected1()) {
        if (firstSelectedCollision < 0) {
          firstSelectedCollision = nColl;
        }
        bcRanges.push_back(getIRFrame(collision));
        nSelected++;
      }
      nColl++;
      filt++;
    }

    if (bcRanges.empty()) {
      LOGF(warning, "No BCs selected!");
      return;
    }

    float fractionSelected{static_cast<float>(nSelected) / nColl};
    int nMB{std::min(static_cast<int>(fillFac * nColl) - nSelected, nColl - 1)};
    LOGF(info, "Selected %d collisions (%.2f%%) and %d MB events", nSelected, fractionSelected * 100, nMB);
    int maxCollisionId = std::max(nMB, firstSelectedCollision);
    int minCollisionId = (maxCollisionId == nMB) ? 0 : firstSelectedCollision - nMB;
    auto minCollision = cols.begin() + minCollisionId;
    IRFrame minFrame{getIRFrame(minCollision)};
    bcRanges[0].getMin() = std::min(bcRanges[0].getMin(), minFrame.getMin());
    if (maxCollisionId == nMB) {
      auto maxCollision = cols.begin() + nMB;
      IRFrame maxFrame{getIRFrame(maxCollision)};
      bcRanges[0].getMax() = std::max(bcRanges[0].getMax(), maxFrame.getMax());
    }

    /// We cannot merge the ranges in the previous loop because while collisions are sorted by time, the corresponding minBCs can be unsorted as the collision time resolution is not constant
    std::sort(bcRanges.begin(), bcRanges.end(), [](const IRFrame& a, const IRFrame& b) {
      return a.getMin() < b.getMin();
    });

    std::vector<IRFrame> bcRangesMerged(1, bcRanges[0]);
    for (uint64_t iR{1}; iR < bcRanges.size(); ++iR) {
      if (bcRangesMerged.back().getMax() >= bcRanges[iR].getMin()) {
        bcRangesMerged.back().getMax() = std::max(bcRangesMerged.back().getMax(), bcRanges[iR].getMax());
      } else {
        bcRangesMerged.push_back(bcRanges[iR]);
      }
    }
    bcRanges.swap(bcRangesMerged);

    for (auto& range : bcRanges) {
      tags(range.getMin().toLong(), range.getMax().toLong());
    }
  }

  // need a trivial process method: the parameters determine the tables available in the input
  void process(CCs const& /*collisions*/, aod::BCs const&, aod::CefpDecisions const&)
  {
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  DataProcessorSpec spec{adaptAnalysisTask<BCRangeSelector>(cfgc, TaskName{"bc-ranges-selector-task"})};
  return WorkflowSpec{spec};
}
