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
///
/// \brief A filter task for diffractive events
///        requires: EvSels, o2-analysis-event-selection
///                  TrackSelection, o2-analysis-trackselection
///                  TracksExtended, o2-analysis-trackextension
///                  pidTOF*, o2-analysis-pid-tof
///                  pidTPC*, o2-analysis-pid-tpc
///                  Timestamps, o2-analysis-timestamp
///        usage: o2-analysis-timestamp --aod-file AO2D.root |
///               o2-analysis-event-selection --processRun3 1 --isMC 1 |
///               o2-analysis-trackextension |
///               o2-analysis-trackselection --isRun3 |
///               o2-analysis-pid-tof |
///               o2-analysis-pid-tpc |
///               o2-analysis-diffraction-filter
/// \author P. Buehler , paul.buehler@oeaw.ac.at
/// \since June 1, 2021

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "CommonConstants/LHCConstants.h"

#include "cutHolder.h"
#include "diffractionSelectors.h"
#include "../filterTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// The associations between collsisions and BCs can be ambiguous.
// By default a collision is associated with the BC closest in time.
// The collision time t_coll is determined by the tracks which are used to
// reconstruct the vertex. t_coll has an uncertainty dt_coll.
// Any BC with a BC time t_BC falling within a time window of +- ndt*dt_coll
// around t_coll could potentially be the true BC. ndt is typically 4.

template <typename T>
T compatibleBCs(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, int ndt, T const& bcs)
{
  LOGF(debug, "Collision time / resolution [ns]: %f / %f", collision.collisionTime(), collision.collisionTimeRes());

  auto bcIter = collision.bc_as<T>();

  // due to the filling scheme the most probably BC may not be the one estimated from the collision time
  uint64_t mostProbableBC = bcIter.globalBC();
  uint64_t meanBC = mostProbableBC - std::lround(collision.collisionTime() / o2::constants::lhc::LHCBunchSpacingNS);
  int deltaBC = std::ceil(collision.collisionTimeRes() / o2::constants::lhc::LHCBunchSpacingNS * ndt);
  int64_t minBC = meanBC - deltaBC;
  uint64_t maxBC = meanBC + deltaBC;
  if (minBC < 0) {
    minBC = 0;
  }

  // find slice of BCs table with BC in [minBC, maxBC]
  int64_t maxBCId = bcIter.globalIndex();
  int moveCount = 0; // optimize to avoid to re-create the iterator
  while (bcIter != bcs.end() && bcIter.globalBC() <= maxBC && (int64_t)bcIter.globalBC() >= minBC) {
    LOGF(debug, "Table id %d BC %llu", bcIter.globalIndex(), bcIter.globalBC());
    maxBCId = bcIter.globalIndex();
    ++bcIter;
    ++moveCount;
  }

  bcIter.moveByIndex(-moveCount); // Move back to original position
  int64_t minBCId = collision.bcId();
  while (bcIter != bcs.begin() && bcIter.globalBC() <= maxBC && (int64_t)bcIter.globalBC() >= minBC) {
    LOGF(debug, "Table id %d BC %llu", bcIter.globalIndex(), bcIter.globalBC());
    minBCId = bcIter.globalIndex();
    --bcIter;
  }

  LOGF(debug, "  BC range: %i (%d) - %i (%d)", minBC, minBCId, maxBC, maxBCId);

  T slice{{bcs.asArrowTable()->Slice(minBCId, maxBCId - minBCId + 1)}, (uint64_t)minBCId};
  bcs.copyIndexBindings(slice);
  return slice;
}

// Run 3
struct DGFilterRun3 {

  // Productions
  Produces<aod::DiffractionFilters> filterTable;

  // configurable cutHolder
  MutableConfigurable<cutHolder> diffCuts{"DiffCuts", {}, "Diffractive events cuts"};

  // DG selector
  DGSelector dgSelector;

  // some general Collisions and Tracks filter

  using CCs = soa::Join<aod::Collisions, aod::EvSels>;
  using CC = CCs::iterator;
  using BCs = soa::Join<aod::BCs, aod::BcSels, aod::Run3MatchedToBCSparse>;
  using BC = BCs::iterator;
  using TCs = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection,
                        aod::pidTPCEl, aod::pidTPCMu, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::TOFSignal, aod::pidTOFEl, aod::pidTOFMu, aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr>;
  // using MFs = aod::MFTTracks;
  using FWs = aod::FwdTracks;

  void process(CC const& collision,
               BCs const& bcs,
               TCs& tracks,
               //               MFs& mfttracks,
               FWs& fwdtracks,
               aod::Zdcs& zdcs,
               aod::FT0s& ft0s,
               aod::FV0As& fv0as,
               aod::FV0Cs& fv0cs,
               aod::FDDs& fdds)
  {
    // nominal BC
    auto bc = collision.bc_as<BCs>();

    // obtain slice of compatible BCs
    auto bcRange = compatibleBCs(collision, diffCuts->NDtcoll(), bcs);
    LOGF(debug, "  Number of compatible BCs in +- %i dtcoll: %i", diffCuts->NDtcoll(), bcRange.size());

    // check that there are no FIT signals in any of the compatible BCs
    // Double Gap (DG) condition
    auto isDGEvent = true;
    for (auto& bc : bcRange) {
      if (bc.has_ft0() || bc.has_fv0a() || bc.has_fdd()) {
        isDGEvent = false;
        break;
      }
    }

    // additional cuts
    if (isDGEvent) {
      isDGEvent = dgSelector.IsSelected(diffCuts, collision, bc, bcRange, tracks, fwdtracks);
    }

    // fill filterTable
    if (isDGEvent) {
      LOGF(info, "This collision is a DG candidate!");
    }
    filterTable(isDGEvent);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DGFilterRun3>(cfgc, TaskName{"DGfilterRun3"}),
  };
}
