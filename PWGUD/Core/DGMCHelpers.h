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
/// \brief
/// \author Paul Buehler, paul.buehler@oeaw.ac.at
/// \since  01.10.2021

#ifndef O2_ANALYSIS_DGMCHELPER_H_
#define O2_ANALYSIS_DGMCHELPER_H_

#include "Framework/Logger.h"
#include "CommonConstants/LHCConstants.h"
#include "Common/DataModel/EventSelection.h"

using namespace o2;
using namespace o2::framework;

// -----------------------------------------------------------------------------
// The associations between collisions and BCs can be ambiguous.
// By default a collision is associated with the BC closest in time.
// The collision time t_coll is determined by the tracks which are used to
// reconstruct the vertex. t_coll has an uncertainty dt_coll.
// Any BC with a BC time t_BC falling within a time window of +- ndt*dt_coll
// around t_coll could potentially be the true BC. ndt is typically 4. The
// total width of the time window is required to be at least 2*nMinBCs* LHCBunchSpacingNS

template <typename T>
T MCcompatibleBCs(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>::iterator const& collision, int ndt, T const& bcs, int nMinBCs = 7)
{
  LOGF(debug, "Collision time / resolution [ns]: %f / %f", collision.collisionTime(), collision.collisionTimeRes());

  // return if collisions has no associated BC
  if (!collision.has_foundBC()) {
    LOGF(info, "Collision %i - no BC found!", collision.globalIndex());
    return T{{bcs.asArrowTable()->Slice(0, 0)}, (uint64_t)0};
  }

  // get associated BC
  auto bcIter = collision.foundBC_as<T>();

  // due to the filling scheme the most probable BC may not be the one estimated from the collision time
  uint64_t mostProbableBC = bcIter.globalBC();
  uint64_t meanBC = mostProbableBC + std::lround(collision.collisionTime() / o2::constants::lhc::LHCBunchSpacingNS);

  // enforce minimum number for deltaBC
  int deltaBC = std::ceil(collision.collisionTimeRes() / o2::constants::lhc::LHCBunchSpacingNS * ndt);
  if (deltaBC < nMinBCs) {
    deltaBC = nMinBCs;
  }

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

// -----------------------------------------------------------------------------
// In PYTHIA a central diffractive produced (CD) particle has the ID
// 9900110. Check the particles of a MC event to contain a CD particle.
template <typename T>
bool isPythiaCDE(T MCparts)
{
  for (auto mcpart : MCparts) {
    if (mcpart.pdgCode() == 9900110) {
      return true;
    }
  }
  return false;
}

// -----------------------------------------------------------------------------
// In pp events produced with GRANIITTI the stack starts with
// 22212/22212/99/22212/2212/99/90
template <typename T>
bool isGraniittiCDE(T MCparts)
{
  if (MCparts.size() < 7) {
    return false;
  } else {
    if (MCparts.iteratorAt(0).pdgCode() != 2212)
      return false;
    if (MCparts.iteratorAt(1).pdgCode() != 2212)
      return false;
    if (MCparts.iteratorAt(2).pdgCode() != 99)
      return false;
    if (MCparts.iteratorAt(3).pdgCode() != 2212)
      return false;
    if (MCparts.iteratorAt(4).pdgCode() != 2212)
      return false;
    if (MCparts.iteratorAt(5).pdgCode() != 99)
      return false;
    if (MCparts.iteratorAt(6).pdgCode() != 90)
      return false;
  }

  return true;
}

// -----------------------------------------------------------------------------

#endif // O2_ANALYSIS_DGMCHELPER_H_
