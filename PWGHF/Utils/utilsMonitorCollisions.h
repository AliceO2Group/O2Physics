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

/// \file utilsMonitoCollisions.h
/// \brief Utility to monitor collisions passing the event selections for HF analyses
/// \author Mattia Faggin <mfaggin@cern.ch>, CERN

#ifndef PWGHF_UTILS_MONITORCOLLISIONS_H_
#define PWGHF_UTILS_MONITORCOLLISIONS_H_

/// @brief function to monitor the event selection satisfied by collisions used for HF analyses
/// \param collisions are the collisions in the DF
/// \param hCollisions is a histogram to keep track of the satisfied event selections
/// \param hPosZBeforeEvSel is PV position Z for all analysed collisions
/// \param hPosZAfterEvSel is PV position Z only for collisions satisfying the event selections
template <typename Coll, typename Hist>
void monitorCollision(Coll const& collision, const uint16_t statusCollision, Hist const& hCollisions, Hist const& hPosZBeforeEvSel, Hist const& hPosZAfterEvSel)
{

  hCollisions->Fill(0); // all collisions
  const float posZ = collision.posZ();
  hPosZBeforeEvSel->Fill(posZ);

  /// centrality
  if (TESTBIT(statusCollision, EventRejection::Centrality)) {
    return;
  }
  hCollisions->Fill(1); // Centrality ok

  /// sel8()
  if (TESTBIT(statusCollision, EventRejection::Trigger)) {
    return;
  }
  hCollisions->Fill(2); // Centrality + sel8 ok

  /// PV position Z
  if (TESTBIT(statusCollision, EventRejection::PositionZ)) {
    return;
  }
  hCollisions->Fill(3); // Centrality + sel8 + posZ ok

  /// Time Frame border cut
  if (TESTBIT(statusCollision, EventRejection::TimeFrameBorderCut)) {
    return;
  }
  hCollisions->Fill(4); // Centrality + sel8 + posZ + TF border ok
  hPosZAfterEvSel->Fill(posZ);
}

#endif