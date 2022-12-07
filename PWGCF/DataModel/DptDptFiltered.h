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
#ifndef O2_ANALYSIS_DPTDPTFILTERED_H
#define O2_ANALYSIS_DPTDPTFILTERED_H

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"

namespace o2
{
namespace aod
{
/* we have to change from int to bool when bool columns work properly */
namespace dptdptfilter
{
DECLARE_SOA_COLUMN(DptDptCFCollisionAccepted, collisionaccepted, uint8_t); //! If the collision/event has been accepted or not
DECLARE_SOA_COLUMN(DptDptCFCollisionCentMult, centmult, float);            //! The centrality/multiplicity pecentile
} // namespace dptdptfilter
DECLARE_SOA_TABLE(DptDptCFAcceptedCollisions, "AOD", "DPTDPTCFACCCOLL", //! Accepted reconstructed collisions/events filtered table
                  o2::soa::Index<>,
                  collision::BCId,
                  collision::PosZ,
                  dptdptfilter::DptDptCFCollisionAccepted,
                  dptdptfilter::DptDptCFCollisionCentMult);
using DptDptCFAcceptedCollision = DptDptCFAcceptedCollisions::iterator;
DECLARE_SOA_TABLE(DptDptCFAcceptedTrueCollisions, "AOD", "DPTCFACCGENCOLL", //! Accepted generated collisions/events filtered table
                  o2::soa::Index<>,
                  collision::BCId,
                  mccollision::PosZ,
                  dptdptfilter::DptDptCFCollisionAccepted,
                  dptdptfilter::DptDptCFCollisionCentMult);
using DptDptCFAcceptedTrueCollision = DptDptCFAcceptedTrueCollisions::iterator;
namespace dptdptfilter
{
DECLARE_SOA_INDEX_COLUMN(DptDptCFAcceptedCollision, event);          //! Reconstructed collision/event
DECLARE_SOA_INDEX_COLUMN(DptDptCFAcceptedTrueCollision, mcevent);    //! Generated collision/event
DECLARE_SOA_COLUMN(TrackacceptedAsOne, trackacceptedasone, uint8_t); //! Track accepted as type one
DECLARE_SOA_COLUMN(TrackacceptedAsTwo, trackacceptedastwo, uint8_t); //! Track accepted as type two
DECLARE_SOA_COLUMN(Pt, pt, float);                                   //! The track transverse momentum
DECLARE_SOA_COLUMN(Eta, eta, float);                                 //! The track pseudorapidity
DECLARE_SOA_COLUMN(Phi, phi, float);                                 //! The track azimuthal angle
DECLARE_SOA_DYNAMIC_COLUMN(Sign, sign,                               //! Charge: positive: 1, negative: -1
                           [](uint8_t trackacceptedasone, uint8_t trackacceptedastwo) -> short {
                             return (trackacceptedasone == uint8_t(true) ? 1 : (trackacceptedastwo ? -1 : 0));
                           });
} // namespace dptdptfilter
DECLARE_SOA_TABLE(ScannedTracks, "AOD", "SCANNEDTRACKS", //! The reconstructed tracks filtered table
                  dptdptfilter::DptDptCFAcceptedCollisionId,
                  dptdptfilter::TrackacceptedAsOne,
                  dptdptfilter::TrackacceptedAsTwo,
                  dptdptfilter::Pt,
                  dptdptfilter::Eta,
                  dptdptfilter::Phi,
                  dptdptfilter::Sign<dptdptfilter::TrackacceptedAsOne, dptdptfilter::TrackacceptedAsTwo>);
DECLARE_SOA_TABLE(ScannedTrueTracks, "AOD", "SCANTRUETRACKS", //! The generated particles filtered table
                  dptdptfilter::DptDptCFAcceptedTrueCollisionId,
                  dptdptfilter::TrackacceptedAsOne,
                  dptdptfilter::TrackacceptedAsTwo,
                  dptdptfilter::Pt,
                  dptdptfilter::Eta,
                  dptdptfilter::Phi,
                  dptdptfilter::Sign<dptdptfilter::TrackacceptedAsOne, dptdptfilter::TrackacceptedAsTwo>);
} // namespace aod
} // namespace o2

#endif // O2_ANALYSIS_DPTDPTFILTERED_H
