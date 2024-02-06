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
#ifndef PWGCF_TWOPARTICLECORRELATIONS_DATAMODEL_IDENTIFIEDBFFILTERED_H_
#define PWGCF_TWOPARTICLECORRELATIONS_DATAMODEL_IDENTIFIEDBFFILTERED_H_

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"

namespace o2
{
namespace aod
{
/* we have to change from int to bool when bool columns work properly */
namespace identifiedbffilter
{
DECLARE_SOA_COLUMN(IdentifiedBfCFCollisionAccepted, collisionaccepted, uint8_t); //! If the collision/event has been accepted or not
DECLARE_SOA_COLUMN(IdentifiedBfCFCollisionCentMult, centmult, float);            //! The centrality/multiplicity pecentile
DECLARE_SOA_DYNAMIC_COLUMN(IsCollisionAccepted,                                  //! Is the collision/event accepted
                           iscollisionaccepted,
                           [](uint8_t _collisionaccepted) -> uint8_t { return _collisionaccepted; });
DECLARE_SOA_DYNAMIC_COLUMN(IsGenCollisionAccepted, //! Is the generated collision/event accepted
                           isgencollisionaccepted,
                           [](uint8_t _collisionaccepted) -> uint8_t { return _collisionaccepted; });
} // namespace identifiedbffilter
DECLARE_SOA_TABLE(IdentifiedBfCFAcceptedCollisions,
                  "AOD",
                  "IBFCFACCCOLL", //! Accepted reconstructed collisions/events filtered table
                  o2::soa::Index<>,
                  collision::BCId,
                  collision::PosZ,
                  identifiedbffilter::IdentifiedBfCFCollisionAccepted,
                  identifiedbffilter::IdentifiedBfCFCollisionCentMult,
                  identifiedbffilter::IsCollisionAccepted<identifiedbffilter::IdentifiedBfCFCollisionAccepted>);
using IdentifiedBfCFAcceptedCollision = IdentifiedBfCFAcceptedCollisions::iterator;
DECLARE_SOA_TABLE(IdentifiedBfCFAcceptedTrueCollisions,
                  "AOD",
                  "IBFCFACCGENCOLL", //! Accepted generated collisions/events filtered table
                  o2::soa::Index<>,
                  collision::BCId,
                  mccollision::PosZ,
                  identifiedbffilter::IdentifiedBfCFCollisionAccepted,
                  identifiedbffilter::IdentifiedBfCFCollisionCentMult,
                  identifiedbffilter::IsGenCollisionAccepted<identifiedbffilter::IdentifiedBfCFCollisionAccepted>);
using IdentifiedBfCFAcceptedTrueCollision = IdentifiedBfCFAcceptedTrueCollisions::iterator;
DECLARE_SOA_TABLE(
  IdentifiedBfCFCollisionsInfo,
  "AOD",
  "IBFCFCOLLINF", //! Accepted reconstructed collisions/events Collisions joinable table
  identifiedbffilter::IdentifiedBfCFCollisionAccepted,
  identifiedbffilter::IdentifiedBfCFCollisionCentMult,
  identifiedbffilter::IsCollisionAccepted<identifiedbffilter::IdentifiedBfCFCollisionAccepted>);
DECLARE_SOA_TABLE(
  IdentifiedBfCFGenCollisionsInfo,
  "AOD",
  "IBFGENCFCOLLINF", //! Accepted generated collisions/events mcCollisions joinable table
  identifiedbffilter::IdentifiedBfCFCollisionAccepted,
  identifiedbffilter::IdentifiedBfCFCollisionCentMult,
  identifiedbffilter::IsGenCollisionAccepted<identifiedbffilter::IdentifiedBfCFCollisionAccepted>);
namespace identifiedbffilter
{
DECLARE_SOA_INDEX_COLUMN(IdentifiedBfCFAcceptedCollision, event);       //! Reconstructed collision/event
DECLARE_SOA_INDEX_COLUMN(IdentifiedBfCFAcceptedTrueCollision, mcevent); //! Generated collision/event
DECLARE_SOA_COLUMN(TrackacceptedId,
                   trackacceptedid,
                   int8_t);          //! Id of accepted track, criteria: even (+) odd (-)
DECLARE_SOA_COLUMN(Pt, pt, float);   //! The track transverse momentum
DECLARE_SOA_COLUMN(Eta, eta, float); //! The track pseudorapidity
DECLARE_SOA_COLUMN(Phi, phi, float); //! The track azimuthal angle
DECLARE_SOA_DYNAMIC_COLUMN(Sign,
                           sign, //! Charge: positive: 1, negative: -1
                           [](int8_t trackacceptedid) -> int8_t {
                             return (trackacceptedid % 2 == 0
                                       ? 1
                                       : (trackacceptedid % 2 == 1 ? -1 : 0));
                           });
DECLARE_SOA_DYNAMIC_COLUMN(TrkID,
                           trkid, //! The track id
                           [](int8_t trackacceptedid) -> int8_t { return trackacceptedid; });
DECLARE_SOA_DYNAMIC_COLUMN(PartID,
                           partid, //! The generated particle id
                           [](int8_t trackacceptedid) -> int8_t { return trackacceptedid; });
} // namespace identifiedbffilter
DECLARE_SOA_TABLE(ScannedTracks, "AOD", "SCANNEDTRACKS", //! The reconstructed tracks filtered table
                  identifiedbffilter::IdentifiedBfCFAcceptedCollisionId,
                  identifiedbffilter::TrackacceptedId,
                  identifiedbffilter::Pt,
                  identifiedbffilter::Eta,
                  identifiedbffilter::Phi,
                  identifiedbffilter::Sign<identifiedbffilter::TrackacceptedId>);
DECLARE_SOA_TABLE(ScannedTrueTracks, "AOD", "SCANTRUETRACKS", //! The generated particles filtered table
                  identifiedbffilter::IdentifiedBfCFAcceptedTrueCollisionId,
                  identifiedbffilter::TrackacceptedId,
                  identifiedbffilter::Pt,
                  identifiedbffilter::Eta,
                  identifiedbffilter::Phi,
                  identifiedbffilter::Sign<identifiedbffilter::TrackacceptedId>);
DECLARE_SOA_TABLE(IdentifiedBfCFTracksInfo, "AOD", "SCANDTRCKINF", //! The additional information Tracks joinable table
                  identifiedbffilter::TrackacceptedId,
                  identifiedbffilter::TrkID<identifiedbffilter::TrackacceptedId>);
DECLARE_SOA_TABLE(IdentifiedBfCFGenTracksInfo, "AOD", "SCANDGENTRCKINF", //! The additional information mcParticle joinable table
                  identifiedbffilter::TrackacceptedId,
                  identifiedbffilter::Sign<identifiedbffilter::TrackacceptedId>,
                  identifiedbffilter::PartID<identifiedbffilter::TrackacceptedId>);
} // namespace aod
} // namespace o2

#endif // PWGCF_TWOPARTICLECORRELATIONS_DATAMODEL_IDENTIFIEDBFFILTERED_H_
