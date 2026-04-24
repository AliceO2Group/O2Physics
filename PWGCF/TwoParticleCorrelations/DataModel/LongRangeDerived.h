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
/// \file LongRangeDerived.h
///
/// \brief task derived table definition for long range correlation
/// \author Abhi Modak (abhi.modak@cern.ch)
/// \since October 28, 2025

#ifndef PWGCF_TWOPARTICLECORRELATIONS_DATAMODEL_LONGRANGEDERIVED_H_
#define PWGCF_TWOPARTICLECORRELATIONS_DATAMODEL_LONGRANGEDERIVED_H_

#include "Common/DataModel/Multiplicity.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>

#include <cstdint>

namespace o2::aod
{
namespace lrcorrmccolltable
{
DECLARE_SOA_COLUMN(Multiplicity, multiplicity, float);
} // namespace lrcorrmccolltable
DECLARE_SOA_TABLE(LRMcCollisions, "AOD", "LRMCCOLLISION",
                  o2::soa::Index<>,
                  mccollision::PosZ,
                  lrcorrmccolltable::Multiplicity,
                  mult::MultMCFT0A,
                  mult::MultMCFT0C);
using LRMcCollision = LRMcCollisions::iterator;

namespace lrcorrmctrktable
{
DECLARE_SOA_INDEX_COLUMN(LRMcCollision, lrMcCollision);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
} // namespace lrcorrmctrktable

DECLARE_SOA_TABLE(LRMidMcTracks, "AOD", "LRMIDMCTRACK",
                  o2::soa::Index<>,
                  lrcorrmctrktable::LRMcCollisionId,
                  lrcorrmctrktable::Pt,
                  lrcorrmctrktable::Eta,
                  lrcorrmctrktable::Phi,
                  mcparticle::PdgCode,
                  mcparticle::Flags,
                  mcparticle::IsPhysicalPrimary<mcparticle::Flags>);
using LRMidMcTrack = LRMidMcTracks::iterator;

DECLARE_SOA_TABLE(LRFt0aMcTracks, "AOD", "LRFT0AMCTRACK",
                  o2::soa::Index<>,
                  lrcorrmctrktable::LRMcCollisionId,
                  lrcorrmctrktable::Pt,
                  lrcorrmctrktable::Eta,
                  lrcorrmctrktable::Phi);
using LRFt0aMcTrack = LRFt0aMcTracks::iterator;

DECLARE_SOA_TABLE(LRFt0cMcTracks, "AOD", "LRFT0CMCTRACK",
                  o2::soa::Index<>,
                  lrcorrmctrktable::LRMcCollisionId,
                  lrcorrmctrktable::Pt,
                  lrcorrmctrktable::Eta,
                  lrcorrmctrktable::Phi);
using LRFt0cMcTrack = LRFt0cMcTracks::iterator;

DECLARE_SOA_TABLE(LRMftMcTracks, "AOD", "LRMFTMCTRACK",
                  o2::soa::Index<>,
                  lrcorrmctrktable::LRMcCollisionId,
                  lrcorrmctrktable::Pt,
                  lrcorrmctrktable::Eta,
                  lrcorrmctrktable::Phi);
using LRMftMcTrack = LRMftMcTracks::iterator;

namespace lrcorrcolltable
{
DECLARE_SOA_INDEX_COLUMN(LRMcCollision, lrMcCollision);
DECLARE_SOA_COLUMN(Multiplicity, multiplicity, float);
DECLARE_SOA_COLUMN(Centrality, centrality, float);
DECLARE_SOA_COLUMN(TotalFT0AmplitudeA, totalFT0AmplitudeA, float); //! sum of amplitudes on A side of FT0
DECLARE_SOA_COLUMN(TotalFT0AmplitudeC, totalFT0AmplitudeC, float); //! sum of amplitudes on C side of FT0
DECLARE_SOA_COLUMN(TotalFV0AmplitudeA, totalFV0AmplitudeA, float); //! sum of amplitudes on A side of FDD
DECLARE_SOA_COLUMN(GapSide, gapSide, uint8_t);                     // 0 for side A, 1 for side C, 2 for both sides
} // namespace lrcorrcolltable

DECLARE_SOA_TABLE(LRCollisions, "AOD", "LRCOLLISION",
                  o2::soa::Index<>,
                  bc::RunNumber,
                  collision::PosZ,
                  lrcorrcolltable::Multiplicity,
                  lrcorrcolltable::Centrality,
                  timestamp::Timestamp);
DECLARE_SOA_TABLE(LRCollLabels, "AOD", "LRCOLLLABEL",
                  lrcorrcolltable::LRMcCollisionId);
using LRCollision = LRCollisions::iterator;
using LRCollLabel = LRCollLabels::iterator;
using LRCollisionsWithLabel = soa::Join<LRCollisions, LRCollLabels>;
using LRCollisionWithLabel = LRCollisionsWithLabel::iterator;

DECLARE_SOA_TABLE(UpcLRCollisions, "AOD", "UPCLRCOLLISION",
                  o2::soa::Index<>,
                  bc::GlobalBC,
                  bc::RunNumber,
                  collision::PosZ,
                  lrcorrcolltable::Multiplicity,
                  lrcorrcolltable::TotalFT0AmplitudeA,
                  lrcorrcolltable::TotalFT0AmplitudeC,
                  lrcorrcolltable::TotalFV0AmplitudeA);
using UpcLRCollision = UpcLRCollisions::iterator;

DECLARE_SOA_TABLE(UpcSgLRCollisions, "AOD", "UPCSGLRCOLLISION",
                  lrcorrcolltable::GapSide);
using UpcSgLRCollision = UpcSgLRCollisions::iterator;

namespace lrcorrzdctable
{
DECLARE_SOA_INDEX_COLUMN(UpcLRCollision, upcLRCollision);
DECLARE_SOA_COLUMN(EnergyCommonZNA, energyCommonZNA, float);
DECLARE_SOA_COLUMN(EnergyCommonZNC, energyCommonZNC, float);
} // namespace lrcorrzdctable

DECLARE_SOA_TABLE(LRZdcs, "AOD", "LRZDC",
                  o2::soa::Index<>,
                  lrcorrzdctable::UpcLRCollisionId,
                  lrcorrzdctable::EnergyCommonZNA,
                  lrcorrzdctable::EnergyCommonZNC);
using LRZdc = LRZdcs::iterator;

namespace lrcorrtrktable
{
DECLARE_SOA_INDEX_COLUMN(LRCollision, lrCollision);
DECLARE_SOA_INDEX_COLUMN(UpcLRCollision, upcLRCollision);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(ChannelID, channelID, int);
DECLARE_SOA_COLUMN(Amplitude, amplitude, float);
DECLARE_SOA_COLUMN(GainAmplitude, gainAmplitude, float);
DECLARE_SOA_COLUMN(InvMass, invMass, float);
DECLARE_SOA_COLUMN(IdPos, idPos, int64_t);
DECLARE_SOA_COLUMN(IdNeg, idNeg, int64_t);
DECLARE_SOA_COLUMN(TrackType, trackType, uint8_t);
DECLARE_SOA_COLUMN(V0Type, v0Type, uint8_t);
enum TrackPid {
  kSpCharge,
  kSpPion,
  kSpKaon,
  kSpProton,
  kNoPid
};
enum V0TrackPid {
  kSpK0short,
  kSpLambda,
  kSpALambda
};
} // namespace lrcorrtrktable

DECLARE_SOA_TABLE(LRMidTracks, "AOD", "LRMIDTRACK",
                  o2::soa::Index<>,
                  lrcorrtrktable::LRCollisionId,
                  lrcorrtrktable::Pt,
                  lrcorrtrktable::Eta,
                  lrcorrtrktable::Phi,
                  lrcorrtrktable::TrackType);
using LRMidTrack = LRMidTracks::iterator;

DECLARE_SOA_TABLE(LRFt0aTracks, "AOD", "LRFT0ATRACK",
                  o2::soa::Index<>,
                  lrcorrtrktable::LRCollisionId,
                  lrcorrtrktable::ChannelID,
                  lrcorrtrktable::Amplitude,
                  lrcorrtrktable::GainAmplitude,
                  lrcorrtrktable::Eta,
                  lrcorrtrktable::Phi);
using LRFt0aTrack = LRFt0aTracks::iterator;

DECLARE_SOA_TABLE(LRFt0cTracks, "AOD", "LRFT0CTRACK",
                  o2::soa::Index<>,
                  lrcorrtrktable::LRCollisionId,
                  lrcorrtrktable::ChannelID,
                  lrcorrtrktable::Amplitude,
                  lrcorrtrktable::GainAmplitude,
                  lrcorrtrktable::Eta,
                  lrcorrtrktable::Phi);
using LRFt0cTrack = LRFt0cTracks::iterator;

DECLARE_SOA_TABLE(LRV0Tracks, "AOD", "LRV0TRACK",
                  o2::soa::Index<>,
                  lrcorrtrktable::LRCollisionId,
                  lrcorrtrktable::IdPos,
                  lrcorrtrktable::IdNeg,
                  lrcorrtrktable::Pt,
                  lrcorrtrktable::Eta,
                  lrcorrtrktable::Phi,
                  lrcorrtrktable::InvMass,
                  lrcorrtrktable::V0Type);
using LRV0Track = LRV0Tracks::iterator;

DECLARE_SOA_TABLE(LRMftTracks, "AOD", "LRMFTTRACK",
                  o2::soa::Index<>,
                  lrcorrtrktable::LRCollisionId,
                  lrcorrtrktable::Pt,
                  lrcorrtrktable::Eta,
                  lrcorrtrktable::Phi);
using LRMftTrack = LRMftTracks::iterator;

DECLARE_SOA_TABLE(LRMftBestTracks, "AOD", "LRMFTBESTTRACK",
                  o2::soa::Index<>,
                  lrcorrtrktable::LRCollisionId,
                  lrcorrtrktable::Pt,
                  lrcorrtrktable::Eta,
                  lrcorrtrktable::Phi);
using LRMftBestTrack = LRMftBestTracks::iterator;

DECLARE_SOA_TABLE(UpcLRMidTracks, "AOD", "UPCLRMIDTRACK",
                  o2::soa::Index<>,
                  lrcorrtrktable::UpcLRCollisionId,
                  lrcorrtrktable::Pt,
                  lrcorrtrktable::Eta,
                  lrcorrtrktable::Phi,
                  lrcorrtrktable::TrackType);
using UpcLRMidTrack = UpcLRMidTracks::iterator;

DECLARE_SOA_TABLE(UpcLRFt0aTracks, "AOD", "UPCLRFT0ATRACK",
                  o2::soa::Index<>,
                  lrcorrtrktable::UpcLRCollisionId,
                  lrcorrtrktable::ChannelID,
                  lrcorrtrktable::Amplitude,
                  lrcorrtrktable::GainAmplitude,
                  lrcorrtrktable::Eta,
                  lrcorrtrktable::Phi);
using UpcLRFt0aTrack = UpcLRFt0aTracks::iterator;

DECLARE_SOA_TABLE(UpcLRFt0cTracks, "AOD", "UPCLRFT0CTRACK",
                  o2::soa::Index<>,
                  lrcorrtrktable::UpcLRCollisionId,
                  lrcorrtrktable::ChannelID,
                  lrcorrtrktable::Amplitude,
                  lrcorrtrktable::GainAmplitude,
                  lrcorrtrktable::Eta,
                  lrcorrtrktable::Phi);
using UpcLRFt0cTrack = UpcLRFt0cTracks::iterator;

DECLARE_SOA_TABLE(UpcLRV0Tracks, "AOD", "UPCLRV0TRACK",
                  o2::soa::Index<>,
                  lrcorrtrktable::UpcLRCollisionId,
                  lrcorrtrktable::IdPos,
                  lrcorrtrktable::IdNeg,
                  lrcorrtrktable::Pt,
                  lrcorrtrktable::Eta,
                  lrcorrtrktable::Phi,
                  lrcorrtrktable::InvMass,
                  lrcorrtrktable::V0Type);
using UpcLRV0Track = UpcLRV0Tracks::iterator;

DECLARE_SOA_TABLE(UpcLRMftTracks, "AOD", "UPCLRMFTTRACK",
                  o2::soa::Index<>,
                  lrcorrtrktable::UpcLRCollisionId,
                  lrcorrtrktable::Pt,
                  lrcorrtrktable::Eta,
                  lrcorrtrktable::Phi);
using UpcLRMftTrack = UpcLRMftTracks::iterator;

DECLARE_SOA_TABLE(UpcLRMftBestTracks, "AOD", "UPCLRMFTBESTTRACK",
                  o2::soa::Index<>,
                  lrcorrtrktable::UpcLRCollisionId,
                  lrcorrtrktable::Pt,
                  lrcorrtrktable::Eta,
                  lrcorrtrktable::Phi);
using UpcLRMftBestTrack = UpcLRMftBestTracks::iterator;
} // namespace o2::aod

#endif // PWGCF_TWOPARTICLECORRELATIONS_DATAMODEL_LONGRANGEDERIVED_H_
