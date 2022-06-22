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

#ifndef O2PHYSICS_UDTABLES_H
#define O2PHYSICS_UDTABLES_H

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/DataTypes.h"
#include "MathUtils/Utils.h"
#include <cmath>

namespace o2::aod
{

// only MCH-MID tracks
namespace skimmuontrack
{
DECLARE_SOA_COLUMN(Px, px, float);                     //!
DECLARE_SOA_COLUMN(Py, py, float);                     //!
DECLARE_SOA_COLUMN(Pz, pz, float);                     //!
DECLARE_SOA_COLUMN(Sign, sign, int);                   //!
DECLARE_SOA_COLUMN(TrackTime, trackTime, double);      //! absolute time in ns
DECLARE_SOA_COLUMN(TrackTimeRes, trackTimeRes, float); //! time resolution

DECLARE_SOA_DYNAMIC_COLUMN(MIDBoardCh1, midBoardCh1, //!
                           [](uint32_t midBoards) -> int { return static_cast<int>(midBoards & 0xFF); });
DECLARE_SOA_DYNAMIC_COLUMN(MIDBoardCh2, midBoardCh2, //!
                           [](uint32_t midBoards) -> int { return static_cast<int>((midBoards >> 8) & 0xFF); });
DECLARE_SOA_DYNAMIC_COLUMN(MIDBoardCh3, midBoardCh3, //!
                           [](uint32_t midBoards) -> int { return static_cast<int>((midBoards >> 16) & 0xFF); });
DECLARE_SOA_DYNAMIC_COLUMN(MIDBoardCh4, midBoardCh4, //!
                           [](uint32_t midBoards) -> int { return static_cast<int>((midBoards >> 24) & 0xFF); });
} // namespace skimmuontrack

// Muon track kinematics
DECLARE_SOA_TABLE(SkimmedMuons, "AOD", "SKIMMUONTRACK",
                  o2::soa::Index<>,
                  skimmuontrack::Px,
                  skimmuontrack::Py,
                  skimmuontrack::Pz,
                  skimmuontrack::Sign,
                  skimmuontrack::TrackTime,
                  skimmuontrack::TrackTimeRes);

// Muon track quality details
DECLARE_SOA_TABLE(SkimmedMuonsExtra, "AOD", "SKIMMUONEXTRA",
                  fwdtrack::NClusters,
                  fwdtrack::PDca,
                  fwdtrack::RAtAbsorberEnd,
                  fwdtrack::Chi2,
                  fwdtrack::Chi2MatchMCHMID,
                  fwdtrack::MCHBitMap,
                  fwdtrack::MIDBitMap,
                  fwdtrack::MIDBoards);

// Muon covariance
DECLARE_SOA_TABLE(SkimmedMuonsCov, "AOD", "SKIMMUONCOV",
                  fwdtrack::X,
                  fwdtrack::Y,
                  fwdtrack::Z,
                  fwdtrack::Tgl,
                  fwdtrack::Signed1Pt,
                  fwdtrack::CXX,
                  fwdtrack::CXY,
                  fwdtrack::CYY,
                  fwdtrack::CPhiX,
                  fwdtrack::CPhiY,
                  fwdtrack::CPhiPhi,
                  fwdtrack::CTglX,
                  fwdtrack::CTglY,
                  fwdtrack::CTglPhi,
                  fwdtrack::CTglTgl,
                  fwdtrack::C1PtX,
                  fwdtrack::C1PtY,
                  fwdtrack::C1PtPhi,
                  fwdtrack::C1PtTgl,
                  fwdtrack::C1Pt21Pt2);

using SkimmedMuon = SkimmedMuons::iterator;
using SkimmedMuonExtra = SkimmedMuonsExtra::iterator;
using SkimmedMuonCov = SkimmedMuonsCov::iterator;

DECLARE_SOA_TABLE(SkimmedMCEvents, "AOD", "SKMCEVENTS",
                  o2::soa::Index<>,
                  mccollision::GeneratorsID,
                  mccollision::PosX,
                  mccollision::PosY,
                  mccollision::PosZ,
                  mccollision::T,
                  mccollision::Weight,
                  mccollision::ImpactParameter);

namespace skimmcpart
{
DECLARE_SOA_INDEX_COLUMN(SkimmedMCEvent, skimmedMCEvent);  //!
DECLARE_SOA_SELF_ARRAY_INDEX_COLUMN(Mothers, mothers);     //! Mother tracks (possible empty) array. Iterate over mcParticle.mothers_as<aod::McParticles>())
DECLARE_SOA_SELF_SLICE_INDEX_COLUMN(Daughters, daughters); //! Daughter tracks (possibly empty) slice. Check for non-zero with mcParticle.has_daughters(). Iterate over mcParticle.daughters_as<aod::McParticles>())
DECLARE_SOA_COLUMN(Px, px, float);                         //!
DECLARE_SOA_COLUMN(Py, py, float);                         //!
DECLARE_SOA_COLUMN(Pz, pz, float);                         //!
DECLARE_SOA_COLUMN(E, e, float);                           //!
} // namespace skimmcpart

DECLARE_SOA_TABLE_FULL(SkimmedMCParticles, "SkimmedMCParticles", "AOD", "SKMCPARTICLES", //!  MC track information (on disk)
                       o2::soa::Index<>, skimmcpart::SkimmedMCEventId,
                       mcparticle::PdgCode,
                       mcparticle::StatusCode,
                       mcparticle::Flags,
                       skimmcpart::MothersIds,
                       skimmcpart::DaughtersIdSlice,
                       mcparticle::Weight,
                       skimmcpart::Px,
                       skimmcpart::Py,
                       skimmcpart::Pz,
                       skimmcpart::E,
                       mcparticle::ProducedByGenerator<mcparticle::Flags>,
                       mcparticle::FromBackgroundEvent<mcparticle::Flags>,
                       mcparticle::GetGenStatusCode<mcparticle::Flags, mcparticle::StatusCode>,
                       mcparticle::GetProcess<mcparticle::Flags, mcparticle::StatusCode>,
                       mcparticle::IsPhysicalPrimary<mcparticle::Flags>);

using SkimmedMCParticle = SkimmedMCParticles::iterator;

namespace skimmuontracklabel
{
DECLARE_SOA_INDEX_COLUMN(McParticle, mcParticle);
DECLARE_SOA_COLUMN(McMask, mcMask, uint16_t);
} // namespace skimmuontracklabel

DECLARE_SOA_TABLE(SkimmedMuonTrackLabels, "AOD", "SKMUONTRLABEL",
                  skimmuontracklabel::McParticleId,
                  skimmuontracklabel::McMask);

using SkimmedMuonTrackLabel = SkimmedMuonTrackLabels::iterator;

namespace eventcand
{
DECLARE_SOA_COLUMN(RunNumber, runNumber, int32_t);                             //! run number
DECLARE_SOA_COLUMN(GlobalBC, globalBC, uint64_t);                              //! global BC instead of BC ID since candidate may not have a corresponding record in BCs table
//
DECLARE_SOA_COLUMN(TotalFT0AmplitudeA, totalFT0AmplitudeA, float);             //! sum of amplitudes on A side of FT0
DECLARE_SOA_COLUMN(TotalFT0AmplitudeC, totalFT0AmplitudeC, float);             //! sum of amplitudes on C side of FT0
///DECLARE_SOA_DYNAMIC_COLUMN(HasFT0, hasFT0, //!
///                           [](uint32_t midBoards) -> int { return static_cast<int>(midBoards & 0xFF); });
//
DECLARE_SOA_SELF_ARRAY_INDEX_COLUMN(MatchedFwdTracks, matchedFwdTracks);       //! array of matched forward tracks
DECLARE_SOA_SELF_ARRAY_INDEX_COLUMN(MatchedBarrelTracks, matchedBarrelTracks); //! array of matched barrel tracks
} // namespace eventcand

DECLARE_SOA_TABLE(EventCandidates, "AOD", "EVENTCAND",
                  o2::soa::Index<>,
                  eventcand::GlobalBC,
                  eventcand::RunNumber,
                  eventcand::MatchedFwdTracksIds,
                  eventcand::MatchedBarrelTracksIds,
                  //
                  eventcand::TotalFT0AmplitudeA,
                  eventcand::TotalFT0AmplitudeC,
                  ft0::TimeA,
                  ft0::TimeC,
                  ft0::TriggerMask);

using EventCanditate = EventCandidates::iterator;

} // namespace o2::aod

#endif // O2PHYSICS_UDTABLES_H
