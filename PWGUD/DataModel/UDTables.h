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

namespace skimbartrack
{
DECLARE_SOA_COLUMN(Px, px, float);                     //!
DECLARE_SOA_COLUMN(Py, py, float);                     //!
DECLARE_SOA_COLUMN(Pz, pz, float);                     //!
DECLARE_SOA_COLUMN(Sign, sign, int);                   //!
DECLARE_SOA_COLUMN(TrackTime, trackTime, double);      //! absolute time in ns
DECLARE_SOA_COLUMN(TrackTimeRes, trackTimeRes, float); //! time resolution
//
DECLARE_SOA_COLUMN(Flags, flags, uint32_t);                                                   //! Track flags. Run 2: see TrackFlagsRun2Enum | Run 3: see TrackFlags
DECLARE_SOA_COLUMN(ITSClusterMap, itsClusterMap, uint8_t);                                    //! ITS cluster map, one bit per a layer, starting from the innermost
DECLARE_SOA_COLUMN(TPCNClsFindable, tpcNClsFindable, uint8_t);                                //! Findable TPC clusters for this track geometry
DECLARE_SOA_COLUMN(TPCNClsFindableMinusFound, tpcNClsFindableMinusFound, int8_t);             //! TPC Clusters: Findable - Found
DECLARE_SOA_COLUMN(TPCNClsFindableMinusCrossedRows, tpcNClsFindableMinusCrossedRows, int8_t); //! TPC Clusters: Findable - crossed rows
DECLARE_SOA_COLUMN(TPCNClsShared, tpcNClsShared, uint8_t);                                    //! Number of shared TPC clusters
DECLARE_SOA_COLUMN(ITSChi2NCl, itsChi2NCl, float);                                            //! Chi2 / cluster for the ITS track segment
DECLARE_SOA_COLUMN(TPCChi2NCl, tpcChi2NCl, float);                                            //! Chi2 / cluster for the TPC track segment
DECLARE_SOA_COLUMN(TOFChi2, tofChi2, float);                                                  //! Chi2 for the TOF track segment
DECLARE_SOA_COLUMN(TPCSignal, tpcSignal, float);                                              //! dE/dx signal in the TPC
DECLARE_SOA_COLUMN(Length, length, float);                                                    //! Track length
DECLARE_SOA_COLUMN(TOFExpMom, tofExpMom, float);                                              //! TOF expected momentum obtained in tracking, used to compute the expected times
} // namespace skimbartrack

// Barrel track kinematics
DECLARE_SOA_TABLE(SkimmedBarTracks, "AOD", "SKIMBARTRACK",
                  o2::soa::Index<>,
                  skimbartrack::Px,
                  skimbartrack::Py,
                  skimbartrack::Pz,
                  skimbartrack::Sign,
                  skimbartrack::TrackTime,
                  skimbartrack::TrackTimeRes);

DECLARE_SOA_TABLE(SkimmedBarTracksCov, "AOD", "SKIMBARTRCOV", //!
                  track::X, track::Alpha,
                  track::Y, track::Z, track::Snp, track::Tgl, track::Signed1Pt,
                  track::CYY, track::CZY, track::CZZ, track::CSnpY, track::CSnpZ,
                  track::CSnpSnp, track::CTglY, track::CTglZ, track::CTglSnp, track::CTglTgl,
                  track::C1PtY, track::C1PtZ, track::C1PtSnp, track::C1PtTgl, track::C1Pt21Pt2);

DECLARE_SOA_TABLE(SkimmedBarTracksExtra, "AOD", "SKIMBARTREXTRA",
                  skimbartrack::Flags,
                  skimbartrack::ITSClusterMap,
                  skimbartrack::TPCNClsFindable,
                  skimbartrack::TPCNClsFindableMinusFound,
                  skimbartrack::TPCNClsFindableMinusCrossedRows,
                  skimbartrack::TPCNClsShared,
                  skimbartrack::ITSChi2NCl,
                  skimbartrack::TPCChi2NCl,
                  skimbartrack::TOFChi2,
                  skimbartrack::TPCSignal,
                  skimbartrack::Length,
                  skimbartrack::TOFExpMom);

using SkimmedBarTrack = SkimmedBarTracks::iterator;
using SkimmedBarTrackCov = SkimmedBarTracksCov::iterator;
using SkimmedBarTrackExtra = SkimmedBarTracksExtra::iterator;

namespace skimbartracklabel
{
DECLARE_SOA_INDEX_COLUMN(McParticle, mcParticle);
DECLARE_SOA_COLUMN(McMask, mcMask, uint16_t);
} // namespace skimbartracklabel

DECLARE_SOA_TABLE(SkimmedBarTrackLabels, "AOD", "SKBARTRLABEL",
                  skimbartracklabel::McParticleId,
                  skimbartracklabel::McMask);

using SkimmedBarTrackLabel = SkimmedBarTrackLabels::iterator;

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

namespace skimmcevent
{
DECLARE_SOA_COLUMN(GlobalBC, globalBC, uint64_t);
}

DECLARE_SOA_TABLE(SkimmedMCEvents, "AOD", "SKMCEVENTS",
                  o2::soa::Index<>,
                  skimmcevent::GlobalBC,
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
// general information
DECLARE_SOA_COLUMN(RunNumber, runNumber, int32_t); //! run number
DECLARE_SOA_COLUMN(GlobalBC, globalBC, uint64_t);  //! global BC instead of BC ID since candidate may not have a corresponding record in BCs table
// FT0 information
DECLARE_SOA_COLUMN(TotalAmplitudeAFT0, totalAmplitudeAFT0, float); //! sum of amplitudes on A side of FT0
DECLARE_SOA_COLUMN(TotalAmplitudeCFT0, totalAmplitudeCFT0, float); //! sum of amplitudes on C side of FT0
DECLARE_SOA_COLUMN(TimeAFT0, timeAFT0, float);                     //! FT0A average time
DECLARE_SOA_COLUMN(TimeCFT0, timeCFT0, float);                     //! FT0C average time
DECLARE_SOA_COLUMN(TriggerMaskFT0, triggerMaskFT0, uint8_t);       //! FT0 trigger mask
DECLARE_SOA_DYNAMIC_COLUMN(HasFT0, hasFT0,                         //! has FT0 signal in the same BC
                           [](float TimeAFT0, float TimeCFT0) -> bool { return TimeAFT0 > -999. && TimeCFT0 > -999.; });
// matched IDs
DECLARE_SOA_SELF_ARRAY_INDEX_COLUMN(MatchedFwdTracks, matchedFwdTracks); //! array of matched forward tracks
DECLARE_SOA_SELF_ARRAY_INDEX_COLUMN(MatchedBarTracks, matchedBarTracks); //! array of matched barrel tracks
} // namespace eventcand

DECLARE_SOA_TABLE(EventCandidates, "AOD", "EVENTCAND",
                  o2::soa::Index<>,
                  eventcand::GlobalBC,
                  eventcand::RunNumber,
                  eventcand::MatchedFwdTracksIds,
                  eventcand::MatchedBarTracksIds,
                  //
                  eventcand::TotalAmplitudeAFT0,
                  eventcand::TotalAmplitudeCFT0,
                  eventcand::TimeAFT0,
                  eventcand::TimeCFT0,
                  eventcand::TriggerMaskFT0,
                  eventcand::HasFT0<eventcand::TimeAFT0, eventcand::TimeCFT0>);

using EventCanditate = EventCandidates::iterator;

} // namespace o2::aod

#endif // O2PHYSICS_UDTABLES_H
