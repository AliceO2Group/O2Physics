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
//
// Contact: iarsene@cern.ch, i.c.arsene@fys.uio.no
//

#ifndef O2_Analysis_ReducedInfoTables_H_
#define O2_Analysis_ReducedInfoTables_H_

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "MathUtils/Utils.h"
#include <cmath>

namespace o2::aod
{

namespace dqppfilter
{
DECLARE_SOA_COLUMN(EventFilter, eventFilter, uint64_t); //! Bit-field used for the high level event triggering
}

DECLARE_SOA_TABLE(DQEventFilter, "AOD", "EVENTFILTER", //! Store event-level decisions (DQ high level triggers)
                  dqppfilter::EventFilter);

namespace reducedevent
{

// basic event information
DECLARE_SOA_COLUMN(Tag, tag, uint64_t);                   //!  Bit-field for storing event information (e.g. high level info, cut decisions)
DECLARE_SOA_COLUMN(TriggerAlias, triggerAlias, uint32_t); //!  Trigger aliases bit field
DECLARE_SOA_COLUMN(Q2X0A, q2x0a, float);                  //!  Q-vector x component, with event eta gap A (harmonic 2 and power 0)
DECLARE_SOA_COLUMN(Q2Y0A, q2y0a, float);                  //!  Q-vector y component, with event eta gap A (harmonic 2 and power 0)
DECLARE_SOA_COLUMN(Q2X0B, q2x0b, float);                  //!  Q-vector x component, with event eta gap B (harmonic 2 and power 0)
DECLARE_SOA_COLUMN(Q2Y0B, q2y0b, float);                  //!  Q-vector y component, with event eta gap B (harmonic 2 and power 0)
DECLARE_SOA_COLUMN(Q2X0C, q2x0c, float);                  //!  Q-vector x component, with event eta gap C (harmonic 2 and power 0)
DECLARE_SOA_COLUMN(Q2Y0C, q2y0c, float);                  //!  Q-vector y component, with event eta gap C (harmonic 2 and power 0)
DECLARE_SOA_COLUMN(MultA, multa, float);                  //!  Event multiplicity eta gap A
DECLARE_SOA_COLUMN(MultB, multb, float);                  //!  Event multiplicity eta gap B
DECLARE_SOA_COLUMN(MultC, multc, float);                  //!  Event multiplicity eta gap C
DECLARE_SOA_COLUMN(Q3X0A, q3x0a, float);                  //!  Q-vector x component, with event eta gap A (harmonic 3 and power 0)
DECLARE_SOA_COLUMN(Q3Y0A, q3y0a, float);                  //!
DECLARE_SOA_COLUMN(Q3X0B, q3x0b, float);                  //!
DECLARE_SOA_COLUMN(Q3Y0B, q3y0b, float);                  //!
DECLARE_SOA_COLUMN(Q3X0C, q3x0c, float);                  //!
DECLARE_SOA_COLUMN(Q3Y0C, q3y0c, float);                  //!
DECLARE_SOA_COLUMN(MCPosX, mcPosX, float);                //!
DECLARE_SOA_COLUMN(MCPosY, mcPosY, float);                //!
DECLARE_SOA_COLUMN(MCPosZ, mcPosZ, float);                //!
} // namespace reducedevent

DECLARE_SOA_TABLE(ReducedEvents, "AOD", "REDUCEDEVENT", //!   Main event information table
                  o2::soa::Index<>,
                  reducedevent::Tag, bc::RunNumber,
                  collision::PosX, collision::PosY, collision::PosZ, collision::NumContrib,
                  collision::CollisionTime, collision::CollisionTimeRes);

DECLARE_SOA_TABLE(ReducedEventsExtended, "AOD", "REEXTENDED", //!  Extended event information
                  bc::GlobalBC, bc::TriggerMask, timestamp::Timestamp, reducedevent::TriggerAlias, cent::CentRun2V0M);

DECLARE_SOA_TABLE(ReducedEventsVtxCov, "AOD", "REVTXCOV", //!    Event vertex covariance matrix
                  collision::CovXX, collision::CovXY, collision::CovXZ,
                  collision::CovYY, collision::CovYZ, collision::CovZZ, collision::Chi2);

DECLARE_SOA_TABLE(ReducedEventsQvector, "AOD", "REQVECTOR", //!    Event Q-vector information
                  reducedevent::Q2X0A, reducedevent::Q2Y0A, reducedevent::Q2X0B, reducedevent::Q2Y0B,
                  reducedevent::Q2X0C, reducedevent::Q2Y0C, reducedevent::MultA, reducedevent::MultB, reducedevent::MultC,
                  reducedevent::Q3X0A, reducedevent::Q3Y0A, reducedevent::Q3X0B, reducedevent::Q3Y0B,
                  reducedevent::Q3X0C, reducedevent::Q3Y0C);

// TODO and NOTE: This table is just an extension of the ReducedEvents table
//       There is no explicit accounting for MC events which were not reconstructed!!!
//       However, for analysis which will require these events, a special skimming process function
//           can be constructed and the same data model could be used
DECLARE_SOA_TABLE(ReducedMCEvents, "AOD", "REMC", //!   Event level MC truth information
                  o2::soa::Index<>,
                  mccollision::GeneratorsID, reducedevent::MCPosX, reducedevent::MCPosY, reducedevent::MCPosZ,
                  mccollision::T, mccollision::Weight, mccollision::ImpactParameter);

using ReducedEvent = ReducedEvents::iterator;
using ReducedEventExtended = ReducedEventsExtended::iterator;
using ReducedEventVtxCov = ReducedEventsVtxCov::iterator;
using ReducedEventQvector = ReducedEventsQvector::iterator;
using ReducedMCEvent = ReducedMCEvents::iterator;

namespace reducedeventlabel
{
DECLARE_SOA_INDEX_COLUMN(ReducedMCEvent, reducedMCevent); //! MC collision
DECLARE_SOA_COLUMN(McMask, mcMask, uint16_t);             //! Bit mask to indicate collision mismatches (bit ON means mismatch). Bit 15: indicates negative label
} // namespace reducedeventlabel

DECLARE_SOA_TABLE(ReducedMCEventLabels, "AOD", "REMCCOLLBL", //! Table joined to the ReducedEvents table containing the MC index
                  reducedeventlabel::ReducedMCEventId, reducedeventlabel::McMask);
using ReducedMCEventLabel = ReducedMCEventLabels::iterator;

namespace reducedtrack
{
// basic track information
DECLARE_SOA_INDEX_COLUMN(ReducedEvent, reducedevent); //!
// ----  flags reserved for storing various information during filtering
DECLARE_SOA_COLUMN(FilteringFlags, filteringFlags, uint64_t); //!
// -----------------------------------------------------
DECLARE_SOA_COLUMN(Pt, pt, float);       //!
DECLARE_SOA_COLUMN(Eta, eta, float);     //!
DECLARE_SOA_COLUMN(Phi, phi, float);     //!
DECLARE_SOA_COLUMN(Sign, sign, int);     //!
DECLARE_SOA_COLUMN(IsAmbiguous, isAmbiguous, int); //!
DECLARE_SOA_COLUMN(DcaXY, dcaXY, float); //!
DECLARE_SOA_COLUMN(DcaZ, dcaZ, float);   //!
DECLARE_SOA_DYNAMIC_COLUMN(Px, px,       //!
                           [](float pt, float phi) -> float { return pt * std::cos(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py, //!
                           [](float pt, float phi) -> float { return pt * std::sin(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz, //!
                           [](float pt, float eta) -> float { return pt * std::sinh(eta); });
DECLARE_SOA_DYNAMIC_COLUMN(P, p, //!
                           [](float pt, float eta) -> float { return pt * std::cosh(eta); });
} //namespace reducedtrack

// basic track information
DECLARE_SOA_TABLE(ReducedTracks, "AOD", "REDUCEDTRACK", //!
                  o2::soa::Index<>, reducedtrack::ReducedEventId, reducedtrack::FilteringFlags,
                  reducedtrack::Pt, reducedtrack::Eta, reducedtrack::Phi, reducedtrack::Sign, reducedtrack::IsAmbiguous,
                  reducedtrack::Px<reducedtrack::Pt, reducedtrack::Phi>,
                  reducedtrack::Py<reducedtrack::Pt, reducedtrack::Phi>,
                  reducedtrack::Pz<reducedtrack::Pt, reducedtrack::Eta>,
                  reducedtrack::P<reducedtrack::Pt, reducedtrack::Eta>);

// barrel track information
DECLARE_SOA_TABLE(ReducedTracksBarrel, "AOD", "RTBARREL", //!
                  track::TPCInnerParam, track::Flags,     // tracking status flags
                  track::ITSClusterMap, track::ITSChi2NCl,
                  track::TPCNClsFindable, track::TPCNClsFindableMinusFound, track::TPCNClsFindableMinusCrossedRows,
                  track::TPCNClsShared, track::TPCChi2NCl,
                  track::TRDChi2, track::TRDPattern, track::TOFChi2, track::Length, reducedtrack::DcaXY, reducedtrack::DcaZ,
                  track::TPCNClsFound<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  track::TPCNClsCrossedRows<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>);

// barrel covariance matrix  TODO: add all the elements required for secondary vertexing
DECLARE_SOA_TABLE(ReducedTracksBarrelCov, "AOD", "RTBARRELCOV", //!
                  track::X, track::Alpha,
                  track::Y, track::Z, track::Snp, track::Tgl, track::Signed1Pt,
                  track::CYY, track::CZY, track::CZZ, track::CSnpY, track::CSnpZ,
                  track::CSnpSnp, track::CTglY, track::CTglZ, track::CTglSnp, track::CTglTgl,
                  track::C1PtY, track::C1PtZ, track::C1PtSnp, track::C1PtTgl, track::C1Pt21Pt2);

// barrel PID information
DECLARE_SOA_TABLE(ReducedTracksBarrelPID, "AOD", "RTBARRELPID", //!
                  track::TPCSignal,
                  pidtpc::TPCNSigmaEl, pidtpc::TPCNSigmaMu,
                  pidtpc::TPCNSigmaPi, pidtpc::TPCNSigmaKa, pidtpc::TPCNSigmaPr,
                  pidtofbeta::Beta,
                  pidtof::TOFNSigmaEl, pidtof::TOFNSigmaMu,
                  pidtof::TOFNSigmaPi, pidtof::TOFNSigmaKa, pidtof::TOFNSigmaPr,
                  track::TRDSignal);

using ReducedTrack = ReducedTracks::iterator;
using ReducedTrackBarrel = ReducedTracksBarrel::iterator;
using ReducedTrackBarrelCov = ReducedTracksBarrelCov::iterator;
using ReducedTrackBarrelPID = ReducedTracksBarrelPID::iterator;

namespace reducedtrackMC
{
DECLARE_SOA_INDEX_COLUMN(ReducedMCEvent, reducedMCevent);                                   //!
DECLARE_SOA_COLUMN(McReducedFlags, mcReducedFlags, uint16_t);                               //! Flags to hold compressed MC selection information
DECLARE_SOA_SELF_INDEX_COLUMN_FULL(Mother0, mother0, int, "ReducedMCTracks_Mother0");       //! Track index of the first mother
DECLARE_SOA_SELF_INDEX_COLUMN_FULL(Mother1, mother1, int, "ReducedMCTracks_Mother1");       //! Track index of the last mother
DECLARE_SOA_SELF_INDEX_COLUMN_FULL(Daughter0, daughter0, int, "ReducedMCTracks_Daughter0"); //! Track index of the first daughter
DECLARE_SOA_SELF_INDEX_COLUMN_FULL(Daughter1, daughter1, int, "ReducedMCTracks_Daughter1"); //! Track index of the last daughter
DECLARE_SOA_SELF_ARRAY_INDEX_COLUMN(Mothers, mothers);                                      //! Mother tracks (possible empty) array. Iterate over mcParticle.mothers_as<aod::McParticles>())
DECLARE_SOA_SELF_SLICE_INDEX_COLUMN(Daughters, daughters);                                  //! Daughter tracks (possibly empty) slice. Check for non-zero with mcParticle.has_daughters(). Iterate over mcParticle.daughters_as<aod::McParticles>())
DECLARE_SOA_COLUMN(Pt, pt, float);                                                          //!
DECLARE_SOA_COLUMN(Eta, eta, float);                                                        //!
DECLARE_SOA_COLUMN(Phi, phi, float);                                                        //!
DECLARE_SOA_COLUMN(E, e, float);                                                            //!
DECLARE_SOA_DYNAMIC_COLUMN(Px, px,                                                          //!
                           [](float pt, float phi) -> float { return pt * std::cos(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py, //!
                           [](float pt, float phi) -> float { return pt * std::sin(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz, //!
                           [](float pt, float eta) -> float { return pt * std::sinh(eta); });
DECLARE_SOA_DYNAMIC_COLUMN(P, p, //!
                           [](float pt, float eta) -> float { return pt * std::cosh(eta); });
DECLARE_SOA_DYNAMIC_COLUMN(Y, y, //! Particle rapidity
                           [](float pt, float eta, float e) -> float {
                             float pz = pt * std::sinh(eta);
                             if ((e - pz) > static_cast<float>(1e-7)) {
                               return 0.5f * std::log((e + pz) / (e - pz));
                             } else {
                               return -999.0f;
                             }
                           });
} // namespace reducedtrackMC
// NOTE: This table is nearly identical to the one from Framework (except that it points to the event ID, not the BC id)
//       This table contains all MC truth tracks (both barrel and muon)
DECLARE_SOA_TABLE_FULL(ReducedMCTracks, "ReducedMCTracks", "AOD", "RTMC", //!  MC track information (on disk)
                       o2::soa::Index<>, reducedtrackMC::ReducedMCEventId,
                       mcparticle::PdgCode, mcparticle::StatusCode, mcparticle::Flags,
                       reducedtrackMC::MothersIds, reducedtrackMC::DaughtersIdSlice,
                       mcparticle::Weight,
                       reducedtrackMC::Pt, reducedtrackMC::Eta, reducedtrackMC::Phi, reducedtrackMC::E,
                       mcparticle::Vx, mcparticle::Vy, mcparticle::Vz, mcparticle::Vt,
                       reducedtrackMC::McReducedFlags,
                       reducedtrackMC::Px<reducedtrackMC::Pt, reducedtrackMC::Phi>,
                       reducedtrackMC::Py<reducedtrackMC::Pt, reducedtrackMC::Phi>,
                       reducedtrackMC::Pz<reducedtrackMC::Pt, reducedtrackMC::Eta>,
                       reducedtrackMC::P<reducedtrackMC::Pt, reducedtrackMC::Eta>,
                       reducedtrackMC::Y<reducedtrackMC::Pt, reducedtrackMC::Eta, reducedtrackMC::E>,
                       mcparticle::ProducedByGenerator<mcparticle::Flags>,
                       mcparticle::FromBackgroundEvent<mcparticle::Flags>,
                       mcparticle::GetGenStatusCode<mcparticle::Flags, mcparticle::StatusCode>,
                       mcparticle::GetProcess<mcparticle::Flags, mcparticle::StatusCode>,
                       mcparticle::IsPhysicalPrimary<mcparticle::Flags>);

using ReducedMCTrack = ReducedMCTracks::iterator;

namespace reducedbarreltracklabel
{
DECLARE_SOA_INDEX_COLUMN(ReducedMCTrack, reducedMCTrack); //!
DECLARE_SOA_COLUMN(McMask, mcMask, uint16_t);
} // namespace reducedbarreltracklabel

// NOTE: MC labels. This table has one entry for each reconstructed track (joinable with the track tables)
//          The McParticleId points to the position of the MC truth track from the ReducedTracksMC table
DECLARE_SOA_TABLE(ReducedTracksBarrelLabels, "AOD", "RTBARRELLABELS", //!
                  reducedbarreltracklabel::ReducedMCTrackId, reducedbarreltracklabel::McMask, reducedtrackMC::McReducedFlags);

using ReducedTrackBarrelLabel = ReducedTracksBarrelLabels::iterator;

// muon quantities
namespace reducedmuon
{
DECLARE_SOA_INDEX_COLUMN(ReducedEvent, reducedevent);        //!
DECLARE_SOA_COLUMN(FilteringFlags, filteringFlags, uint8_t); //!
// the (pt,eta,phi,sign) will be computed in the skimming task //!
DECLARE_SOA_COLUMN(Pt, pt, float);   //!
DECLARE_SOA_COLUMN(Eta, eta, float); //!
DECLARE_SOA_COLUMN(Phi, phi, float); //!
DECLARE_SOA_COLUMN(Sign, sign, int); //!
DECLARE_SOA_COLUMN(FwdDcaX, fwdDcaX, float);       //!  Impact parameter in X of forward track to the primary vertex
DECLARE_SOA_COLUMN(FwdDcaY, fwdDcaY, float);       //!  Impact parameter in Y of forward track to the primary vertex
DECLARE_SOA_COLUMN(IsAmbiguous, isAmbiguous, int); //!
DECLARE_SOA_DYNAMIC_COLUMN(Px, px,   //!
                           [](float pt, float phi) -> float { return pt * std::cos(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py, //!
                           [](float pt, float phi) -> float { return pt * std::sin(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz, //!
                           [](float pt, float eta) -> float { return pt * std::sinh(eta); });
DECLARE_SOA_DYNAMIC_COLUMN(P, p, //!
                           [](float pt, float eta) -> float { return pt * std::cosh(eta); });
DECLARE_SOA_COLUMN(RawPhi, rawPhi, float);           //!
DECLARE_SOA_DYNAMIC_COLUMN(MIDBoardCh1, midBoardCh1, //!
                           [](uint32_t midBoards) -> int { return static_cast<int>(midBoards & 0xFF); });
DECLARE_SOA_DYNAMIC_COLUMN(MIDBoardCh2, midBoardCh2, //!
                           [](uint32_t midBoards) -> int { return static_cast<int>((midBoards >> 8) & 0xFF); });
DECLARE_SOA_DYNAMIC_COLUMN(MIDBoardCh3, midBoardCh3, //!
                           [](uint32_t midBoards) -> int { return static_cast<int>((midBoards >> 16) & 0xFF); });
DECLARE_SOA_DYNAMIC_COLUMN(MIDBoardCh4, midBoardCh4, //!
                           [](uint32_t midBoards) -> int { return static_cast<int>((midBoards >> 24) & 0xFF); });
DECLARE_SOA_SELF_INDEX_COLUMN_FULL(MCHTrack, matchMCHTrack, int, "Muons_MatchMCHTrack");
} // namespace reducedmuon

// Muon track kinematics
DECLARE_SOA_TABLE(ReducedMuons, "AOD", "RTMUON", //!
                  o2::soa::Index<>, reducedmuon::ReducedEventId, reducedmuon::FilteringFlags,
                  reducedmuon::Pt, reducedmuon::Eta, reducedmuon::Phi, reducedmuon::Sign, reducedmuon::IsAmbiguous,
                  reducedmuon::Px<reducedmuon::Pt, reducedmuon::Phi>,
                  reducedmuon::Py<reducedmuon::Pt, reducedmuon::Phi>,
                  reducedmuon::Pz<reducedmuon::Pt, reducedmuon::Eta>,
                  reducedmuon::P<reducedmuon::Pt, reducedmuon::Eta>);

// Muon track quality details
DECLARE_SOA_TABLE(ReducedMuonsExtra, "AOD", "RTMUONEXTRA", //!
                  fwdtrack::NClusters, fwdtrack::PDca, fwdtrack::RAtAbsorberEnd,
                  fwdtrack::Chi2, fwdtrack::Chi2MatchMCHMID, fwdtrack::Chi2MatchMCHMFT,
                  fwdtrack::MatchScoreMCHMFT, reducedmuon::MCHTrackId,
                  fwdtrack::MCHBitMap, fwdtrack::MIDBitMap, fwdtrack::MIDBoards, fwdtrack::TrackType,
                  reducedmuon::FwdDcaX, reducedmuon::FwdDcaY);

// Muon covariance, TODO: the rest of the matrix should be added when needed
DECLARE_SOA_TABLE(ReducedMuonsCov, "AOD", "RTMUONCOV",
                  fwdtrack::X, fwdtrack::Y, fwdtrack::Z, reducedmuon::RawPhi, fwdtrack::Tgl, fwdtrack::Signed1Pt,
                  fwdtrack::CXX, fwdtrack::CXY, fwdtrack::CYY, fwdtrack::CPhiX, fwdtrack::CPhiY, fwdtrack::CPhiPhi,
                  fwdtrack::CTglX, fwdtrack::CTglY, fwdtrack::CTglPhi, fwdtrack::CTglTgl, fwdtrack::C1PtX,
                  fwdtrack::C1PtY, fwdtrack::C1PtPhi, fwdtrack::C1PtTgl, fwdtrack::C1Pt21Pt2);

// iterators
using ReducedMuon = ReducedMuons::iterator;
using ReducedMuonExtra = ReducedMuonsExtra::iterator;
using ReducedMuonCov = ReducedMuonsCov::iterator;

namespace reducedmuonlabel
{
DECLARE_SOA_INDEX_COLUMN(ReducedMCTrack, reducedMCTrack); //!
DECLARE_SOA_COLUMN(McMask, mcMask, uint16_t);
DECLARE_SOA_COLUMN(McReducedFlags, mcReducedFlags, uint16_t);
} // namespace reducedmuonlabel
// NOTE: MC labels. This table has one entry for each reconstructed muon (joinable with the muon tables)
//          The McParticleId points to the position of the MC truth track from the ReducedTracksMC table
DECLARE_SOA_TABLE(ReducedMuonsLabels, "AOD", "RTMUONSLABELS", //!
                  reducedmuonlabel::ReducedMCTrackId, reducedmuonlabel::McMask, reducedtrackMC::McReducedFlags);

using ReducedMuonsLabel = ReducedMuonsLabels::iterator;

namespace dilepton_track_index
{
DECLARE_SOA_INDEX_COLUMN_FULL(Index0, index0, int, ReducedMuons, "_0"); //! Index to first prong
DECLARE_SOA_INDEX_COLUMN_FULL(Index1, index1, int, ReducedMuons, "_1"); //! Index to second prong
DECLARE_SOA_COLUMN(Pt1, pt1, float);                                    //! Pt of the first prong
DECLARE_SOA_COLUMN(Eta1, eta1, float);                                  //! Eta of the first prong
DECLARE_SOA_COLUMN(Phi1, phi1, float);                                  //! Phi of the first prong
DECLARE_SOA_COLUMN(Sign1, sign1, int);                                  //! Sign of the first prong

DECLARE_SOA_COLUMN(Pt2, pt2, float);   //! Pt of the second prong
DECLARE_SOA_COLUMN(Eta2, eta2, float); //! Eta of the second prong
DECLARE_SOA_COLUMN(Phi2, phi2, float); //! Phi of the second prong
DECLARE_SOA_COLUMN(Sign2, sign2, int); //! Sign of the second prong

DECLARE_SOA_COLUMN(McMask1, mcMask1, uint16_t); //! MC mask of the MCLabel of the first prong
DECLARE_SOA_COLUMN(McMask2, mcMask2, uint16_t); //! MC mask of the MCLabel of the second prong

DECLARE_SOA_COLUMN(Chi2MatchMCHMID1, chi2MatchMCHMID1, float); //! MCH-MID Match Chi2 for MUONStandalone tracks
DECLARE_SOA_COLUMN(Chi2MatchMCHMFT1, chi2MatchMCHMFT1, float); //! MCH-MFT Match Chi2 for GlobalMuonTracks

DECLARE_SOA_COLUMN(Chi2MatchMCHMID2, chi2MatchMCHMID2, float); //! MCH-MID Match Chi2 for MUONStandalone tracks
DECLARE_SOA_COLUMN(Chi2MatchMCHMFT2, chi2MatchMCHMFT2, float); //! MCH-MFT Match Chi2 for GlobalMuonTracks

DECLARE_SOA_COLUMN(PtMC1, ptMC1, float);   //! MC Pt of the first prong
DECLARE_SOA_COLUMN(EtaMC1, etaMC1, float); //! MC Eta of the first prong
DECLARE_SOA_COLUMN(PhiMC1, phiMC1, float); //! MC Phi of the first prong
DECLARE_SOA_COLUMN(EMC1, eMC1, float);     //! MC Energy of the first prong

DECLARE_SOA_COLUMN(PtMC2, ptMC2, float);   //! MC Pt of the second prong
DECLARE_SOA_COLUMN(EtaMC2, etaMC2, float); //! MC Eta of the second prong
DECLARE_SOA_COLUMN(PhiMC2, phiMC2, float); //! MC Phi of the second prong
DECLARE_SOA_COLUMN(EMC2, eMC2, float);     //! MC Energy of the second prong

DECLARE_SOA_COLUMN(Vx1, vx1, float); //! X production vertex in cm
DECLARE_SOA_COLUMN(Vy1, vy1, float); //! Y production vertex in cm
DECLARE_SOA_COLUMN(Vz1, vz1, float); //! Z production vertex in cm
DECLARE_SOA_COLUMN(Vt1, vt1, float); //! Production vertex time

DECLARE_SOA_COLUMN(Vx2, vx2, float); //! X production vertex in cm
DECLARE_SOA_COLUMN(Vy2, vy2, float); //! Y production vertex in cm
DECLARE_SOA_COLUMN(Vz2, vz2, float); //! Z production vertex in cm
DECLARE_SOA_COLUMN(Vt2, vt2, float); //! Production vertex time

} // namespace dilepton_track_index

// pair information
namespace reducedpair
{
DECLARE_SOA_INDEX_COLUMN(ReducedEvent, reducedevent); //!
DECLARE_SOA_COLUMN(Mass, mass, float);                //!
DECLARE_SOA_COLUMN(Pt, pt, float);                    //!
DECLARE_SOA_COLUMN(Eta, eta, float);                  //!
DECLARE_SOA_COLUMN(Phi, phi, float);                  //!
DECLARE_SOA_COLUMN(Sign, sign, int);                  //!
DECLARE_SOA_COLUMN(FilterMap, filterMap, uint32_t);   //!
DECLARE_SOA_COLUMN(McDecision, mcDecision, uint32_t); //!
DECLARE_SOA_COLUMN(Tauz, tauz, float);                //! Longitudinal pseudo-proper time of lepton pair (in ns)
DECLARE_SOA_COLUMN(TauzErr, tauzErr, float);          //! Error on longitudinal pseudo-proper time of lepton pair (in ns)
DECLARE_SOA_COLUMN(Tauxy, tauxy, float);              //! Transverse pseudo-proper time of lepton pair (in ns)
DECLARE_SOA_COLUMN(TauxyErr, tauxyErr, float);        //! Error on transverse pseudo-proper time of lepton pair (in ns)
DECLARE_SOA_COLUMN(Lz, lz, float);                    //! Longitudinal projection of decay length
DECLARE_SOA_COLUMN(Lxy, lxy, float);                  //! Transverse projection of decay length
DECLARE_SOA_COLUMN(U2Q2, u2q2, float);                //! Scalar product between unitary vector with event flow vector (harmonic 2)
DECLARE_SOA_COLUMN(U3Q3, u3q3, float);                //! Scalar product between unitary vector with event flow vector (harmonic 3)
DECLARE_SOA_COLUMN(Cos2DeltaPhi, cos2deltaphi, float); //! Cosinus term using event plane angle (harmonic 2)
DECLARE_SOA_COLUMN(Cos3DeltaPhi, cos3deltaphi, float); //! Cosinus term using event plane angle (harmonic 3)
// DECLARE_SOA_INDEX_COLUMN(ReducedMuon, reducedmuon2); //!
DECLARE_SOA_DYNAMIC_COLUMN(Px, px, //!
                           [](float pt, float phi) -> float { return pt * std::cos(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py, //!
                           [](float pt, float phi) -> float { return pt * std::sin(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz, //!
                           [](float pt, float eta) -> float { return pt * std::sinh(eta); });
DECLARE_SOA_DYNAMIC_COLUMN(P, p, //!
                           [](float pt, float eta) -> float { return pt * std::cosh(eta); });
DECLARE_SOA_DYNAMIC_COLUMN(Rap, rap, //!
                           [](float pt, float eta, float m) -> float { return std::log((std::sqrt(m * m + pt * pt * std::cosh(eta) * std::cosh(eta)) + pt * std::sinh(eta)) / std::sqrt(m * m + pt * pt)); });
} // namespace reducedpair

DECLARE_SOA_TABLE(Dileptons, "AOD", "RTDILEPTON", //!
                  reducedpair::ReducedEventId,
                  reducedpair::Mass,
                  reducedpair::Pt, reducedpair::Eta, reducedpair::Phi, reducedpair::Sign,
                  reducedpair::FilterMap,
                  reducedpair::McDecision,
                  reducedpair::Px<reducedpair::Pt, reducedpair::Phi>,
                  reducedpair::Py<reducedpair::Pt, reducedpair::Phi>,
                  reducedpair::Pz<reducedpair::Pt, reducedpair::Eta>,
                  reducedpair::P<reducedpair::Pt, reducedpair::Eta>,
                  reducedpair::Rap<reducedpair::Pt, reducedpair::Eta, reducedpair::Mass>);

DECLARE_SOA_TABLE(DileptonsExtra, "AOD", "RTDILEPTONEXTRA", //!
                  dilepton_track_index::Index0Id, dilepton_track_index::Index1Id,
                  reducedpair::Tauz,
                  reducedpair::Lz,
                  reducedpair::Lxy);

DECLARE_SOA_TABLE(DileptonsFlow, "AOD", "RTDILEPTONFLOW", //!
                  reducedpair::U2Q2,
                  reducedpair::U3Q3,
                  reducedpair::Cos2DeltaPhi,
                  reducedpair::Cos3DeltaPhi);

DECLARE_SOA_TABLE(DimuonsAll, "AOD", "RTDIMUONALL", //!
                  collision::PosX, collision::PosY, collision::PosZ,
                  reducedevent::MCPosX, reducedevent::MCPosY, reducedevent::MCPosZ,
                  reducedpair::Mass,
                  reducedpair::McDecision,
                  reducedpair::Pt, reducedpair::Eta, reducedpair::Phi, reducedpair::Sign,
                  reducedpair::Tauz, reducedpair::TauzErr,
                  reducedpair::Tauxy, reducedpair::TauxyErr,
                  dilepton_track_index::Pt1, dilepton_track_index::Eta1, dilepton_track_index::Phi1, dilepton_track_index::Sign1,
                  dilepton_track_index::Pt2, dilepton_track_index::Eta2, dilepton_track_index::Phi2, dilepton_track_index::Sign2,
                  dilepton_track_index::McMask1, dilepton_track_index::McMask2,
                  dilepton_track_index::Chi2MatchMCHMID1, dilepton_track_index::Chi2MatchMCHMID2,
                  dilepton_track_index::Chi2MatchMCHMFT1, dilepton_track_index::Chi2MatchMCHMFT2,
                  dilepton_track_index::PtMC1, dilepton_track_index::EtaMC1, dilepton_track_index::PhiMC1, dilepton_track_index::EMC1,
                  dilepton_track_index::PtMC2, dilepton_track_index::EtaMC2, dilepton_track_index::PhiMC2, dilepton_track_index::EMC2,
                  dilepton_track_index::Vx1, dilepton_track_index::Vy1, dilepton_track_index::Vz1, dilepton_track_index::Vt1,
                  dilepton_track_index::Vx2, dilepton_track_index::Vy2, dilepton_track_index::Vz2, dilepton_track_index::Vt2);

using Dilepton = Dileptons::iterator;
using DileptonExtra = DileptonsExtra::iterator;
using DileptonFlow = DileptonsFlow::iterator;
using DimuonAll = DimuonsAll::iterator;

// candidate information
namespace dileptonTrackCandidate
{
DECLARE_SOA_INDEX_COLUMN(ReducedEvent, reducedevent); //!
DECLARE_SOA_COLUMN(McDecision, mcDecision, uint32_t); //!
DECLARE_SOA_COLUMN(Mass, mass, float);                //!
DECLARE_SOA_COLUMN(Pt, pt, float);                    //!
DECLARE_SOA_COLUMN(Eta, eta, float);                  //!
DECLARE_SOA_COLUMN(Tauz, tauz, float);                //!
DECLARE_SOA_COLUMN(Tauxy, tauxy, float);              //!
DECLARE_SOA_COLUMN(Lz, lz, float);                    //! Longitudinal projection of decay length
DECLARE_SOA_COLUMN(Lxy, lxy, float);                  //! Transverse projection of decay length
} // namespace dileptonTrackCandidate

DECLARE_SOA_TABLE(DileptonTrackCandidates, "AOD", "RTDILEPTONTRACK", //!
                  dileptonTrackCandidate::McDecision,
                  dileptonTrackCandidate::Mass,
                  dileptonTrackCandidate::Pt,
                  dileptonTrackCandidate::Eta,
                  dileptonTrackCandidate::Tauz,
                  dileptonTrackCandidate::Tauxy,
                  dileptonTrackCandidate::Lz,
                  dileptonTrackCandidate::Lxy);

using DileptonTrackCandidate = DileptonTrackCandidates::iterator;

namespace v0bits
{
DECLARE_SOA_COLUMN(PIDBit, pidbit, uint8_t); //!
} // namespace v0bits

// bit information for particle species.
DECLARE_SOA_TABLE(V0Bits, "AOD", "V0BITS", //!
                  v0bits::PIDBit);

// iterators
using V0Bit = V0Bits::iterator;

namespace DalBits
{
DECLARE_SOA_COLUMN(DALITZBits, dalitzBits, int); //!
} // namespace DalBits

// bit information for particle species.
DECLARE_SOA_TABLE(DalitzBits, "AOD", "DALITZBITS", DalBits::DALITZBits);
} // namespace o2::aod

#endif // O2_Analysis_ReducedInfoTables_H_
