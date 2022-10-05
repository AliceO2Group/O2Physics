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
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include <cmath>

namespace o2::aod
{

namespace udmccollision
{
DECLARE_SOA_COLUMN(GlobalBC, globalBC, uint64_t);
} // namespace udmccollision

DECLARE_SOA_TABLE(UDMcCollisions, "AOD", "UDMCCOLLISIONS",
                  o2::soa::Index<>,
                  udmccollision::GlobalBC,
                  mccollision::GeneratorsID,
                  mccollision::PosX,
                  mccollision::PosY,
                  mccollision::PosZ,
                  mccollision::T,
                  mccollision::Weight,
                  mccollision::ImpactParameter);

using UDMcCollision = UDMcCollisions::iterator;

namespace udmcparticle
{
DECLARE_SOA_INDEX_COLUMN(UDMcCollision, udMcCollision);    //!
DECLARE_SOA_SELF_ARRAY_INDEX_COLUMN(Mothers, mothers);     //! Mother tracks (possible empty) array. Iterate over mcParticle.mothers_as<aod::McParticles>())
DECLARE_SOA_SELF_SLICE_INDEX_COLUMN(Daughters, daughters); //! Daughter tracks (possibly empty) slice. Check for non-zero with mcParticle.has_daughters(). Iterate over mcParticle.daughters_as<aod::McParticles>())
DECLARE_SOA_COLUMN(Px, px, float);                         //!
DECLARE_SOA_COLUMN(Py, py, float);                         //!
DECLARE_SOA_COLUMN(Pz, pz, float);                         //!
DECLARE_SOA_COLUMN(E, e, float);                           //!
} // namespace udmcparticle

DECLARE_SOA_TABLE_FULL(UDMcParticles, "UDMcParticles", "AOD", "UDMCPARTICLES", //!
                       o2::soa::Index<>, udmcparticle::UDMcCollisionId,
                       mcparticle::PdgCode,
                       mcparticle::StatusCode,
                       mcparticle::Flags,
                       udmcparticle::MothersIds,
                       udmcparticle::DaughtersIdSlice,
                       mcparticle::Weight,
                       udmcparticle::Px,
                       udmcparticle::Py,
                       udmcparticle::Pz,
                       udmcparticle::E,
                       mcparticle::ProducedByGenerator<mcparticle::Flags>,
                       mcparticle::FromBackgroundEvent<mcparticle::Flags>,
                       mcparticle::GetGenStatusCode<mcparticle::Flags, mcparticle::StatusCode>,
                       mcparticle::GetProcess<mcparticle::Flags, mcparticle::StatusCode>,
                       mcparticle::IsPhysicalPrimary<mcparticle::Flags>);

using UDMcParticle = UDMcParticles::iterator;

namespace udcollision
{
// general information
DECLARE_SOA_COLUMN(RunNumber, runNumber, int32_t); //! run number
DECLARE_SOA_COLUMN(GlobalBC, globalBC, uint64_t);  //! global BC instead of BC ID since candidate may not have a corresponding record in BCs table
DECLARE_SOA_COLUMN(NetCharge, netCharge, int8_t);  //! Sum of track signs
DECLARE_SOA_COLUMN(RgtrwTOF, rgtrwTOF, float);     //! Fraction of global tracks with TOF hit
// FT0 information
DECLARE_SOA_COLUMN(TotalFT0AmplitudeA, totalFT0AmplitudeA, float); //! sum of amplitudes on A side of FT0
DECLARE_SOA_COLUMN(TotalFT0AmplitudeC, totalFT0AmplitudeC, float); //! sum of amplitudes on C side of FT0
DECLARE_SOA_COLUMN(TimeFT0A, timeFT0A, float);                     //! FT0A average time
DECLARE_SOA_COLUMN(TimeFT0C, timeFT0C, float);                     //! FT0C average time
DECLARE_SOA_COLUMN(TriggerMaskFT0, triggerMaskFT0, uint8_t);       //! FT0 trigger mask
// FDD information
DECLARE_SOA_COLUMN(TotalFDDAmplitudeA, totalFDDAmplitudeA, float); //! sum of amplitudes on A side of FDD
DECLARE_SOA_COLUMN(TotalFDDAmplitudeC, totalFDDAmplitudeC, float); //! sum of amplitudes on C side of FDD
DECLARE_SOA_COLUMN(TimeFDDA, timeFDDA, float);                     //! FDDA average time
DECLARE_SOA_COLUMN(TimeFDDC, timeFDDC, float);                     //! FDDC average time
DECLARE_SOA_COLUMN(TriggerMaskFDD, triggerMaskFDD, uint8_t);       //! FDD trigger mask
// FV0A information
DECLARE_SOA_COLUMN(TotalFV0AmplitudeA, totalFV0AmplitudeA, float); //! sum of amplitudes on A side of FDD
DECLARE_SOA_COLUMN(TimeFV0A, timeFV0A, float);                     //! FV0A average time
DECLARE_SOA_COLUMN(TriggerMaskFV0A, triggerMaskFV0A, uint8_t);     //! FV0 trigger mask
// FIT selection flags
DECLARE_SOA_COLUMN(BBFT0A, bbFT0A, bool); //! Beam-beam time in FT0A
DECLARE_SOA_COLUMN(BBFT0C, bbFT0C, bool); //! Beam-beam time in FT0C
DECLARE_SOA_COLUMN(BGFT0A, bgFT0A, bool); //! Beam-gas time in FT0A
DECLARE_SOA_COLUMN(BGFT0C, bgFT0C, bool); //! Beam-gas time in FT0C
DECLARE_SOA_COLUMN(BBFV0A, bbFV0A, bool); //! Beam-beam time in V0A
DECLARE_SOA_COLUMN(BGFV0A, bgFV0A, bool); //! Beam-gas time in V0A
DECLARE_SOA_COLUMN(BBFDDA, bbFDDA, bool); //! Beam-beam time in FDA
DECLARE_SOA_COLUMN(BBFDDC, bbFDDC, bool); //! Beam-beam time in FDC
DECLARE_SOA_COLUMN(BGFDDA, bgFDDA, bool); //! Beam-gas time in FDA
DECLARE_SOA_COLUMN(BGFDDC, bgFDDC, bool); //! Beam-gas time in FDC

} // namespace udcollision

DECLARE_SOA_TABLE(UDCollisions, "AOD", "UDCOLLISION",
                  o2::soa::Index<>,
                  udcollision::GlobalBC,
                  udcollision::RunNumber,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  collision::NumContrib,
                  udcollision::NetCharge,
                  udcollision::RgtrwTOF);

DECLARE_SOA_TABLE(UDCollisionsSels, "AOD", "UDCOLLISIONSEL",
                  udcollision::TotalFT0AmplitudeA,
                  udcollision::TotalFT0AmplitudeC,
                  udcollision::TimeFT0A,
                  udcollision::TimeFT0C,
                  udcollision::TriggerMaskFT0,
                  udcollision::TotalFDDAmplitudeA,
                  udcollision::TotalFDDAmplitudeC,
                  udcollision::TimeFDDA,
                  udcollision::TimeFDDC,
                  udcollision::TriggerMaskFDD,
                  udcollision::TotalFV0AmplitudeA,
                  udcollision::TimeFV0A,
                  udcollision::TriggerMaskFV0A,
                  udcollision::BBFT0A, udcollision::BBFT0C, udcollision::BGFT0A, udcollision::BGFT0C,
                  udcollision::BBFV0A, udcollision::BGFV0A,
                  udcollision::BBFDDA, udcollision::BBFDDC, udcollision::BGFDDA, udcollision::BGFDDC);

using UDCollision = UDCollisions::iterator;
using UDCollisionsSel = UDCollisionsSels::iterator;

namespace udtrack
{
DECLARE_SOA_INDEX_COLUMN(UDCollision, udCollision);    //!
DECLARE_SOA_COLUMN(Px, px, float);                     //!
DECLARE_SOA_COLUMN(Py, py, float);                     //!
DECLARE_SOA_COLUMN(Pz, pz, float);                     //!
DECLARE_SOA_COLUMN(Sign, sign, int);                   //!
DECLARE_SOA_COLUMN(GlobalBC, globalBC, uint64_t);      //!
DECLARE_SOA_COLUMN(TrackTime, trackTime, double);      //!
DECLARE_SOA_COLUMN(TrackTimeRes, trackTimeRes, float); //! time resolution
DECLARE_SOA_COLUMN(DetectorMap, detectorMap, uint8_t); //!
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt,                     //!
                           [](float px, float py) -> float {
                             return std::sqrt(px * px + py * py);
                           });
} // namespace udtrack

// Barrel track kinematics
DECLARE_SOA_TABLE(UDTracks, "AOD", "UDTRACK",
                  o2::soa::Index<>,
                  udtrack::Px,
                  udtrack::Py,
                  udtrack::Pz,
                  udtrack::Sign,
                  udtrack::GlobalBC,
                  udtrack::TrackTime,
                  udtrack::TrackTimeRes,
                  udtrack::Pt<udtrack::Px, udtrack::Py>);

DECLARE_SOA_TABLE(UDTracksCov, "AOD", "UDTRACKCOV",
                  track::X, track::Y, track::Z,
                  track::SigmaY, track::SigmaZ);

DECLARE_SOA_TABLE(UDTracksPID, "AOD", "UDTRACKPID",
                  pidtpc::TPCNSigmaEl, pidtpc::TPCNSigmaMu, pidtpc::TPCNSigmaPi, pidtpc::TPCNSigmaKa, pidtpc::TPCNSigmaPr,
                  pidtof::TOFNSigmaEl, pidtof::TOFNSigmaMu, pidtof::TOFNSigmaPi, pidtof::TOFNSigmaKa, pidtof::TOFNSigmaPr);

DECLARE_SOA_TABLE(UDTrackCollisionIDs, "AOD", "UDTRCOLID",
                  udtrack::UDCollisionId);

DECLARE_SOA_TABLE(UDTracksExtra, "AOD", "UDTRACKEXTRA",
                  track::ITSClusterMap,
                  track::TPCNClsFindable,
                  track::TPCNClsFindableMinusFound,
                  track::TPCNClsFindableMinusCrossedRows,
                  track::TPCNClsShared,
                  track::TRDPattern,
                  track::ITSChi2NCl,
                  track::TPCChi2NCl,
                  track::TRDChi2,
                  track::TOFChi2,
                  track::TPCSignal,
                  pidtofsignal::TOFSignal,
                  track::TRDSignal,
                  track::Length,
                  track::TOFExpMom,
                  udtrack::DetectorMap,
                  track::HasITS<udtrack::DetectorMap>,
                  track::HasTPC<udtrack::DetectorMap>,
                  track::HasTRD<udtrack::DetectorMap>,
                  track::HasTOF<udtrack::DetectorMap>,
                  track::ITSNCls<track::ITSClusterMap>,
                  track::TPCNClsCrossedRows<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>);

DECLARE_SOA_TABLE(UDTracksDCA, "AOD", "UDTRACKDCA",
                  track::DcaZ,
                  track::DcaXY)

using UDTrack = UDTracks::iterator;
using UDTrackCov = UDTracksCov::iterator;
using UDTrackExtra = UDTracksExtra::iterator;
using UDTrackDCA = UDTracksDCA::iterator;
using UDTrackCollisionID = UDTrackCollisionIDs::iterator;

namespace udmctracklabel
{
DECLARE_SOA_INDEX_COLUMN(UDMcParticle, udMcParticle);
DECLARE_SOA_COLUMN(McMask, mcMask, uint16_t);
} // namespace udmctracklabel

DECLARE_SOA_TABLE(UDMcTrackLabels, "AOD", "UDMCTRACKLABEL",
                  udmctracklabel::UDMcParticleId,
                  udmctracklabel::McMask);

using UDMcTrackLabel = UDMcTrackLabels::iterator;

// only MCH-MID tracks
namespace udfwdtrack
{
DECLARE_SOA_INDEX_COLUMN(UDCollision, udCollision);    //!
DECLARE_SOA_COLUMN(Px, px, float);                     //!
DECLARE_SOA_COLUMN(Py, py, float);                     //!
DECLARE_SOA_COLUMN(Pz, pz, float);                     //!
DECLARE_SOA_COLUMN(Sign, sign, int);                   //!
DECLARE_SOA_COLUMN(GlobalBC, globalBC, uint64_t);      //!
DECLARE_SOA_COLUMN(TrackTime, trackTime, double);      //!
DECLARE_SOA_COLUMN(TrackTimeRes, trackTimeRes, float); //! time resolution

DECLARE_SOA_DYNAMIC_COLUMN(MIDBoardCh1, midBoardCh1, //!
                           [](uint32_t midBoards) -> int { return static_cast<int>(midBoards & 0xFF); });
DECLARE_SOA_DYNAMIC_COLUMN(MIDBoardCh2, midBoardCh2, //!
                           [](uint32_t midBoards) -> int { return static_cast<int>((midBoards >> 8) & 0xFF); });
DECLARE_SOA_DYNAMIC_COLUMN(MIDBoardCh3, midBoardCh3, //!
                           [](uint32_t midBoards) -> int { return static_cast<int>((midBoards >> 16) & 0xFF); });
DECLARE_SOA_DYNAMIC_COLUMN(MIDBoardCh4, midBoardCh4, //!
                           [](uint32_t midBoards) -> int { return static_cast<int>((midBoards >> 24) & 0xFF); });
} // namespace udfwdtrack

// Muon track kinematics
DECLARE_SOA_TABLE(UDFwdTracks, "AOD", "UDFWDTRACK",
                  o2::soa::Index<>,
                  udfwdtrack::Px,
                  udfwdtrack::Py,
                  udfwdtrack::Pz,
                  udfwdtrack::Sign,
                  udfwdtrack::GlobalBC,
                  udfwdtrack::TrackTime,
                  udfwdtrack::TrackTimeRes);

DECLARE_SOA_TABLE(UDFwdTrackCollisionIDs, "AOD", "UDFWDTRCOLID",
                  udfwdtrack::UDCollisionId);

// Muon track quality details
DECLARE_SOA_TABLE(UDFwdTracksExtra, "AOD", "UDFWDTRACKEXTRA",
                  fwdtrack::NClusters,
                  fwdtrack::PDca,
                  fwdtrack::RAtAbsorberEnd,
                  fwdtrack::Chi2,
                  fwdtrack::Chi2MatchMCHMID,
                  fwdtrack::MCHBitMap,
                  fwdtrack::MIDBitMap,
                  fwdtrack::MIDBoards);

using UDFwdTrack = UDFwdTracks::iterator;
using UDFwdTrackExtra = UDFwdTracksExtra::iterator;
using UDFwdTrackCollisionID = UDFwdTrackCollisionIDs::iterator;

namespace udmcfwdtracklabel
{
DECLARE_SOA_INDEX_COLUMN(UDMcParticle, udMcParticle);
DECLARE_SOA_COLUMN(McMask, mcMask, uint16_t);
} // namespace udmcfwdtracklabel

DECLARE_SOA_TABLE(UDMcFwdTrackLabels, "AOD", "UDMCFWDTRLABEL",
                  udmcfwdtracklabel::UDMcParticleId,
                  udmcfwdtracklabel::McMask);

using UDMcFwdTrackLabel = UDMcFwdTrackLabels::iterator;

} // namespace o2::aod

#endif // O2PHYSICS_UDTABLES_H
