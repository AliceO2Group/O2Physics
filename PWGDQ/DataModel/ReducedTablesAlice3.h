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
/// \file ReducedTablesAlice3.h
///
/// \brief Reduced DQ table definitions for ALICE 3
///
/// \author Alexander Tiekoetter (atiekoet@cern.ch) University of Muenster

#ifndef PWGDQ_DATAMODEL_REDUCEDTABLESALICE3_H_
#define PWGDQ_DATAMODEL_REDUCEDTABLESALICE3_H_

#include "PWGDQ/DataModel/ReducedInfoTables.h"

#include "ALICE3/DataModel/OTFPIDTrk.h"
#include "ALICE3/DataModel/OTFRICH.h"
#include "ALICE3/DataModel/OTFTOF.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace reducedeventalice3
{
DECLARE_SOA_COLUMN(MultDensity, multDensity, float);
DECLARE_SOA_COLUMN(MCPosX, mcPosX, float); //!  MC event position X
DECLARE_SOA_COLUMN(MCPosY, mcPosY, float); //!  MC event position Y
DECLARE_SOA_COLUMN(MCPosZ, mcPosZ, float); //!  MC event position Z
} // namespace reducedeventalice3

namespace reducedeventmcalice3
{
DECLARE_SOA_COLUMN(MultMCNParticlesEta05, multMCNParticlesEta05, float);
DECLARE_SOA_COLUMN(MultMCNParticlesEta08, multMCNParticlesEta08, float);
DECLARE_SOA_COLUMN(MultMCNParticlesEta10, multMCNParticlesEta10, float);
DECLARE_SOA_COLUMN(MultMCNParticlesEta20, multMCNParticlesEta20, float);
DECLARE_SOA_COLUMN(MultMCNParticlesEta40, multMCNParticlesEta40, float);
} // namespace reducedeventmcalice3

DECLARE_SOA_TABLE_STAGED(ReducedA3Events, "REA3EVENTS", //!   Main event information table
                         o2::soa::Index<>,
                         collision::PosX, collision::PosY, collision::PosZ, collision::NumContrib,
                         collision::CollisionTime, collision::CollisionTimeRes, reducedeventalice3::MultDensity);

DECLARE_SOA_TABLE(ReducedA3EventsVtxCov, "AOD", "REA3VTXCOV", //!    Event vertex covariance matrix
                  collision::CovXX, collision::CovXY, collision::CovXZ,
                  collision::CovYY, collision::CovYZ, collision::CovZZ, collision::Chi2);

DECLARE_SOA_TABLE(ReducedA3EventsInfo, "AOD", "REA3EVENTINFO", //!   Main event index table
                  reducedevent::CollisionId);

DECLARE_SOA_TABLE(ReducedA3MCEvents, "AOD", "REA3MCEVTS", //!   Event level MC truth information
                  o2::soa::Index<>,
                  mccollision::GeneratorsID, reducedeventalice3::MCPosX, reducedeventalice3::MCPosY, reducedeventalice3::MCPosZ,
                  mccollision::T, mccollision::Weight, mccollision::ImpactParameter);

using ReducedA3MCEvent = ReducedA3MCEvents::iterator;
using ReducedA3Event = ReducedA3Events::iterator;

namespace reducedtrackalice3
{
// basic track information
DECLARE_SOA_INDEX_COLUMN(ReducedA3Event, reduceda3event); //!
DECLARE_SOA_INDEX_COLUMN(Track, track);                   //!
// ----  flags reserved for storing various information during filtering
DECLARE_SOA_BITMAP_COLUMN(FilteringFlags, filteringFlags, 64); //!
// -----------------------------------------------------

DECLARE_SOA_COLUMN(IsReconstructed, isReconstructed, bool);
DECLARE_SOA_COLUMN(NSiliconHits, nSiliconHits, int);
DECLARE_SOA_COLUMN(NTPCHits, nTPCHits, int);

DECLARE_SOA_COLUMN(Pt, pt, float);                     //!
DECLARE_SOA_COLUMN(Eta, eta, float);                   //!
DECLARE_SOA_COLUMN(Phi, phi, float);                   //!
DECLARE_SOA_COLUMN(Sign, sign, int);                   //!
DECLARE_SOA_COLUMN(IsAmbiguous, isAmbiguous, int);     //!
DECLARE_SOA_COLUMN(DcaXY, dcaXY, float);               //!
DECLARE_SOA_COLUMN(DcaZ, dcaZ, float);                 //!
DECLARE_SOA_COLUMN(DetectorMap, detectorMap, uint8_t); //! Detector map: see enum DetectorMapEnum
DECLARE_SOA_INDEX_COLUMN(Collision, collision);        //!
DECLARE_SOA_DYNAMIC_COLUMN(HasITS, hasITS,             //! Flag to check if track has a ITS match
                           [](uint8_t detectorMap) -> bool { return detectorMap & o2::aod::track::ITS; });
DECLARE_SOA_DYNAMIC_COLUMN(HasTPC, hasTPC, //! Flag to check if track has a TPC match
                           [](uint8_t detectorMap) -> bool { return detectorMap & o2::aod::track::TPC; });
DECLARE_SOA_DYNAMIC_COLUMN(HasTRD, hasTRD, //! Flag to check if track has a TRD match
                           [](uint8_t detectorMap) -> bool { return detectorMap & o2::aod::track::TRD; });
DECLARE_SOA_DYNAMIC_COLUMN(HasTOF, hasTOF, //! Flag to check if track has a TOF measurement
                           [](uint8_t detectorMap) -> bool { return detectorMap & o2::aod::track::TOF; });
DECLARE_SOA_DYNAMIC_COLUMN(Px, px, //!
                           [](float pt, float phi) -> float { return pt * std::cos(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py, //!
                           [](float pt, float phi) -> float { return pt * std::sin(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz, //!
                           [](float pt, float eta) -> float { return pt * std::sinh(eta); });
DECLARE_SOA_DYNAMIC_COLUMN(P, p, //!
                           [](float pt, float eta) -> float { return pt * std::cosh(eta); });
} // namespace reducedtrackalice3

// basic track information
DECLARE_SOA_TABLE(ReducedA3Tracks, "AOD", "REA3TRACK", //!
                  o2::soa::Index<>, reducedtrackalice3::ReducedA3EventId, reducedtrackalice3::FilteringFlags,
                  reducedtrackalice3::Pt, reducedtrackalice3::Eta, reducedtrackalice3::Phi, reducedtrackalice3::Sign, reducedtrackalice3::IsAmbiguous,
                  reducedtrackalice3::Px<reducedtrackalice3::Pt, reducedtrackalice3::Phi>,
                  reducedtrackalice3::Py<reducedtrackalice3::Pt, reducedtrackalice3::Phi>,
                  reducedtrackalice3::Pz<reducedtrackalice3::Pt, reducedtrackalice3::Eta>,
                  reducedtrackalice3::P<reducedtrackalice3::Pt, reducedtrackalice3::Eta>);

DECLARE_SOA_TABLE(ReducedA3TracksBarrelCov, "AOD", "REA3BARRELCOV", //!
                  track::CYY, track::CZY, track::CZZ, track::CSnpY, track::CSnpZ,
                  track::CSnpSnp, track::CTglY, track::CTglZ, track::CTglSnp, track::CTglTgl,
                  track::C1PtY, track::C1PtZ, track::C1PtSnp, track::C1PtTgl, track::C1Pt21Pt2);

namespace reducedA3trackMC
{
DECLARE_SOA_INDEX_COLUMN(ReducedA3MCEvent, reducedA3MCEvent);                                 //!
DECLARE_SOA_COLUMN(McReducedFlags, mcReducedFlags, uint16_t);                                 //! Flags to hold compressed MC selection information
DECLARE_SOA_SELF_INDEX_COLUMN_FULL(Mother0, mother0, int, "ReducedA3MCTracks_Mother0");       //! Track index of the first mother
DECLARE_SOA_SELF_INDEX_COLUMN_FULL(Mother1, mother1, int, "ReducedA3MCTracks_Mother1");       //! Track index of the last mother
DECLARE_SOA_SELF_INDEX_COLUMN_FULL(Daughter0, daughter0, int, "ReducedA3MCTracks_Daughter0"); //! Track index of the first daughter
DECLARE_SOA_SELF_INDEX_COLUMN_FULL(Daughter1, daughter1, int, "ReducedA3MCTracks_Daughter1"); //! Track index of the last daughter
DECLARE_SOA_SELF_ARRAY_INDEX_COLUMN(Mothers, mothers);                                        //! Mother tracks (possible empty) array. Iterate over mcParticle.mothers_as<aod::McParticles>())
DECLARE_SOA_SELF_SLICE_INDEX_COLUMN(Daughters, daughters);                                    //! Daughter tracks (possibly empty) slice. Check for non-zero with mcParticle.has_daughters(). Iterate over mcParticle.daughters_as<aod::McParticles>())
DECLARE_SOA_COLUMN(Pt, pt, float);                                                            //!
DECLARE_SOA_COLUMN(Eta, eta, float);                                                          //!
DECLARE_SOA_COLUMN(Phi, phi, float);                                                          //!
DECLARE_SOA_COLUMN(E, e, float);                                                              //!
DECLARE_SOA_DYNAMIC_COLUMN(Px, px,                                                            //!
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
} // namespace reducedA3trackMC

// NOTE: This table is nearly identical to the one from Framework (except that it points to the event ID, not the BC id)
//       This table contains all MC truth tracks (both barrel and muon)
DECLARE_SOA_TABLE(ReducedA3MCTracks, "AOD", "REA3MCTRACK", //!  MC track information (on disk)
                  o2::soa::Index<>, reducedA3trackMC::ReducedA3MCEventId,
                  mcparticle::PdgCode, mcparticle::StatusCode, mcparticle::Flags,
                  reducedA3trackMC::MothersIds, reducedA3trackMC::DaughtersIdSlice,
                  mcparticle::Weight,
                  reducedA3trackMC::Pt, reducedA3trackMC::Eta, reducedA3trackMC::Phi, reducedA3trackMC::E,
                  mcparticle::Vx, mcparticle::Vy, mcparticle::Vz, mcparticle::Vt,
                  reducedA3trackMC::McReducedFlags,
                  reducedA3trackMC::Px<reducedA3trackMC::Pt, reducedA3trackMC::Phi>,
                  reducedA3trackMC::Py<reducedA3trackMC::Pt, reducedA3trackMC::Phi>,
                  reducedA3trackMC::Pz<reducedA3trackMC::Pt, reducedA3trackMC::Eta>,
                  reducedA3trackMC::P<reducedA3trackMC::Pt, reducedA3trackMC::Eta>,
                  reducedA3trackMC::Y<reducedA3trackMC::Pt, reducedA3trackMC::Eta, reducedA3trackMC::E>,
                  mcparticle::ProducedByGenerator<mcparticle::Flags>,
                  mcparticle::FromBackgroundEvent<mcparticle::Flags>,
                  mcparticle::GetGenStatusCode<mcparticle::Flags, mcparticle::StatusCode>,
                  mcparticle::GetProcess<mcparticle::Flags, mcparticle::StatusCode>,
                  mcparticle::GetHepMCStatusCode<mcparticle::Flags, mcparticle::StatusCode>,
                  mcparticle::IsPhysicalPrimary<mcparticle::Flags>);

using ReducedA3MCTrack = ReducedA3MCTracks::iterator;

namespace reduceda3barreltracklabel
{
DECLARE_SOA_INDEX_COLUMN(ReducedA3MCTrack, reducedA3MCTrack); //!
DECLARE_SOA_COLUMN(McMask, mcMask, uint16_t);
} // namespace reduceda3barreltracklabel

// NOTE: MC labels. This table has one entry for each reconstructed track (joinable with the track tables)
//          The McParticleId points to the position of the MC truth track from the ReducedTracksMC table
DECLARE_SOA_TABLE(ReducedA3TracksBarrelLabels, "AOD", "REA3BARLA", //!
                  reduceda3barreltracklabel::ReducedA3MCTrackId, reduceda3barreltracklabel::McMask, reducedA3trackMC::McReducedFlags);

using ReducedA3TrackBarrelLabel = ReducedA3TracksBarrelLabels::iterator;

DECLARE_SOA_TABLE(ReducedA3TracksBarrel, "AOD", "REA3BARREL",
                  track::X, track::Alpha, track::IsWithinBeamPipe<track::X>,
                  track::Y, track::Z, track::Snp, track::Tgl, track::Signed1Pt,
                  track::Flags, track::ITSClusterMap, track::ITSChi2NCl,
                  reducedtrackalice3::IsReconstructed, reducedtrackalice3::NSiliconHits,
                  reducedtrackalice3::NTPCHits, track::Length, reducedtrack::DcaXY, reducedtrack::DcaZ,
                  track::IsPVContributor<track::Flags>);

// barrel collision information (joined with ReducedTracks) allowing to connect different tables (cross PWGs)
DECLARE_SOA_TABLE(ReducedA3TracksBarrelInfo, "AOD", "REA3BARRELINFO",
                  reducedtrackalice3::CollisionId, collision::PosX, collision::PosY, collision::PosZ, reducedtrackalice3::TrackId);

using ReducedA3Track = ReducedA3Tracks::iterator;
using ReducedA3TrackBarrel = ReducedA3TracksBarrel::iterator;
using ReducedA3TrackBarrelCov = ReducedA3TracksBarrelCov::iterator;
using ReducedA3TrackBarrelInfo = ReducedA3TracksBarrelInfo::iterator;

namespace reducedeventlabela3
{
DECLARE_SOA_INDEX_COLUMN(ReducedA3MCEvent, reducedA3MCEvent); //! MC collision
DECLARE_SOA_COLUMN(McMask, mcMask, uint16_t);                 //! Bit mask to indicate collision mismatches (bit ON means mismatch). Bit 15: indicates negative label
} // namespace reducedeventlabela3

DECLARE_SOA_TABLE(ReducedA3MCEventLabels, "AOD", "REA3MCCOLLBL", //! Table joined to the ReducedEvents table containing the MC index
                  reducedeventlabela3::ReducedA3MCEventId, reducedeventlabela3::McMask);

using ReducedA3MCEventLabel = ReducedA3MCEventLabels::iterator;

namespace reducedA3track_association
{
DECLARE_SOA_INDEX_COLUMN(ReducedA3Event, reducedA3event); //! ReducedEvent index
DECLARE_SOA_INDEX_COLUMN(ReducedA3Track, reducedA3track); //! ReducedTrack index
} // namespace reducedA3track_association

DECLARE_SOA_TABLE(ReducedA3TracksAssoc, "AOD", "REA3ASSOC", //! Table for reducedtrack-to-reducedcollision association
                  reducedA3track_association::ReducedA3EventId,
                  reducedA3track_association::ReducedA3TrackId);

DECLARE_SOA_TABLE(ReducedA3PIDTOF, "AOD", "REA3PIDTOF",
                  upgrade_tof::TOFEventTime,
                  upgrade_tof::TOFEventTimeErr,
                  upgrade_tof::NSigmaElectronInnerTOF,
                  upgrade_tof::NSigmaMuonInnerTOF,
                  upgrade_tof::NSigmaPionInnerTOF,
                  upgrade_tof::NSigmaKaonInnerTOF,
                  upgrade_tof::NSigmaProtonInnerTOF,
                  upgrade_tof::NSigmaDeuteronInnerTOF,
                  upgrade_tof::NSigmaTritonInnerTOF,
                  upgrade_tof::NSigmaHelium3InnerTOF,
                  upgrade_tof::NSigmaAlphaInnerTOF,
                  upgrade_tof::InnerTOFTrackTimeReco,
                  upgrade_tof::InnerTOFTrackLengthReco,
                  upgrade_tof::NSigmaElectronOuterTOF,
                  upgrade_tof::NSigmaMuonOuterTOF,
                  upgrade_tof::NSigmaPionOuterTOF,
                  upgrade_tof::NSigmaKaonOuterTOF,
                  upgrade_tof::NSigmaProtonOuterTOF,
                  upgrade_tof::NSigmaDeuteronOuterTOF,
                  upgrade_tof::NSigmaTritonOuterTOF,
                  upgrade_tof::NSigmaHelium3OuterTOF,
                  upgrade_tof::NSigmaAlphaOuterTOF,
                  upgrade_tof::OuterTOFTrackTimeReco,
                  upgrade_tof::OuterTOFTrackLengthReco,
                  upgrade_tof::NSigmaInnerTOF<upgrade_tof::NSigmaElectronInnerTOF,
                                              upgrade_tof::NSigmaMuonInnerTOF,
                                              upgrade_tof::NSigmaPionInnerTOF,
                                              upgrade_tof::NSigmaKaonInnerTOF,
                                              upgrade_tof::NSigmaProtonInnerTOF,
                                              upgrade_tof::NSigmaDeuteronInnerTOF,
                                              upgrade_tof::NSigmaTritonInnerTOF,
                                              upgrade_tof::NSigmaHelium3InnerTOF,
                                              upgrade_tof::NSigmaAlphaInnerTOF>,
                  upgrade_tof::NSigmaOuterTOF<upgrade_tof::NSigmaElectronOuterTOF,
                                              upgrade_tof::NSigmaMuonOuterTOF,
                                              upgrade_tof::NSigmaPionOuterTOF,
                                              upgrade_tof::NSigmaKaonOuterTOF,
                                              upgrade_tof::NSigmaProtonOuterTOF,
                                              upgrade_tof::NSigmaDeuteronOuterTOF,
                                              upgrade_tof::NSigmaTritonOuterTOF,
                                              upgrade_tof::NSigmaHelium3OuterTOF,
                                              upgrade_tof::NSigmaAlphaOuterTOF>);

DECLARE_SOA_TABLE(ReducedA3PIDRich, "AOD", "REA3PIDRICH",
                  upgrade_rich::NSigmaElectronRich,
                  upgrade_rich::NSigmaMuonRich,
                  upgrade_rich::NSigmaPionRich,
                  upgrade_rich::NSigmaKaonRich,
                  upgrade_rich::NSigmaProtonRich,
                  upgrade_rich::NSigmaDeuteronRich,
                  upgrade_rich::NSigmaTritonRich,
                  upgrade_rich::NSigmaHelium3Rich,
                  upgrade_rich::NSigmaAlphaRich,
                  upgrade_rich::NSigmaRich<upgrade_rich::NSigmaElectronRich,
                                           upgrade_rich::NSigmaMuonRich,
                                           upgrade_rich::NSigmaPionRich,
                                           upgrade_rich::NSigmaKaonRich,
                                           upgrade_rich::NSigmaProtonRich,
                                           upgrade_rich::NSigmaDeuteronRich,
                                           upgrade_rich::NSigmaTritonRich,
                                           upgrade_rich::NSigmaHelium3Rich,
                                           upgrade_rich::NSigmaAlphaRich>);

DECLARE_SOA_TABLE(ReducedA3PIDRichSignals, "AOD", "REA3PIDRICHSIG",
                  upgrade_rich::HasSig,
                  upgrade_rich::HasSigInGas,
                  upgrade_rich::HasSigEl,
                  upgrade_rich::HasSigMu,
                  upgrade_rich::HasSigPi,
                  upgrade_rich::HasSigKa,
                  upgrade_rich::HasSigPr,
                  upgrade_rich::HasSigDe,
                  upgrade_rich::HasSigTr,
                  upgrade_rich::HasSigHe3,
                  upgrade_rich::HasSigAl);

DECLARE_SOA_TABLE(ReducedA3PIDOT, "AOD", "REA3PIDOT",
                  upgrade::trk::TimeOverThresholdBarrel,
                  upgrade::trk::NSigmaTrkEl,
                  upgrade::trk::NSigmaTrkMu,
                  upgrade::trk::NSigmaTrkPi,
                  upgrade::trk::NSigmaTrkKa,
                  upgrade::trk::NSigmaTrkPr,
                  upgrade::trk::NSigmaTrkDe,
                  upgrade::trk::NSigmaTrkTr,
                  upgrade::trk::NSigmaTrkHe,
                  upgrade::trk::NSigmaTrkAl,
                  upgrade::trk::NSigmaTrk<upgrade::trk::NSigmaTrkEl,
                                          upgrade::trk::NSigmaTrkMu,
                                          upgrade::trk::NSigmaTrkPi,
                                          upgrade::trk::NSigmaTrkKa,
                                          upgrade::trk::NSigmaTrkPr,
                                          upgrade::trk::NSigmaTrkDe,
                                          upgrade::trk::NSigmaTrkTr,
                                          upgrade::trk::NSigmaTrkHe,
                                          upgrade::trk::NSigmaTrkAl>);

} // namespace o2::aod

#endif // PWGDQ_DATAMODEL_REDUCEDTABLESALICE3_H_
