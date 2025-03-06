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
#ifndef PWGLF_DATAMODEL_LFSTRANGENESSTABLES_H_
#define PWGLF_DATAMODEL_LFSTRANGENESSTABLES_H_

#include <cmath>
#include <vector>
#include "Framework/AnalysisDataModel.h"
#include "Common/Core/RecoDecay.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/McCollisionExtra.h"
#include "PWGLF/DataModel/EPCalibrationTables.h"
#include "PWGLF/DataModel/SPCalibrationTables.h"
#include "PWGUD/DataModel/UDTables.h"

namespace o2::aod
{

// for DF name follow-up and debug
namespace straorigin
{
DECLARE_SOA_COLUMN(DataframeID, dataframeID, uint64_t); //! Data frame ID (what is usually found in directory name in the AO2D.root, i.e.
} // namespace straorigin

DECLARE_SOA_TABLE(StraOrigins, "AOD", "STRAORIGIN", //! Table which contains the IDs of all dataframes merged into this dataframe
                  o2::soa::Index<>, straorigin::DataframeID);

namespace stracollision
{
DECLARE_SOA_DYNAMIC_COLUMN(IsUPC, isUPC, //! check whether this is a UPC or hadronic collision
                           [](int value) -> bool { return value <= 2 ? true : false; });
DECLARE_SOA_DYNAMIC_COLUMN(TotalFV0AmplitudeA, totalFV0AmplitudeA, //! get the total sum of the FV0 A amplitudes
                           [](float value) -> float { return value; });
DECLARE_SOA_DYNAMIC_COLUMN(TotalFT0AmplitudeA, totalFT0AmplitudeA, //! get the total sum of the FT0 A amplitudes
                           [](float value) -> float { return value; });
DECLARE_SOA_DYNAMIC_COLUMN(TotalFT0AmplitudeC, totalFT0AmplitudeC, //! get the total sum of the FT0 C amplitudes
                           [](float value) -> float { return value; });
DECLARE_SOA_DYNAMIC_COLUMN(TotalFDDAmplitudeA, totalFDDAmplitudeA, //! get the total sum of the FDD A amplitudes
                           [](float value) -> float { return value; });
DECLARE_SOA_DYNAMIC_COLUMN(TotalFDDAmplitudeC, totalFDDAmplitudeC, //! get the total sum of the FDD C amplitudes
                           [](float value) -> float { return value; });
DECLARE_SOA_DYNAMIC_COLUMN(EnergyCommonZNA, energyCommonZNA, //! get the total sum of the ZN A amplitudes
                           [](float value) -> float { return value; });
DECLARE_SOA_DYNAMIC_COLUMN(EnergyCommonZNC, energyCommonZNC, //! get the total sum of the ZN A amplitudes
                           [](float value) -> float { return value; });
} // namespace stracollision

//______________________________________________________
// Collision declarations for derived data analysis
// this is optional but will ensure full flexibility
// if required (for 2pc, etc)
DECLARE_SOA_TABLE(StraCollisions, "AOD", "STRACOLLISION", //! basic collision properties: position
                  o2::soa::Index<>, collision::PosX, collision::PosY, collision::PosZ);
DECLARE_SOA_TABLE(StraCents_000, "AOD", "STRACENTS", //! centrality percentiles
                  cent::CentFT0M, cent::CentFT0A,
                  cent::CentFT0C, cent::CentFV0A);
DECLARE_SOA_TABLE_VERSIONED(StraCents_001, "AOD", "STRACENTS", 1, //! centrality percentiles in Run 3
                            cent::CentFT0M, cent::CentFT0A,
                            cent::CentFT0C, cent::CentFV0A,
                            cent::CentFT0CVariant1, cent::CentMFT,
                            cent::CentNGlobal);

DECLARE_SOA_TABLE(StraCentsRun2, "AOD", "STRACENTSRUN2", //! centrality percentiles in Run 2
                  cent::CentRun2V0M, cent::CentRun2V0A,
                  cent::CentRun2SPDTracklets, cent::CentRun2SPDClusters);

// !!! DEPRECATED TABLE: StraRawCents_000 !!! All info in StraEvSels_001, in order to group all event characteristics in a unique table. Please use StraEvSels_001
DECLARE_SOA_TABLE(StraRawCents_000, "AOD", "STRARAWCENTS", //! debug information
                  mult::MultFT0A, mult::MultFT0C, mult::MultFV0A, mult::MultNTracksPVeta1);
// !!! DEPRECATED TABLE: StraRawCents_001 !!! All info in StraEvSels_001, in order to group all event characteristics in a unique table. Please use StraEvSels_001
DECLARE_SOA_TABLE_VERSIONED(StraRawCents_001, "AOD", "STRARAWCENTS", 1,     //! debug information
                            mult::MultFT0A, mult::MultFT0C, mult::MultFV0A, // FIT detectors
                            mult::MultNTracksPVeta1,                        // track multiplicities
                            mult::MultZNA, mult::MultZNC, mult::MultZEM1,   // ZDC signals
                            mult::MultZEM2, mult::MultZPA, mult::MultZPC);
// !!! DEPRECATED TABLE: StraRawCents_002 !!! All info in StraEvSels_001, in order to group all event characteristics in a unique table. Please use StraEvSels_001
DECLARE_SOA_TABLE_VERSIONED(StraRawCents_002, "AOD", "STRARAWCENTS", 2,            //! debug information
                            mult::MultFT0A, mult::MultFT0C, mult::MultFV0A,        // FIT detectors
                            mult::MultNTracksPVeta1,                               // track multiplicities with eta cut for INEL>0
                            mult::MultNTracksITSTPC,                               // track multiplicities, PV contribs, no eta cut
                            mult::MultAllTracksTPCOnly, mult::MultAllTracksITSTPC, // track multiplicities, all, no eta cut
                            mult::MultZNA, mult::MultZNC, mult::MultZEM1,          // ZDC signals
                            mult::MultZEM2, mult::MultZPA, mult::MultZPC);
// !!! DEPRECATED TABLE: StraRawCents_003 !!! All info in StraEvSels_001, in order to group all event characteristics in a unique table. Please use StraEvSels_001
DECLARE_SOA_TABLE_VERSIONED(StraRawCents_003, "AOD", "STRARAWCENTS", 3,     //! debug information
                            mult::MultFT0A, mult::MultFT0C, mult::MultFV0A, // FIT detectors
                            mult::MultNTracksPVeta1,                        // track multiplicities with eta cut for INEL>0
                            mult::MultPVTotalContributors,                  // number of PV contribs total
                            mult::MultNTracksGlobal,                        // global track multiplicities
                            mult::MultNTracksITSTPC,                        // track multiplicities, PV contribs, no eta cut
                            mult::MultAllTracksTPCOnly,                     // TPConly track multiplicities, all, no eta cut
                            mult::MultAllTracksITSTPC,                      // ITSTPC track multiplicities, all, no eta cut
                            mult::MultZNA, mult::MultZNC, mult::MultZEM1,   // ZDC signals
                            mult::MultZEM2, mult::MultZPA, mult::MultZPC);
// !!! DEPRECATED TABLE: StraRawCents_004 !!! All info in StraEvSels_001, in order to group all event characteristics in a unique table. Please use StraEvSels_001
DECLARE_SOA_TABLE_VERSIONED(StraRawCents_004, "AOD", "STRARAWCENTS", 4,     //! debug information
                            mult::MultFT0A, mult::MultFT0C, mult::MultFV0A, // FIT detectors
                            mult::MultNTracksPVeta1,                        // track multiplicities with eta cut for INEL>0
                            mult::MultPVTotalContributors,                  // number of PV contribs total
                            mult::MultNTracksGlobal,                        // global track multiplicities
                            mult::MultNTracksITSTPC,                        // track multiplicities, PV contribs, no eta cut
                            mult::MultAllTracksTPCOnly,                     // TPConly track multiplicities, all, no eta cut
                            mult::MultAllTracksITSTPC,                      // ITSTPC track multiplicities, all, no eta cut
                            mult::MultZNA, mult::MultZNC, mult::MultZEM1,   // ZDC signals
                            mult::MultZEM2, mult::MultZPA, mult::MultZPC,
                            evsel::NumTracksInTimeRange); // add occupancy as extra
DECLARE_SOA_TABLE(StraEvSels_000, "AOD", "STRAEVSELS",    //! event selection: sel8
                  evsel::Sel8, evsel::Selection);
DECLARE_SOA_TABLE_VERSIONED(StraEvSels_001, "AOD", "STRAEVSELS", 1,         //! debug information
                            evsel::Sel8, evsel::Selection,                  //! event selection: sel8
                            mult::MultFT0A, mult::MultFT0C, mult::MultFV0A, // FIT detectors
                            mult::MultFDDA, mult::MultFDDC,
                            mult::MultNTracksPVeta1,                      // track multiplicities with eta cut for INEL>0
                            mult::MultPVTotalContributors,                // number of PV contribs total
                            mult::MultNTracksGlobal,                      // global track multiplicities
                            mult::MultNTracksITSTPC,                      // track multiplicities, PV contribs, no eta cut
                            mult::MultAllTracksTPCOnly,                   // TPConly track multiplicities, all, no eta cut
                            mult::MultAllTracksITSTPC,                    // ITSTPC track multiplicities, all, no eta cut
                            mult::MultZNA, mult::MultZNC, mult::MultZEM1, // ZDC signals
                            mult::MultZEM2, mult::MultZPA, mult::MultZPC,
                            evsel::NumTracksInTimeRange,     // add occupancy as extra
                            udcollision::GapSide,            // UPC info: 0 for side A, 1 for side C, 2 for both sides, 3 neither A or C, 4 not enough or too many pv contributors
                            udcollision::TotalFT0AmplitudeA, // UPC info: re-assigned FT0-A amplitude, in case of SG event, from the most active bc
                            udcollision::TotalFT0AmplitudeC, // UPC info: re-assigned FT0-C amplitude, in case of SG event, from the most active bc
                            udcollision::TotalFV0AmplitudeA, // UPC info: re-assigned FV0-A amplitude, in case of SG event, from the most active bc
                            udcollision::TotalFDDAmplitudeA, // UPC info: re-assigned FDD-A amplitude, in case of SG event, from the most active bc
                            udcollision::TotalFDDAmplitudeC, // UPC info: re-assigned FDD-C amplitude, in case of SG event, from the most active bc
                            udzdc::EnergyCommonZNA,          // UPC info: re-assigned ZN-A amplitude, in case of SG event, from the most active bc
                            udzdc::EnergyCommonZNC,          // UPC info: re-assigned ZN-C amplitude, in case of SG event, from the most active bc
                            stracollision::IsUPC<udcollision::GapSide>);

DECLARE_SOA_TABLE_VERSIONED(StraEvSels_002, "AOD", "STRAEVSELS", 2,         //! debug information
                            evsel::Sel8, evsel::Selection,                  //! event selection: sel8
                            mult::MultFT0A, mult::MultFT0C, mult::MultFV0A, // FIT detectors
                            mult::MultFDDA, mult::MultFDDC,
                            mult::MultNTracksPVeta1,                      // track multiplicities with eta cut for INEL>0
                            mult::MultPVTotalContributors,                // number of PV contribs total
                            mult::MultNTracksGlobal,                      // global track multiplicities
                            mult::MultNTracksITSTPC,                      // track multiplicities, PV contribs, no eta cut
                            mult::MultAllTracksTPCOnly,                   // TPConly track multiplicities, all, no eta cut
                            mult::MultAllTracksITSTPC,                    // ITSTPC track multiplicities, all, no eta cut
                            mult::MultZNA, mult::MultZNC, mult::MultZEM1, // ZDC signals
                            mult::MultZEM2, mult::MultZPA, mult::MultZPC,
                            evsel::NumTracksInTimeRange,     // add occupancy in specified time interval by a number of tracks from nearby collisions
                            udcollision::GapSide,            // UPC info: 0 for side A, 1 for side C, 2 for both sides, 3 neither A or C, 4 not enough or too many pv contributors
                            udcollision::TotalFT0AmplitudeA, // UPC info: re-assigned FT0-A amplitude, in case of SG event, from the most active bc
                            udcollision::TotalFT0AmplitudeC, // UPC info: re-assigned FT0-C amplitude, in case of SG event, from the most active bc
                            udcollision::TotalFV0AmplitudeA, // UPC info: re-assigned FV0-A amplitude, in case of SG event, from the most active bc
                            udcollision::TotalFDDAmplitudeA, // UPC info: re-assigned FDD-A amplitude, in case of SG event, from the most active bc
                            udcollision::TotalFDDAmplitudeC, // UPC info: re-assigned FDD-C amplitude, in case of SG event, from the most active bc
                            udzdc::EnergyCommonZNA,          // UPC info: re-assigned ZN-A amplitude, in case of SG event, from the most active bc
                            udzdc::EnergyCommonZNC,          // UPC info: re-assigned ZN-C amplitude, in case of SG event, from the most active bc
                            collision::Flags,                // Contains Vertex::Flags, with most notably the UPCMode to know whether the vertex has been found using UPC settings
                            stracollision::IsUPC<udcollision::GapSide>);

DECLARE_SOA_TABLE_VERSIONED(StraEvSels_003, "AOD", "STRAEVSELS", 3,         //! debug information
                            evsel::Sel8, evsel::Selection,                  //! event selection: sel8
                            mult::MultFT0A, mult::MultFT0C, mult::MultFV0A, // FIT detectors
                            mult::MultFDDA, mult::MultFDDC,
                            mult::MultNTracksPVeta1,                      // track multiplicities with eta cut for INEL>0
                            mult::MultPVTotalContributors,                // number of PV contribs total
                            mult::MultNTracksGlobal,                      // global track multiplicities
                            mult::MultNTracksITSTPC,                      // track multiplicities, PV contribs, no eta cut
                            mult::MultAllTracksTPCOnly,                   // TPConly track multiplicities, all, no eta cut
                            mult::MultAllTracksITSTPC,                    // ITSTPC track multiplicities, all, no eta cut
                            mult::MultZNA, mult::MultZNC, mult::MultZEM1, // ZDC signals
                            mult::MultZEM2, mult::MultZPA, mult::MultZPC,
                            evsel::NumTracksInTimeRange,     // add occupancy in specified time interval by a number of tracks from nearby collisions
                            evsel::SumAmpFT0CInTimeRange,    // add occupancy in specified time interval by a sum of FT0C amplitudes from nearby collisions
                            udcollision::GapSide,            // UPC info: 0 for side A, 1 for side C, 2 for both sides, 3 neither A or C, 4 not enough or too many pv contributors
                            udcollision::TotalFT0AmplitudeA, // UPC info: re-assigned FT0-A amplitude, in case of SG event, from the most active bc
                            udcollision::TotalFT0AmplitudeC, // UPC info: re-assigned FT0-C amplitude, in case of SG event, from the most active bc
                            udcollision::TotalFV0AmplitudeA, // UPC info: re-assigned FV0-A amplitude, in case of SG event, from the most active bc
                            udcollision::TotalFDDAmplitudeA, // UPC info: re-assigned FDD-A amplitude, in case of SG event, from the most active bc
                            udcollision::TotalFDDAmplitudeC, // UPC info: re-assigned FDD-C amplitude, in case of SG event, from the most active bc
                            udzdc::EnergyCommonZNA,          // UPC info: re-assigned ZN-A amplitude, in case of SG event, from the most active bc
                            udzdc::EnergyCommonZNC,          // UPC info: re-assigned ZN-C amplitude, in case of SG event, from the most active bc
                            collision::Flags,                // Contains Vertex::Flags, with most notably the UPCMode to know whether the vertex has been found using UPC settings
                            stracollision::IsUPC<udcollision::GapSide>);

DECLARE_SOA_TABLE_VERSIONED(StraEvSels_004, "AOD", "STRAEVSELS", 4,         //! debug information
                            evsel::Sel8, evsel::Selection,                  //! event selection: sel8
                            mult::MultFT0A, mult::MultFT0C, mult::MultFV0A, // FIT detectors
                            mult::MultFDDA, mult::MultFDDC,
                            mult::MultNTracksPVeta1,                      // track multiplicities with eta cut for INEL>0
                            mult::MultPVTotalContributors,                // number of PV contribs total
                            mult::MultNTracksGlobal,                      // global track multiplicities
                            mult::MultNTracksITSTPC,                      // track multiplicities, PV contribs, no eta cut
                            mult::MultAllTracksTPCOnly,                   // TPConly track multiplicities, all, no eta cut
                            mult::MultAllTracksITSTPC,                    // ITSTPC track multiplicities, all, no eta cut
                            mult::MultZNA, mult::MultZNC, mult::MultZEM1, // ZDC signals
                            mult::MultZEM2, mult::MultZPA, mult::MultZPC,
                            evsel::NumTracksInTimeRange,     // add occupancy in specified time interval by a number of tracks from nearby collisions
                            evsel::SumAmpFT0CInTimeRange,    // add occupancy in specified time interval by a sum of FT0C amplitudes from nearby collisions
                            udcollision::GapSide,            // UPC info: 0 for side A, 1 for side C, 2 for both sides, 3 neither A or C, 4 not enough or too many pv contributors
                            udcollision::TotalFT0AmplitudeA, // UPC info: re-assigned FT0-A amplitude, in case of SG event, from the most active bc
                            udcollision::TotalFT0AmplitudeC, // UPC info: re-assigned FT0-C amplitude, in case of SG event, from the most active bc
                            udcollision::TotalFV0AmplitudeA, // UPC info: re-assigned FV0-A amplitude, in case of SG event, from the most active bc
                            udcollision::TotalFDDAmplitudeA, // UPC info: re-assigned FDD-A amplitude, in case of SG event, from the most active bc
                            udcollision::TotalFDDAmplitudeC, // UPC info: re-assigned FDD-C amplitude, in case of SG event, from the most active bc
                            udzdc::EnergyCommonZNA,          // UPC info: re-assigned ZN-A amplitude, in case of SG event, from the most active bc
                            udzdc::EnergyCommonZNC,          // UPC info: re-assigned ZN-C amplitude, in case of SG event, from the most active bc

                            collision::Flags, // Contains Vertex::Flags, with most notably the UPCMode to know whether the vertex has been found using UPC settings
                            evsel::Alias,     // trigger aliases (e.g. kTVXinTRD for v2)

                            // Dynamic columns for manipulating information
                            // stracollision::TotalFV0AmplitudeA<mult::MultFV0A>,
                            // stracollision::TotalFT0AmplitudeA<mult::MultFT0A>,
                            // stracollision::TotalFT0AmplitudeC<mult::MultFT0C>,
                            // stracollision::TotalFDDAmplitudeA<mult::MultFDDA>,
                            // stracollision::TotalFDDAmplitudeC<mult::MultFDDC>,
                            // stracollision::EnergyCommonZNA<mult::MultZNA>,
                            // stracollision::EnergyCommonZNC<mult::MultZNC>,
                            stracollision::IsUPC<udcollision::GapSide>);

DECLARE_SOA_TABLE(StraEvSelsRun2, "AOD", "STRAEVSELSRUN2",        //! debug information
                  evsel::Sel8, evsel::Sel7, evsel::Selection,     //! event selection: sel8
                  mult::MultFT0A, mult::MultFT0C, mult::MultFV0A, // FIT detectors
                  mult::MultFDDA, mult::MultFDDC,
                  mult::MultNTracksPVeta1,                      // track multiplicities with eta cut for INEL>0
                  mult::MultPVTotalContributors,                // number of PV contribs total
                  mult::MultNTracksGlobal,                      // global track multiplicities
                  mult::MultNTracksITSTPC,                      // track multiplicities, PV contribs, no eta cut
                  mult::MultAllTracksTPCOnly,                   // TPConly track multiplicities, all, no eta cut
                  mult::MultAllTracksITSTPC,                    // ITSTPC track multiplicities, all, no eta cut
                  mult::MultZNA, mult::MultZNC, mult::MultZEM1, // ZDC signals
                  mult::MultZEM2, mult::MultZPA, mult::MultZPC,
                  evsel::Alias); // trigger aliases (e.g. kTVXinTRD for v2)

DECLARE_SOA_TABLE(StraFT0AQVs, "AOD", "STRAFT0AQVS", //! t0a Qvec
                  qvec::QvecFT0ARe, qvec::QvecFT0AIm, qvec::SumAmplFT0A);
DECLARE_SOA_TABLE(StraFT0CQVs, "AOD", "STRAFT0CQVS", //! t0c Qvec
                  qvec::QvecFT0CRe, qvec::QvecFT0CIm, qvec::SumAmplFT0C);
DECLARE_SOA_TABLE(StraFT0MQVs, "AOD", "STRAFT0MQVS", //! t0m Qvec
                  qvec::QvecFT0MRe, qvec::QvecFT0MIm, qvec::SumAmplFT0M);
DECLARE_SOA_TABLE(StraFV0AQVs, "AOD", "STRAFV0AQVS", //! v0a Qvec
                  qvec::QvecFV0ARe, qvec::QvecFV0AIm, qvec::SumAmplFV0A);
DECLARE_SOA_TABLE(StraTPCQVs, "AOD", "STRATPCQVS", //! tpc Qvec
                  qvec::QvecBNegRe, qvec::QvecBNegIm, epcalibrationtable::QTPCL,
                  qvec::QvecBPosRe, qvec::QvecBPosIm, epcalibrationtable::QTPCR);
DECLARE_SOA_TABLE(StraFT0CQVsEv, "AOD", "STRAFT0CQVSEv", //! events used to compute t0c Qvec
                  epcalibrationtable::TriggerEventEP);
DECLARE_SOA_TABLE(StraZDCSP, "AOD", "STRAZDCSP", //! ZDC SP information
                  spcalibrationtable::TriggerEventSP,
                  spcalibrationtable::PsiZDCA, spcalibrationtable::PsiZDCC, spcalibrationtable::QXZDCA, spcalibrationtable::QXZDCC, spcalibrationtable::QYZDCA, spcalibrationtable::QYZDCC);
DECLARE_SOA_TABLE(StraStamps_000, "AOD", "STRASTAMPS", //! information for ID-ing mag field if needed
                  bc::RunNumber, timestamp::Timestamp);
DECLARE_SOA_TABLE_VERSIONED(StraStamps_001, "AOD", "STRASTAMPS", 1, //! information for ID-ing mag field if needed
                            bc::RunNumber, timestamp::Timestamp, bc::GlobalBC);

using StraRawCents = StraRawCents_004;
using StraCents = StraCents_001;
using StraEvSels = StraEvSels_004;
using StraStamps = StraStamps_001;
using StraCollision = StraCollisions::iterator;
using StraCent = StraCents_001::iterator;

namespace stramccollision
{
DECLARE_SOA_COLUMN(TotalMultMCParticles, totalMultMCParticles, int); //! total number of MC particles in a generated collision
} // namespace stramccollision

//______________________________________________________
// for correlating information with MC
// also allows for collision association cross-checks
DECLARE_SOA_TABLE(StraMCCollisions_000, "AOD", "STRAMCCOLLISION", //! MC collision properties
                  o2::soa::Index<>, mccollision::PosX, mccollision::PosY, mccollision::PosZ,
                  mccollision::ImpactParameter);
DECLARE_SOA_TABLE_VERSIONED(StraMCCollisions_001, "AOD", "STRAMCCOLLISION", 1, //! debug information
                            o2::soa::Index<>, mccollision::PosX, mccollision::PosY, mccollision::PosZ,
                            mccollision::ImpactParameter, mccollision::EventPlaneAngle);
using StraMCCollisions = StraMCCollisions_001;
using StraMCCollision = StraMCCollisions::iterator;

DECLARE_SOA_TABLE(StraMCCollMults_000, "AOD", "STRAMCCOLLMULTS", //! MC collision multiplicities
                  mult::MultMCFT0A, mult::MultMCFT0C, mult::MultMCNParticlesEta05, mult::MultMCNParticlesEta08, mult::MultMCNParticlesEta10, o2::soa::Marker<2>);
DECLARE_SOA_TABLE_VERSIONED(StraMCCollMults_001, "AOD", "STRAMCCOLLMULTS", 1, //! MC collision multiplicities
                            mult::MultMCFT0A, mult::MultMCFT0C, mult::MultMCNParticlesEta05, mult::MultMCNParticlesEta08, mult::MultMCNParticlesEta10, stramccollision::TotalMultMCParticles);

using StraMCCollMults = StraMCCollMults_001;

namespace dautrack
{
//______________________________________________________
// Daughter track declarations for derived data analysis
// These definitions are for the first version of the table
// The latest version will inherit most properties from TracksExtra
DECLARE_SOA_COLUMN(ITSChi2PerNcl, itsChi2PerNcl, float);        //! ITS chi2 per N cluster
DECLARE_SOA_COLUMN(DetectorMap, detectorMap, uint8_t);          //! detector map for reference (see DetectorMapEnum)
DECLARE_SOA_COLUMN(ITSClusterSizes, itsClusterSizes, uint32_t); //! ITS cluster sizes per layer
DECLARE_SOA_COLUMN(TPCClusters, tpcClusters, uint8_t);          //! N TPC clusters
DECLARE_SOA_COLUMN(TPCCrossedRows, tpcCrossedRows, uint8_t);    //! N TPC clusters

//______________________________________________________
// Daughter track MC information
DECLARE_SOA_COLUMN(ParticleMCId, particleMCId, int); //! particle MC Id

//______________________________________________________
// for extras: replicated here to ensure ease of manipulating the ITS information
// directly from the V0 extras table in simple ways for derived data as well
DECLARE_SOA_DYNAMIC_COLUMN(ITSClusterMap, itsClusterMap, //! ITS cluster map, one bit per layer, starting from the innermost
                           [](uint32_t itsClusterSizes) -> uint8_t {
                             uint8_t clmap = 0;
                             for (unsigned int layer = 0; layer < 7; layer++) {
                               if ((itsClusterSizes >> (layer * 4)) & 0xf) {
                                 clmap |= (1 << layer);
                               }
                             }
                             return clmap;
                           });
DECLARE_SOA_DYNAMIC_COLUMN(ITSNCls, itsNCls, //! Number of ITS clusters
                           [](uint32_t itsClusterSizes) -> uint8_t {
                             uint8_t itsNcls = 0;
                             for (int layer = 0; layer < 7; layer++) {
                               if ((itsClusterSizes >> (layer * 4)) & 0xf)
                                 itsNcls++;
                             }
                             return itsNcls;
                           });
DECLARE_SOA_DYNAMIC_COLUMN(HasITS, hasITS, //! Flag to check if track has a ITS match
                           [](uint8_t detectorMap) -> bool { return detectorMap & o2::aod::track::ITS; });
DECLARE_SOA_DYNAMIC_COLUMN(HasTPC, hasTPC, //! Flag to check if track has a TPC match
                           [](uint8_t detectorMap) -> bool { return detectorMap & o2::aod::track::TPC; });
DECLARE_SOA_DYNAMIC_COLUMN(HasTRD, hasTRD, //! Flag to check if track has a TRD match
                           [](uint8_t detectorMap) -> bool { return detectorMap & o2::aod::track::TRD; });
DECLARE_SOA_DYNAMIC_COLUMN(HasTOF, hasTOF, //! Flag to check if track has a TOF measurement
                           [](uint8_t detectorMap) -> bool { return detectorMap & o2::aod::track::TOF; });
DECLARE_SOA_DYNAMIC_COLUMN(HasITSTracker, hasITSTracker, //! Flag to check if track is from ITS tracker
                           [](uint8_t detectorMap, float itsChi2PerNcl) -> bool { return (detectorMap & o2::aod::track::ITS) ? (itsChi2PerNcl > -1e-3f) : false; });
DECLARE_SOA_DYNAMIC_COLUMN(HasITSAfterburner, hasITSAfterburner, //! Flag to check if track is from ITS AB
                           [](uint8_t detectorMap, float itsChi2PerNcl) -> bool { return (detectorMap & o2::aod::track::ITS) ? (itsChi2PerNcl < -1e-3f) : false; });

// sub-namespace for compatibility purposes
namespace compatibility
{                                                    // adds dynamics that ensure full backwards compatibility with previous getters
DECLARE_SOA_DYNAMIC_COLUMN(TPCClusters, tpcClusters, //! number of TPC clusters
                           [](uint8_t tpcNClsFindable, int8_t tpcNClsFindableMinusFound) -> int16_t { return (int16_t)tpcNClsFindable - tpcNClsFindableMinusFound; });
DECLARE_SOA_DYNAMIC_COLUMN(TPCCrossedRows, tpcCrossedRows, //! Number of crossed TPC Rows
                           [](uint8_t tpcNClsFindable, int8_t TPCNClsFindableMinusCrossedRows) -> int16_t { return (int16_t)tpcNClsFindable - TPCNClsFindableMinusCrossedRows; });
DECLARE_SOA_DYNAMIC_COLUMN(ITSChi2PerNcl, itsChi2PerNcl, //! simple equivalent return
                           [](float itsChi2NCl) -> float { return (float)itsChi2NCl; });
} // namespace compatibility

} // namespace dautrack

DECLARE_SOA_TABLE(DauTrackExtras_000, "AOD", "DAUTRACKEXTRA", //! detector properties of decay daughters
                  dautrack::DetectorMap, dautrack::ITSClusterSizes,
                  dautrack::TPCClusters, dautrack::TPCCrossedRows,

                  // Dynamic columns for manipulating information
                  dautrack::ITSClusterMap<dautrack::ITSClusterSizes>,
                  dautrack::ITSNCls<dautrack::ITSClusterSizes>,
                  dautrack::HasITS<dautrack::DetectorMap>,
                  dautrack::HasTPC<dautrack::DetectorMap>,
                  dautrack::HasTRD<dautrack::DetectorMap>,
                  dautrack::HasTOF<dautrack::DetectorMap>,
                  dautrack::HasITSTracker<dautrack::DetectorMap, dautrack::ITSChi2PerNcl>,
                  dautrack::HasITSAfterburner<dautrack::DetectorMap, dautrack::ITSChi2PerNcl>);

DECLARE_SOA_TABLE_VERSIONED(DauTrackExtras_001, "AOD", "DAUTRACKEXTRA", 1, //! detector properties of decay daughters
                            dautrack::ITSChi2PerNcl,
                            dautrack::DetectorMap, dautrack::ITSClusterSizes,
                            dautrack::TPCClusters, dautrack::TPCCrossedRows,

                            // Dynamic columns for manipulating information
                            dautrack::ITSClusterMap<dautrack::ITSClusterSizes>,
                            dautrack::ITSNCls<dautrack::ITSClusterSizes>,
                            dautrack::HasITS<dautrack::DetectorMap>,
                            dautrack::HasTPC<dautrack::DetectorMap>,
                            dautrack::HasTRD<dautrack::DetectorMap>,
                            dautrack::HasTOF<dautrack::DetectorMap>,
                            dautrack::HasITSTracker<dautrack::DetectorMap, dautrack::ITSChi2PerNcl>,
                            dautrack::HasITSAfterburner<dautrack::DetectorMap, dautrack::ITSChi2PerNcl>);

DECLARE_SOA_TABLE_VERSIONED(DauTrackExtras_002, "AOD", "DAUTRACKEXTRA", 2, //! detector properties of decay daughters
                            track::ITSChi2NCl,
                            dautrack::DetectorMap, // here we don´t save everything so we simplify this
                            track::ITSClusterSizes,
                            track::TPCNClsFindable,
                            track::TPCNClsFindableMinusFound,
                            track::TPCNClsFindableMinusCrossedRows,

                            // Dynamics for ITS matching TracksExtra
                            track::v001::ITSNClsInnerBarrel<track::ITSClusterSizes>,
                            track::v001::ITSClsSizeInLayer<track::ITSClusterSizes>,
                            track::v001::ITSClusterMap<track::ITSClusterSizes>,
                            track::v001::ITSNCls<track::ITSClusterSizes>,
                            track::v001::IsITSAfterburner<track::v001::DetectorMap, track::ITSChi2NCl>,
                            /*compatibility*/ dautrack::HasITSTracker<dautrack::DetectorMap, track::ITSChi2NCl>,
                            /*compatibility*/ dautrack::HasITSAfterburner<dautrack::DetectorMap, track::ITSChi2NCl>,

                            // dynamics for TPC tracking properties matching main data model
                            track::TPCCrossedRowsOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                            track::TPCFoundOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                            track::TPCNClsFound<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                            track::TPCNClsCrossedRows<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                            /*compatibility*/ dautrack::compatibility::TPCClusters<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                            /*compatibility*/ dautrack::compatibility::TPCCrossedRows<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                            /*compatibility*/ dautrack::compatibility::ITSChi2PerNcl<track::ITSChi2NCl>,

                            // dynamics to identify detectors
                            dautrack::HasITS<dautrack::DetectorMap>,
                            dautrack::HasTPC<dautrack::DetectorMap>,
                            dautrack::HasTRD<dautrack::DetectorMap>,
                            dautrack::HasTOF<dautrack::DetectorMap>);

DECLARE_SOA_TABLE_VERSIONED(DauTrackExtras_003, "AOD", "DAUTRACKEXTRA", 3, //! detector properties of decay daughters
                            track::ITSChi2NCl,
                            track::TPCChi2NCl,
                            dautrack::DetectorMap, // here we don´t save everything so we simplify this
                            track::ITSClusterSizes,
                            track::TPCNClsFindable,
                            track::TPCNClsFindableMinusFound,
                            track::TPCNClsFindableMinusCrossedRows,
                            track::TPCNClsShared,

                            // Dynamics for ITS matching TracksExtra
                            track::v001::ITSNClsInnerBarrel<track::ITSClusterSizes>,
                            track::v001::ITSClsSizeInLayer<track::ITSClusterSizes>,
                            track::v001::ITSClusterMap<track::ITSClusterSizes>,
                            track::v001::ITSNCls<track::ITSClusterSizes>,
                            track::v001::IsITSAfterburner<track::v001::DetectorMap, track::ITSChi2NCl>,
                            /*compatibility*/ dautrack::HasITSTracker<dautrack::DetectorMap, track::ITSChi2NCl>,
                            /*compatibility*/ dautrack::HasITSAfterburner<dautrack::DetectorMap, track::ITSChi2NCl>,

                            // dynamics for TPC tracking properties matching main data model
                            track::TPCCrossedRowsOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                            track::TPCFoundOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                            track::TPCNClsFound<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                            track::TPCNClsCrossedRows<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                            track::TPCFractionSharedCls<track::TPCNClsShared, track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                            /*compatibility*/ dautrack::compatibility::TPCClusters<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                            /*compatibility*/ dautrack::compatibility::TPCCrossedRows<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                            /*compatibility*/ dautrack::compatibility::ITSChi2PerNcl<track::ITSChi2NCl>,

                            // dynamics to identify detectors
                            dautrack::HasITS<dautrack::DetectorMap>,
                            dautrack::HasTPC<dautrack::DetectorMap>,
                            dautrack::HasTRD<dautrack::DetectorMap>,
                            dautrack::HasTOF<dautrack::DetectorMap>);

DECLARE_SOA_TABLE(DauTrackMCIds, "AOD", "DAUTRACKMCID", // index table when using AO2Ds
                  dautrack::ParticleMCId);

using DauTrackExtras = DauTrackExtras_003;
using DauTrackExtra = DauTrackExtras::iterator;

namespace motherParticle
{
DECLARE_SOA_COLUMN(Px, px, float);                              //! px
DECLARE_SOA_COLUMN(Py, py, float);                              //! py
DECLARE_SOA_COLUMN(Pz, pz, float);                              //! pz
DECLARE_SOA_COLUMN(PDGCode, pdgCode, int);                      //! pdg code
DECLARE_SOA_COLUMN(IsPhysicalPrimary, isPhysicalPrimary, bool); //! primary criterion
} // namespace motherParticle

DECLARE_SOA_TABLE(MotherMCParts, "AOD", "MOTHERMCPART", //! mother MC information, abbreviated name due to size limit
                  motherParticle::Px, motherParticle::Py, motherParticle::Pz,
                  motherParticle::PDGCode, motherParticle::IsPhysicalPrimary);

using MotherMCPart = MotherMCParts::iterator;

namespace v0data
{
//______________________________________________________
// REGULAR COLUMNS FOR INDEXING
DECLARE_SOA_INDEX_COLUMN_FULL(PosTrack, posTrack, int, Tracks, "_Pos"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(NegTrack, negTrack, int, Tracks, "_Neg"); //!
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                         //!
DECLARE_SOA_INDEX_COLUMN(V0, v0);                                       //!
// FOR DERIVED
DECLARE_SOA_INDEX_COLUMN_FULL(PosTrackExtra, posTrackExtra, int, DauTrackExtras, "_PosExtra"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(NegTrackExtra, negTrackExtra, int, DauTrackExtras, "_NegExtra"); //!
DECLARE_SOA_INDEX_COLUMN(StraCollision, straCollision);                                        //!
DECLARE_SOA_INDEX_COLUMN(StraMCCollision, straMCCollision);                                    //!
DECLARE_SOA_INDEX_COLUMN(MotherMCPart, motherMCPart);                                          //!
} // namespace v0data
//______________________________________________________
// intermezzo: define StraCollRefs, then carry on
DECLARE_SOA_TABLE(StraCollLabels, "AOD", "STRACOLLLABELS", //! optional table to refer back to a MC collision
                  o2::soa::Index<>, v0data::StraMCCollisionId, o2::soa::Marker<1>);
//______________________________________________________
namespace v0data
{
//______________________________________________________
// REGULAR COLUMNS FOR V0CORES
// General V0 properties: position, momentum
DECLARE_SOA_COLUMN(PosX, posX, float);   //! positive track X at min
DECLARE_SOA_COLUMN(NegX, negX, float);   //! negative track X at min
DECLARE_SOA_COLUMN(PxPos, pxpos, float); //! positive track px at min
DECLARE_SOA_COLUMN(PyPos, pypos, float); //! positive track py at min
DECLARE_SOA_COLUMN(PzPos, pzpos, float); //! positive track pz at min
DECLARE_SOA_COLUMN(PxNeg, pxneg, float); //! negative track px at min
DECLARE_SOA_COLUMN(PyNeg, pyneg, float); //! negative track py at min
DECLARE_SOA_COLUMN(PzNeg, pzneg, float); //! negative track pz at min
DECLARE_SOA_COLUMN(X, x, float);         //! decay position X
DECLARE_SOA_COLUMN(Y, y, float);         //! decay position Y
DECLARE_SOA_COLUMN(Z, z, float);         //! decay position Z

// decay daughter positions for refit studies (specific purpose)
DECLARE_SOA_COLUMN(XPosAtDCA, xPosAtDCA, float); //! decay position X
DECLARE_SOA_COLUMN(YPosAtDCA, yPosAtDCA, float); //! decay position Y
DECLARE_SOA_COLUMN(ZPosAtDCA, zPosAtDCA, float); //! decay position Z
DECLARE_SOA_COLUMN(XNegAtDCA, xNegAtDCA, float); //! decay position X
DECLARE_SOA_COLUMN(YNegAtDCA, yNegAtDCA, float); //! decay position Y
DECLARE_SOA_COLUMN(ZNegAtDCA, zNegAtDCA, float); //! decay position Z
DECLARE_SOA_COLUMN(XPosAtIU, xPosAtIU, float);   //! decay position X
DECLARE_SOA_COLUMN(YPosAtIU, yPosAtIU, float);   //! decay position Y
DECLARE_SOA_COLUMN(ZPosAtIU, zPosAtIU, float);   //! decay position Z
DECLARE_SOA_COLUMN(XNegAtIU, xNegAtIU, float);   //! decay position X
DECLARE_SOA_COLUMN(YNegAtIU, yNegAtIU, float);   //! decay position Y
DECLARE_SOA_COLUMN(ZNegAtIU, zNegAtIU, float);   //! decay position Z

// ivanov scaling
DECLARE_SOA_COLUMN(IvanovMap, ivanovMap, int); //! coded downscale bits

// Saved from finding: DCAs
DECLARE_SOA_COLUMN(DCAV0Daughters, dcaV0daughters, float); //! DCA between V0 daughters
DECLARE_SOA_COLUMN(DCAPosToPV, dcapostopv, float);         //! DCA positive prong to PV
DECLARE_SOA_COLUMN(DCANegToPV, dcanegtopv, float);         //! DCA negative prong to PV
DECLARE_SOA_COLUMN(V0CosPA, v0cosPA, float);               //! V0 CosPA
DECLARE_SOA_COLUMN(DCAV0ToPV, dcav0topv, float);           //! DCA V0 to PV

// Type of V0 from the svertexer (photon, regular, from cascade)
DECLARE_SOA_COLUMN(V0Type, v0Type, uint8_t); //! type of V0. 0: built solely for cascades (does not pass standard V0 cuts), 1: standard 2, 3: photon-like with TPC-only use. Regular analysis should always use type 1.

// Saved from finding: covariance matrix of parent track (on request)
DECLARE_SOA_COLUMN(PositionCovMat, positionCovMat, float[6]);  //! covariance matrix elements
DECLARE_SOA_COLUMN(MomentumCovMat, momentumCovMat, float[6]);  //! covariance matrix elements
DECLARE_SOA_COLUMN(CovMatPosDau, covMatPosDau, float[21]);     //! covariance matrix elements positive daughter track
DECLARE_SOA_COLUMN(CovMatNegDau, covMatNegDau, float[21]);     //! covariance matrix elements negative daughter track
DECLARE_SOA_COLUMN(CovMatPosDauIU, covMatPosDauIU, float[21]); //! covariance matrix elements positive daughter track
DECLARE_SOA_COLUMN(CovMatNegDauIU, covMatNegDauIU, float[21]); //! covariance matrix elements negative daughter track

// Saved from KF particle fit for specic table
DECLARE_SOA_COLUMN(KFV0Chi2, kfV0Chi2, float); //!

//______________________________________________________
// REGULAR COLUMNS FOR V0MCCORES
DECLARE_SOA_COLUMN(ParticleIdMC, particleIdMC, int);            //! V0 Particle ID
DECLARE_SOA_COLUMN(PDGCode, pdgCode, int);                      //! V0 PDG Code
DECLARE_SOA_COLUMN(PDGCodeMother, pdgCodeMother, int);          //! V0 mother PDG code (for feeddown)
DECLARE_SOA_COLUMN(PDGCodePositive, pdgCodePositive, int);      //! V0 positive prong PDG code
DECLARE_SOA_COLUMN(PDGCodeNegative, pdgCodeNegative, int);      //! V0 negative prong PDG code
DECLARE_SOA_COLUMN(IsPhysicalPrimary, isPhysicalPrimary, bool); //! is V0 physical primary
DECLARE_SOA_COLUMN(XMC, xMC, float);                            //! V0 decay position X (cm)
DECLARE_SOA_COLUMN(YMC, yMC, float);                            //! V0 decay position Y (cm)
DECLARE_SOA_COLUMN(ZMC, zMC, float);                            //! V0 decay position Z (cm)
DECLARE_SOA_COLUMN(PxPosMC, pxPosMC, float);                    //! V0 positive daughter px (GeV/c)
DECLARE_SOA_COLUMN(PyPosMC, pyPosMC, float);                    //! V0 positive daughter py (GeV/c)
DECLARE_SOA_COLUMN(PzPosMC, pzPosMC, float);                    //! V0 positive daughter pz (GeV/c)
DECLARE_SOA_COLUMN(PxNegMC, pxNegMC, float);                    //! V0 positive daughter px (GeV/c)
DECLARE_SOA_COLUMN(PyNegMC, pyNegMC, float);                    //! V0 positive daughter py (GeV/c)
DECLARE_SOA_COLUMN(PzNegMC, pzNegMC, float);                    //! V0 positive daughter pz (GeV/c)
DECLARE_SOA_COLUMN(PxMC, pxMC, float);                          //! V0 px (GeV/c)
DECLARE_SOA_COLUMN(PyMC, pyMC, float);                          //! V0 py (GeV/c)
DECLARE_SOA_COLUMN(PzMC, pzMC, float);                          //! V0 pz (GeV/c)

//______________________________________________________
// Binned content for generated particles: derived data
DECLARE_SOA_COLUMN(GeneratedK0Short, generatedK0Short, std::vector<uint32_t>);       //! K0Short binned generated data
DECLARE_SOA_COLUMN(GeneratedLambda, generatedLambda, std::vector<uint32_t>);         //! Lambda binned generated data
DECLARE_SOA_COLUMN(GeneratedAntiLambda, generatedAntiLambda, std::vector<uint32_t>); //! AntiLambda binned generated data

//______________________________________________________
// EXPRESSION COLUMNS
DECLARE_SOA_EXPRESSION_COLUMN(Px, px, //! V0 px
                              float, 1.f * aod::v0data::pxpos + 1.f * aod::v0data::pxneg);
DECLARE_SOA_EXPRESSION_COLUMN(Py, py, //! V0 py
                              float, 1.f * aod::v0data::pypos + 1.f * aod::v0data::pyneg);
DECLARE_SOA_EXPRESSION_COLUMN(Pz, pz, //! V0 pz
                              float, 1.f * aod::v0data::pzpos + 1.f * aod::v0data::pzneg);
DECLARE_SOA_EXPRESSION_COLUMN(Pt, pt, float, //! Transverse momentum in GeV/c
                              nsqrt((1.f * aod::v0data::pxpos + 1.f * aod::v0data::pxneg) *
                                      (1.f * aod::v0data::pxpos + 1.f * aod::v0data::pxneg) +
                                    (1.f * aod::v0data::pypos + 1.f * aod::v0data::pyneg) * (1.f * aod::v0data::pypos + 1.f * aod::v0data::pyneg)));
DECLARE_SOA_EXPRESSION_COLUMN(P, p, float, //! Total momentum in GeV/c
                              nsqrt((1.f * aod::v0data::pxpos + 1.f * aod::v0data::pxneg) *
                                      (1.f * aod::v0data::pxpos + 1.f * aod::v0data::pxneg) +
                                    (1.f * aod::v0data::pypos + 1.f * aod::v0data::pyneg) * (1.f * aod::v0data::pypos + 1.f * aod::v0data::pyneg) +
                                    (1.f * aod::v0data::pzpos + 1.f * aod::v0data::pzneg) * (1.f * aod::v0data::pzpos + 1.f * aod::v0data::pzneg)));
DECLARE_SOA_EXPRESSION_COLUMN(Phi, phi, float, //! Phi in the range [0, 2pi)
                              o2::constants::math::PI + natan2(-1.0f * (1.f * aod::v0data::pypos + 1.f * aod::v0data::pyneg), -1.0f * (1.f * aod::v0data::pxpos + 1.f * aod::v0data::pxneg)));
DECLARE_SOA_EXPRESSION_COLUMN(Eta, eta, float, //! Pseudorapidity, conditionally defined to avoid FPEs
                              ifnode((nsqrt((1.f * aod::v0data::pxpos + 1.f * aod::v0data::pxneg) * (1.f * aod::v0data::pxpos + 1.f * aod::v0data::pxneg) +
                                            (1.f * aod::v0data::pypos + 1.f * aod::v0data::pyneg) * (1.f * aod::v0data::pypos + 1.f * aod::v0data::pyneg) +
                                            (1.f * aod::v0data::pzpos + 1.f * aod::v0data::pzneg) * (1.f * aod::v0data::pzpos + 1.f * aod::v0data::pzneg)) -
                                      (1.f * aod::v0data::pzpos + 1.f * aod::v0data::pzneg)) < static_cast<float>(1e-7),
                                     ifnode((1.f * aod::v0data::pzpos + 1.f * aod::v0data::pzneg) < 0.f, -100.f, 100.f),
                                     0.5f * nlog((nsqrt((1.f * aod::v0data::pxpos + 1.f * aod::v0data::pxneg) * (1.f * aod::v0data::pxpos + 1.f * aod::v0data::pxneg) +
                                                        (1.f * aod::v0data::pypos + 1.f * aod::v0data::pyneg) * (1.f * aod::v0data::pypos + 1.f * aod::v0data::pyneg) +
                                                        (1.f * aod::v0data::pzpos + 1.f * aod::v0data::pzneg) * (1.f * aod::v0data::pzpos + 1.f * aod::v0data::pzneg)) +
                                                  (1.f * aod::v0data::pzpos + 1.f * aod::v0data::pzneg)) /
                                                 (nsqrt((1.f * aod::v0data::pxpos + 1.f * aod::v0data::pxneg) * (1.f * aod::v0data::pxpos + 1.f * aod::v0data::pxneg) +
                                                        (1.f * aod::v0data::pypos + 1.f * aod::v0data::pyneg) * (1.f * aod::v0data::pypos + 1.f * aod::v0data::pyneg) +
                                                        (1.f * aod::v0data::pzpos + 1.f * aod::v0data::pzneg) * (1.f * aod::v0data::pzpos + 1.f * aod::v0data::pzneg)) -
                                                  (1.f * aod::v0data::pzpos + 1.f * aod::v0data::pzneg)))));

//______________________________________________________
// DYNAMIC COLUMNS
// Account for rigidity in case of hypertriton
DECLARE_SOA_DYNAMIC_COLUMN(PtHypertriton, ptHypertriton, //! V0 pT
                           [](float pxpos, float pypos, float pxneg, float pyneg) -> float { return RecoDecay::sqrtSumOfSquares(2.0f * pxpos + pxneg, 2.0f * pypos + pyneg); });
DECLARE_SOA_DYNAMIC_COLUMN(PtAntiHypertriton, ptAntiHypertriton, //! V0 pT
                           [](float pxpos, float pypos, float pxneg, float pyneg) -> float { return RecoDecay::sqrtSumOfSquares(pxpos + 2.0f * pxneg, pypos + 2.0f * pyneg); });

// Length quantities
DECLARE_SOA_DYNAMIC_COLUMN(V0Radius, v0radius, //! V0 decay radius (2D, centered at zero)
                           [](float x, float y) -> float { return RecoDecay::sqrtSumOfSquares(x, y); });

// Distance Over To Mom
DECLARE_SOA_DYNAMIC_COLUMN(DistOverTotMom, distovertotmom, //! PV to V0decay distance over total momentum
                           [](float X, float Y, float Z, float Px, float Py, float Pz, float pvX, float pvY, float pvZ) {
                             float P = RecoDecay::sqrtSumOfSquares(Px, Py, Pz);
                             return std::sqrt(std::pow(X - pvX, 2) + std::pow(Y - pvY, 2) + std::pow(Z - pvZ, 2)) / (P + 1E-10);
                           });

// Armenteros-Podolanski variables
DECLARE_SOA_DYNAMIC_COLUMN(Alpha, alpha, //! Armenteros Alpha
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) {
                             float momTot = RecoDecay::p(pxpos + pxneg, pypos + pyneg, pzpos + pzneg);
                             float lQlNeg = RecoDecay::dotProd(std::array{pxneg, pyneg, pzneg}, std::array{pxpos + pxneg, pypos + pyneg, pzpos + pzneg}) / momTot;
                             float lQlPos = RecoDecay::dotProd(std::array{pxpos, pypos, pzpos}, std::array{pxpos + pxneg, pypos + pyneg, pzpos + pzneg}) / momTot;
                             return (lQlPos - lQlNeg) / (lQlPos + lQlNeg); // alphav0
                           });

DECLARE_SOA_DYNAMIC_COLUMN(QtArm, qtarm, //! Armenteros Qt
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) {
                             float momTot = RecoDecay::p2(pxpos + pxneg, pypos + pyneg, pzpos + pzneg);
                             float dp = RecoDecay::dotProd(std::array{pxneg, pyneg, pzneg}, std::array{pxpos + pxneg, pypos + pyneg, pzpos + pzneg});
                             return std::sqrt(RecoDecay::p2(pxneg, pyneg, pzneg) - dp * dp / momTot); // qtarm
                           });

// Psi pair angle: angle between the plane defined by the electron and positron momenta and the xy plane
DECLARE_SOA_DYNAMIC_COLUMN(PsiPair, psipair, //! psi pair angle
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) {
                             auto clipToPM1 = [](float x) { return x < -1.f ? -1.f : (x > 1.f ? 1.f : x); };
                             float ptot2 = RecoDecay::p2(pxpos, pypos, pzpos) * RecoDecay::p2(pxneg, pyneg, pzneg);
                             float argcos = RecoDecay::dotProd(std::array{pxpos, pypos, pzpos}, std::array{pxneg, pyneg, pzneg}) / std::sqrt(ptot2);
                             float thetaPos = std::atan2(RecoDecay::sqrtSumOfSquares(pxpos, pypos), pzpos);
                             float thetaNeg = std::atan2(RecoDecay::sqrtSumOfSquares(pxneg, pyneg), pzneg);
                             float argsin = (thetaNeg - thetaPos) / std::acos(clipToPM1(argcos));
                             return std::asin(clipToPM1(argsin));
                           });

// calculate the fraction of the pos/neg momentum of the V0 momentum
DECLARE_SOA_DYNAMIC_COLUMN(PFracPos, pfracpos,
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) {
                             float ppos = RecoDecay::sqrtSumOfSquares(pxpos, pypos, pzpos);
                             float PV0 = RecoDecay::sqrtSumOfSquares(pxpos + pxneg, pypos + pyneg, pzpos + pzneg);
                             return (ppos / PV0);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(PFracNeg, pfracneg,
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) {
                             float pneg = RecoDecay::sqrtSumOfSquares(pxneg, pyneg, pzneg);
                             float PV0 = RecoDecay::sqrtSumOfSquares(pxpos + pxneg, pypos + pyneg, pzpos + pzneg);
                             return (pneg / PV0);
                           });

// Calculated on the fly with mass assumption + dynamic tables
DECLARE_SOA_DYNAMIC_COLUMN(MLambda, mLambda, //! mass under lambda hypothesis
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float { return RecoDecay::m(std::array{std::array{pxpos, pypos, pzpos}, std::array{pxneg, pyneg, pzneg}}, std::array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged}); });
DECLARE_SOA_DYNAMIC_COLUMN(MAntiLambda, mAntiLambda, //! mass under antilambda hypothesis
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float { return RecoDecay::m(std::array{std::array{pxpos, pypos, pzpos}, std::array{pxneg, pyneg, pzneg}}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton}); });
DECLARE_SOA_DYNAMIC_COLUMN(MK0Short, mK0Short, //! mass under K0short hypothesis
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float { return RecoDecay::m(std::array{std::array{pxpos, pypos, pzpos}, std::array{pxneg, pyneg, pzneg}}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassPionCharged}); });
DECLARE_SOA_DYNAMIC_COLUMN(MGamma, mGamma, //! mass under gamma hypothesis
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float { return RecoDecay::m(std::array{std::array{pxpos, pypos, pzpos}, std::array{pxneg, pyneg, pzneg}}, std::array{o2::constants::physics::MassElectron, o2::constants::physics::MassElectron}); });
// Account for rigidity in case of hypertriton
DECLARE_SOA_DYNAMIC_COLUMN(MHypertriton, mHypertriton, //! mass under hypertriton  hypothesis
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float { return RecoDecay::m(std::array{std::array{2.0f * pxpos, 2.0f * pypos, 2.0f * pzpos}, std::array{pxneg, pyneg, pzneg}}, std::array{o2::constants::physics::MassHelium3, o2::constants::physics::MassPionCharged}); });
DECLARE_SOA_DYNAMIC_COLUMN(MAntiHypertriton, mAntiHypertriton, //! mass under antihypertriton hypothesis
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float { return RecoDecay::m(std::array{std::array{pxpos, pypos, pzpos}, std::array{2.0f * pxneg, 2.0f * pyneg, 2.0f * pzneg}}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassHelium3}); });
DECLARE_SOA_DYNAMIC_COLUMN(M, m, //! mass under a certain hypothesis (0:K0, 1:L, 2:Lbar, 3:gamma, 4:hyp, 5:ahyp)
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg, int value) -> float {
                             if (value == 0)
                               return RecoDecay::m(std::array{std::array{pxpos, pypos, pzpos}, std::array{pxneg, pyneg, pzneg}}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassPionCharged});
                             if (value == 1)
                               return RecoDecay::m(std::array{std::array{pxpos, pypos, pzpos}, std::array{pxneg, pyneg, pzneg}}, std::array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged});
                             if (value == 2)
                               return RecoDecay::m(std::array{std::array{pxpos, pypos, pzpos}, std::array{pxneg, pyneg, pzneg}}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton});
                             if (value == 3)
                               return RecoDecay::m(std::array{std::array{pxpos, pypos, pzpos}, std::array{pxneg, pyneg, pzneg}}, std::array{o2::constants::physics::MassElectron, o2::constants::physics::MassElectron});
                             if (value == 4)
                               return RecoDecay::m(std::array{std::array{2.0f * pxpos, 2.0f * pypos, 2.0f * pzpos}, std::array{pxneg, pyneg, pzneg}}, std::array{o2::constants::physics::MassHelium3, o2::constants::physics::MassPionCharged});
                             if (value == 5)
                               return RecoDecay::m(std::array{std::array{pxpos, pypos, pzpos}, std::array{2.0f * pxneg, 2.0f * pyneg, 2.0f * pzneg}}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassHelium3});
                             return 0.0f;
                           });

DECLARE_SOA_DYNAMIC_COLUMN(YK0Short, yK0Short, //! V0 y with K0short hypothesis
                           [](float Px, float Py, float Pz) -> float { return RecoDecay::y(std::array{Px, Py, Pz}, o2::constants::physics::MassKaonNeutral); });
DECLARE_SOA_DYNAMIC_COLUMN(YLambda, yLambda, //! V0 y with lambda or antilambda hypothesis
                           [](float Px, float Py, float Pz) -> float { return RecoDecay::y(std::array{Px, Py, Pz}, o2::constants::physics::MassLambda); });
DECLARE_SOA_DYNAMIC_COLUMN(YHypertriton, yHypertriton, //! V0 y with hypertriton hypothesis
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float { return RecoDecay::y(std::array{2.0f * pxpos + pxneg, 2.0f * pypos + pyneg, 2.0f * pzpos + pzneg}, o2::constants::physics::MassHyperTriton); });
DECLARE_SOA_DYNAMIC_COLUMN(YAntiHypertriton, yAntiHypertriton, //! V0 y with antihypertriton hypothesis
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float { return RecoDecay::y(std::array{pxpos + 2.0f * pxneg, pypos + 2.0f * pyneg, pzpos + 2.0f * pzneg}, o2::constants::physics::MassHyperTriton); });
DECLARE_SOA_DYNAMIC_COLUMN(Rapidity, rapidity, //! rapidity (0:K0, 1:L, 2:Lbar)
                           [](float Px, float Py, float Pz, int value) -> float {
                             if (value == 0)
                               return RecoDecay::y(std::array{Px, Py, Pz}, o2::constants::physics::MassKaonNeutral);
                             if (value == 1 || value == 2)
                               return RecoDecay::y(std::array{Px, Py, Pz}, o2::constants::physics::MassLambda);
                             return 0.0f;
                           });

DECLARE_SOA_DYNAMIC_COLUMN(NegativePt, negativept, //! negative daughter pT
                           [](float pxneg, float pyneg) -> float { return RecoDecay::sqrtSumOfSquares(pxneg, pyneg); });
DECLARE_SOA_DYNAMIC_COLUMN(PositivePt, positivept, //! positive daughter pT
                           [](float pxpos, float pypos) -> float { return RecoDecay::sqrtSumOfSquares(pxpos, pypos); });
DECLARE_SOA_DYNAMIC_COLUMN(NegativeEta, negativeeta, //! negative daughter eta
                           [](float PxNeg, float PyNeg, float PzNeg) -> float { return RecoDecay::eta(std::array{PxNeg, PyNeg, PzNeg}); });
DECLARE_SOA_DYNAMIC_COLUMN(NegativePhi, negativephi, //! negative daughter phi
                           [](float PxNeg, float PyNeg) -> float { return RecoDecay::phi(PxNeg, PyNeg); });
DECLARE_SOA_DYNAMIC_COLUMN(PositiveEta, positiveeta, //! positive daughter eta
                           [](float PxPos, float PyPos, float PzPos) -> float { return RecoDecay::eta(std::array{PxPos, PyPos, PzPos}); });
DECLARE_SOA_DYNAMIC_COLUMN(PositivePhi, positivephi, //! positive daughter phi
                           [](float PxPos, float PyPos) -> float { return RecoDecay::phi(PxPos, PyPos); });

DECLARE_SOA_DYNAMIC_COLUMN(IsStandardV0, isStandardV0, //! is standard V0
                           [](uint8_t V0Type) -> bool { return V0Type == 1; });
DECLARE_SOA_DYNAMIC_COLUMN(IsPhotonTPConly, isPhotonTPConly, //! is tpc-only photon V0
                           [](uint8_t V0Type) -> bool { return V0Type & (1 << 1); });
DECLARE_SOA_DYNAMIC_COLUMN(IsCollinear, isCollinear, //! is collinear V0
                           [](uint8_t V0Type) -> bool { return V0Type & (1 << 2); });

DECLARE_SOA_DYNAMIC_COLUMN(RapidityMC, rapidityMC, //! rapidity (0:K0, 1:L, 2:Lbar)
                           [](float PxMC, float PyMC, float PzMC, int value) -> float {
                             if (value == 0)
                               return RecoDecay::y(std::array{PxMC, PyMC, PzMC}, o2::constants::physics::MassKaonNeutral);
                             if (value == 1 || value == 2)
                               return RecoDecay::y(std::array{PxMC, PyMC, PzMC}, o2::constants::physics::MassLambda);
                             return 0.0f;
                           });

DECLARE_SOA_DYNAMIC_COLUMN(NegativePtMC, negativeptMC, //! negative daughter pT
                           [](float pxnegMC, float pynegMC) -> float { return RecoDecay::sqrtSumOfSquares(pxnegMC, pynegMC); });
DECLARE_SOA_DYNAMIC_COLUMN(PositivePtMC, positiveptMC, //! positive daughter pT
                           [](float pxposMC, float pyposMC) -> float { return RecoDecay::sqrtSumOfSquares(pxposMC, pyposMC); });
DECLARE_SOA_DYNAMIC_COLUMN(PtMC, ptMC, //! V0 pT
                           [](float pxMC, float pyMC) -> float { return RecoDecay::sqrtSumOfSquares(pxMC, pyMC); });
} // namespace v0data

DECLARE_SOA_TABLE(V0Indices, "AOD", "V0INDEX", //! index table when using AO2Ds
                  o2::soa::Index<>, v0data::PosTrackId, v0data::NegTrackId, v0data::CollisionId, v0data::V0Id, o2::soa::Marker<1>);

DECLARE_SOA_TABLE(V0CollRefs, "AOD", "V0COLLREF", //! optional table to refer back to a collision
                  o2::soa::Index<>, v0data::StraCollisionId);

DECLARE_SOA_TABLE_STAGED(V0Extras, "V0EXTRA", //! optional table to refer to custom track extras
                         o2::soa::Index<>, v0data::PosTrackExtraId, v0data::NegTrackExtraId);

DECLARE_SOA_TABLE(V0TrackXs, "AOD", "V0TRACKX", //! track X positions at minima when using AO2Ds
                  v0data::PosX, v0data::NegX, o2::soa::Marker<1>);

DECLARE_SOA_TABLE_STAGED(V0CoresBase, "V0CORE", //! core information about decay, viable with AO2Ds or derived
                         o2::soa::Index<>,
                         v0data::X, v0data::Y, v0data::Z,
                         v0data::PxPos, v0data::PyPos, v0data::PzPos,
                         v0data::PxNeg, v0data::PyNeg, v0data::PzNeg,
                         v0data::DCAV0Daughters, v0data::DCAPosToPV, v0data::DCANegToPV,
                         v0data::V0CosPA, v0data::DCAV0ToPV, v0data::V0Type,

                         // Dynamic columns
                         v0data::PtHypertriton<v0data::PxPos, v0data::PyPos, v0data::PxNeg, v0data::PyNeg>,
                         v0data::PtAntiHypertriton<v0data::PxPos, v0data::PyPos, v0data::PxNeg, v0data::PyNeg>,
                         v0data::V0Radius<v0data::X, v0data::Y>,
                         v0data::DistOverTotMom<v0data::X, v0data::Y, v0data::Z, v0data::Px, v0data::Py, v0data::Pz>,
                         v0data::Alpha<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                         v0data::QtArm<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                         v0data::PsiPair<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                         v0data::PFracPos<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                         v0data::PFracNeg<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>, // 24

                         // Invariant masses
                         v0data::MLambda<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                         v0data::MAntiLambda<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                         v0data::MK0Short<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                         v0data::MGamma<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                         v0data::MHypertriton<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                         v0data::MAntiHypertriton<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                         v0data::M<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,

                         // Longitudinal
                         v0data::YK0Short<v0data::Px, v0data::Py, v0data::Pz>,
                         v0data::YLambda<v0data::Px, v0data::Py, v0data::Pz>,
                         v0data::YHypertriton<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                         v0data::YAntiHypertriton<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                         v0data::Rapidity<v0data::Px, v0data::Py, v0data::Pz>,
                         v0data::NegativePt<v0data::PxNeg, v0data::PyNeg>,
                         v0data::PositivePt<v0data::PxPos, v0data::PyPos>,
                         v0data::NegativeEta<v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                         v0data::NegativePhi<v0data::PxNeg, v0data::PyNeg>,
                         v0data::PositiveEta<v0data::PxPos, v0data::PyPos, v0data::PzPos>,
                         v0data::PositivePhi<v0data::PxPos, v0data::PyPos>,
                         v0data::IsStandardV0<v0data::V0Type>,
                         v0data::IsPhotonTPConly<v0data::V0Type>);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(V0Cores, V0CoresBase, "V0COREEXT",                                                    //!
                                v0data::Px, v0data::Py, v0data::Pz, v0data::Pt, v0data::P, v0data::Phi, v0data::Eta); // the table name has here to be the one with EXT which is not nice and under study

// // extended table with expression columns that can be used as arguments of dynamic columns
// DECLARE_SOA_EXTENDED_TABLE_USER(StoredV0Cores, StoredV0CoresBase, "V0COREEXT",                                                            //!
//                                 v0data::Px, v0data::Py, v0data::Pz, v0data::Pt, v0data::P, v0data::Phi, v0data::Eta, o2::soa::Marker<2>); // the table name has here to be the one with EXT which is not nice and under study

DECLARE_SOA_TABLE(V0TraPosAtDCAs, "AOD", "V0TRAPOSATDCAs", //! positions of tracks at their DCA for debug
                  v0data::XPosAtDCA, v0data::YPosAtDCA, v0data::ZPosAtDCA,
                  v0data::XNegAtDCA, v0data::YNegAtDCA, v0data::ZNegAtDCA);
DECLARE_SOA_TABLE(V0TraPosAtIUs, "AOD", "V0TRAPOSATIUs", //! positions of tracks at their IU for debug
                  v0data::XPosAtIU, v0data::YPosAtIU, v0data::ZPosAtIU,
                  v0data::XNegAtIU, v0data::YNegAtIU, v0data::ZNegAtIU);

DECLARE_SOA_TABLE(V0Ivanovs, "AOD", "V0Ivanovs", //! bitmaps for Marian
                  v0data::IvanovMap);

DECLARE_SOA_TABLE_FULL(V0Covs, "V0Covs", "AOD", "V0COVS", //! V0 covariance matrices
                       v0data::PositionCovMat, v0data::MomentumCovMat, o2::soa::Marker<1>);

DECLARE_SOA_TABLE_FULL(V0DauCovs, "V0DauCovs", "AOD", "V0DAUCOVS", //! V0 covariance matrices of the dauther tracks
                       v0data::CovMatPosDau, v0data::CovMatNegDau, o2::soa::Marker<1>);
DECLARE_SOA_TABLE_FULL(V0DauCovIUs, "V0DauCovIUs", "AOD", "V0DAUCOVIUS", //! V0 covariance matrices of the dauther tracks
                       v0data::CovMatPosDauIU, v0data::CovMatNegDauIU, o2::soa::Marker<1>);

DECLARE_SOA_TABLE(V0fCIndices, "AOD", "V0FCINDEX", //! index table when using AO2Ds
                  o2::soa::Index<>, v0data::PosTrackId, v0data::NegTrackId, v0data::CollisionId, v0data::V0Id, o2::soa::Marker<2>);

DECLARE_SOA_TABLE(V0fCTrackXs, "AOD", "V0FCTRACKX", //! track X positions at minima when using AO2Ds
                  v0data::PosX, v0data::NegX, o2::soa::Marker<2>);

DECLARE_SOA_TABLE_FULL(StoredV0fCCores, "V0fCCores", "AOD", "V0FCCORE", //! core information about decay, exclusive to V0s used in cascades but not passing in V0 selections - multiple viable getters skipped since use is protected to cascades only
                       v0data::X, v0data::Y, v0data::Z,
                       v0data::PxPos, v0data::PyPos, v0data::PzPos,
                       v0data::PxNeg, v0data::PyNeg, v0data::PzNeg,
                       v0data::DCAV0Daughters, v0data::DCAPosToPV, v0data::DCANegToPV,
                       v0data::V0CosPA, v0data::DCAV0ToPV, v0data::V0Type,

                       // Dynamic columns
                       v0data::PtHypertriton<v0data::PxPos, v0data::PyPos, v0data::PxNeg, v0data::PyNeg>,
                       v0data::PtAntiHypertriton<v0data::PxPos, v0data::PyPos, v0data::PxNeg, v0data::PyNeg>,
                       v0data::V0Radius<v0data::X, v0data::Y>,
                       v0data::DistOverTotMom<v0data::X, v0data::Y, v0data::Z, v0data::Px, v0data::Py, v0data::Pz>,
                       v0data::Alpha<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                       v0data::QtArm<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                       v0data::PsiPair<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                       v0data::PFracPos<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                       v0data::PFracNeg<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>, // 24

                       // Invariant masses
                       v0data::MLambda<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                       v0data::MAntiLambda<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                       v0data::MK0Short<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                       v0data::MGamma<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                       v0data::MHypertriton<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                       v0data::MAntiHypertriton<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                       v0data::M<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,

                       // Longitudinal
                       v0data::YK0Short<v0data::Px, v0data::Py, v0data::Pz>,
                       v0data::YLambda<v0data::Px, v0data::Py, v0data::Pz>,
                       v0data::YHypertriton<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                       v0data::YAntiHypertriton<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                       v0data::Rapidity<v0data::Px, v0data::Py, v0data::Pz>,
                       v0data::NegativePt<v0data::PxNeg, v0data::PyNeg>,
                       v0data::PositivePt<v0data::PxPos, v0data::PyPos>,
                       v0data::NegativeEta<v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                       v0data::NegativePhi<v0data::PxNeg, v0data::PyNeg>,
                       v0data::PositiveEta<v0data::PxPos, v0data::PyPos, v0data::PzPos>,
                       v0data::PositivePhi<v0data::PxPos, v0data::PyPos>,
                       v0data::IsStandardV0<v0data::V0Type>,
                       v0data::IsPhotonTPConly<v0data::V0Type>,
                       o2::soa::Marker<2>);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(V0fCCores, StoredV0fCCores, "V0FCCOREEXT",                                            //!
                                v0data::Px, v0data::Py, v0data::Pz, v0data::Pt, v0data::P, v0data::Phi, v0data::Eta); // the table name has here to be the one with EXT which is not nice and under study

DECLARE_SOA_TABLE_FULL(V0fCCovs, "V0fCCovs", "AOD", "V0FCCOVS", //! V0 covariance matrices
                       v0data::PositionCovMat, v0data::MomentumCovMat, o2::soa::Marker<2>);

DECLARE_SOA_TABLE_STAGED(V0MCCores_000, "V0MCCORE", //! MC properties of the V0 for posterior analysis
                         v0data::PDGCode, v0data::PDGCodeMother,
                         v0data::PDGCodePositive, v0data::PDGCodeNegative,
                         v0data::IsPhysicalPrimary, v0data::XMC, v0data::YMC, v0data::ZMC,
                         v0data::PxPosMC, v0data::PyPosMC, v0data::PzPosMC,
                         v0data::PxNegMC, v0data::PyNegMC, v0data::PzNegMC);

DECLARE_SOA_TABLE_STAGED_VERSIONED(V0MCCores_001, "V0MCCORE", 1, //! debug information
                                   v0data::ParticleIdMC,         //! MC properties of the V0 for posterior analysis
                                   v0data::PDGCode, v0data::PDGCodeMother,
                                   v0data::PDGCodePositive, v0data::PDGCodeNegative,
                                   v0data::IsPhysicalPrimary, v0data::XMC, v0data::YMC, v0data::ZMC,
                                   v0data::PxPosMC, v0data::PyPosMC, v0data::PzPosMC,
                                   v0data::PxNegMC, v0data::PyNegMC, v0data::PzNegMC);

DECLARE_SOA_TABLE_STAGED_VERSIONED(V0MCCores_002, "V0MCCORE", 2, //! debug information
                                   v0data::ParticleIdMC,         //! MC properties of the V0 for posterior analysis
                                   v0data::PDGCode, v0data::PDGCodeMother,
                                   v0data::PDGCodePositive, v0data::PDGCodeNegative,
                                   v0data::IsPhysicalPrimary, v0data::XMC, v0data::YMC, v0data::ZMC,
                                   v0data::PxPosMC, v0data::PyPosMC, v0data::PzPosMC,
                                   v0data::PxNegMC, v0data::PyNegMC, v0data::PzNegMC,
                                   v0data::PxMC, v0data::PyMC, v0data::PzMC,
                                   v0data::RapidityMC<v0data::PxMC, v0data::PyMC, v0data::PzMC>,
                                   v0data::NegativePtMC<v0data::PxNegMC, v0data::PyNegMC>,
                                   v0data::PositivePtMC<v0data::PxPosMC, v0data::PyPosMC>,
                                   v0data::PtMC<v0data::PxMC, v0data::PyMC>);

// DECLARE_SOA_TABLE(StoredV0MCCores_000, "AOD", "V0MCCORE", //! MC properties of the V0 for posterior analysis
//                   v0data::PDGCode, v0data::PDGCodeMother,
//                   v0data::PDGCodePositive, v0data::PDGCodeNegative,
//                   v0data::IsPhysicalPrimary, v0data::XMC, v0data::YMC, v0data::ZMC,
//                   v0data::PxPosMC, v0data::PyPosMC, v0data::PzPosMC,
//                   v0data::PxNegMC, v0data::PyNegMC, v0data::PzNegMC,
//                   o2::soa::Marker<1>);

// DECLARE_SOA_TABLE_VERSIONED(StoredV0MCCores_001, "AOD", "V0MCCORE", 1, //! debug information
//                             v0data::ParticleIdMC,                      //! MC properties of the V0 for posterior analysis
//                             v0data::PDGCode, v0data::PDGCodeMother,
//                             v0data::PDGCodePositive, v0data::PDGCodeNegative,
//                             v0data::IsPhysicalPrimary, v0data::XMC, v0data::YMC, v0data::ZMC,
//                             v0data::PxPosMC, v0data::PyPosMC, v0data::PzPosMC,
//                             v0data::PxNegMC, v0data::PyNegMC, v0data::PzNegMC,
//                             o2::soa::Marker<1>);

// DECLARE_SOA_TABLE_VERSIONED(StoredV0MCCores_002, "AOD", "V0MCCORE", 2, //! debug information
//                             v0data::ParticleIdMC,                      //! MC properties of the V0 for posterior analysis
//                             v0data::PDGCode, v0data::PDGCodeMother,
//                             v0data::PDGCodePositive, v0data::PDGCodeNegative,
//                             v0data::IsPhysicalPrimary, v0data::XMC, v0data::YMC, v0data::ZMC,
//                             v0data::PxPosMC, v0data::PyPosMC, v0data::PzPosMC,
//                             v0data::PxNegMC, v0data::PyNegMC, v0data::PzNegMC,
//                             v0data::PxMC, v0data::PyMC, v0data::PzMC,
//                             o2::soa::Marker<1>);

DECLARE_SOA_TABLE(V0MCCollRefs, "AOD", "V0MCCOLLREF", //! refers MC candidate back to proper MC Collision
                  o2::soa::Index<>, v0data::StraMCCollisionId, o2::soa::Marker<2>);

DECLARE_SOA_TABLE(GeK0Short, "AOD", "GeK0Short", v0data::GeneratedK0Short);
DECLARE_SOA_TABLE(GeLambda, "AOD", "GeLambda", v0data::GeneratedLambda);
DECLARE_SOA_TABLE(GeAntiLambda, "AOD", "GeAntiLambda", v0data::GeneratedAntiLambda);

DECLARE_SOA_TABLE_STAGED(V0MCMothers, "V0MCMOTHER", //! optional table for MC mothers
                         o2::soa::Index<>, v0data::MotherMCPartId);

using V0MCCores = V0MCCores_002;
using StoredV0MCCores = StoredV0MCCores_002;

using V0Index = V0Indices::iterator;
using V0Core = V0Cores::iterator;
using V0TrackX = V0TrackXs::iterator;
using V0Datas = soa::Join<V0Indices, V0TrackXs, V0Cores>;
using V0Data = V0Datas::iterator;
using V0fCDatas = soa::Join<V0fCIndices, V0fCTrackXs, V0fCCores>;
using V0fCData = V0fCDatas::iterator;
using V0MCDatas = soa::Join<V0MCCores, V0MCMothers>;
using V0MCData = V0MCDatas::iterator;
using V0MCCore = V0MCCores::iterator;

// definitions of indices for interlink tables
namespace v0data
{
DECLARE_SOA_INDEX_COLUMN(V0Data, v0Data);                         //! Index to V0Data entry
DECLARE_SOA_INDEX_COLUMN(V0fCData, v0fCData);                     //! Index to V0Data entry
DECLARE_SOA_INDEX_COLUMN_FULL(V0MC, v0MC, int, V0MCCores, "_MC"); //!
DECLARE_SOA_INDEX_COLUMN(V0MCCore, v0MCCore);
} // namespace v0data

DECLARE_SOA_TABLE(V0DataLink, "AOD", "V0DATALINK", //! Joinable table with V0s which links to V0Data which is not produced for all entries
                  o2::soa::Index<>, v0data::V0DataId, v0data::V0fCDataId);
DECLARE_SOA_TABLE(V0MCRefs, "AOD", "V0MCREF", //! index table when using AO2Ds
                  o2::soa::Index<>, v0data::V0MCId);
DECLARE_SOA_TABLE(V0CoreMCLabels, "AOD", "V0COREMCLABEL", //! optional table to refer to V0MCCores if not joinable
                  o2::soa::Index<>, v0data::V0MCCoreId);

using V0sLinked = soa::Join<V0s, V0DataLink>;
using V0Linked = V0sLinked::iterator;

namespace v0data
{
DECLARE_SOA_COLUMN(IsFound, isFound, bool); //! is this FindableV0 actually in the V0s table?
}

// Major bypass for simultaneous found vs findable study
DECLARE_SOA_TABLE(FindableV0s, "AOD", "FindableV0", //! Will store findable
                  o2::soa::Index<>, v0::CollisionId,
                  v0::PosTrackId, v0::NegTrackId,
                  v0::V0Type,
                  v0::IsStandardV0<v0::V0Type>,
                  v0::IsPhotonV0<v0::V0Type>,
                  v0::IsCollinearV0<v0::V0Type>,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(V0FoundTags, "AOD", "V0FoundTag", //! found or not?
                  v0data::IsFound);

using FindableV0sLinked = soa::Join<FindableV0s, V0DataLink>;
using FindableV0Linked = FindableV0sLinked::iterator;

// helper for building
namespace v0tag
{
// Global bool
DECLARE_SOA_COLUMN(IsInteresting, isInteresting, bool); //! will this be built or not?

// MC association bools
DECLARE_SOA_COLUMN(IsTrueGamma, isTrueGamma, bool);                     //! PDG checked correctly in MC
DECLARE_SOA_COLUMN(IsTrueK0Short, isTrueK0Short, bool);                 //! PDG checked correctly in MC
DECLARE_SOA_COLUMN(IsTrueLambda, isTrueLambda, bool);                   //! PDG checked correctly in MC
DECLARE_SOA_COLUMN(IsTrueAntiLambda, isTrueAntiLambda, bool);           //! PDG checked correctly in MC
DECLARE_SOA_COLUMN(IsTrueHypertriton, isTrueHypertriton, bool);         //! PDG checked correctly in MC
DECLARE_SOA_COLUMN(IsTrueAntiHypertriton, isTrueAntiHypertriton, bool); //! PDG checked correctly in MC
DECLARE_SOA_COLUMN(IsPhysicalPrimary, isPhysicalPrimary, bool);         //! physical primary

// dE/dx compatibility bools
DECLARE_SOA_COLUMN(IsdEdxGamma, isdEdxGamma, bool);                     //! compatible with dE/dx hypotheses
DECLARE_SOA_COLUMN(IsdEdxK0Short, isdEdxK0Short, bool);                 //! compatible with dE/dx hypotheses
DECLARE_SOA_COLUMN(IsdEdxLambda, isdEdxLambda, bool);                   //! compatible with dE/dx hypotheses
DECLARE_SOA_COLUMN(IsdEdxAntiLambda, isdEdxAntiLambda, bool);           //! compatible with dE/dx hypotheses
DECLARE_SOA_COLUMN(IsdEdxHypertriton, isdEdxHypertriton, bool);         //! compatible with dE/dx hypotheses
DECLARE_SOA_COLUMN(IsdEdxAntiHypertriton, isdEdxAntiHypertriton, bool); //! compatible with dE/dx hypotheses

// used in cascades (potentially useful in general, make available as tags)
DECLARE_SOA_COLUMN(IsFromCascade, isFromCascade, bool);               //! compatible with dE/dx hypotheses
DECLARE_SOA_COLUMN(IsFromTrackedCascade, isFromTrackedCascade, bool); //! compatible with dE/dx hypotheses
} // namespace v0tag
DECLARE_SOA_TABLE(V0Tags, "AOD", "V0TAGS",
                  v0tag::IsInteresting,
                  v0tag::IsTrueGamma,
                  v0tag::IsTrueK0Short,
                  v0tag::IsTrueLambda,
                  v0tag::IsTrueAntiLambda,
                  v0tag::IsTrueHypertriton,
                  v0tag::IsTrueAntiHypertriton,
                  v0tag::IsPhysicalPrimary,
                  v0tag::IsdEdxGamma,
                  v0tag::IsdEdxK0Short,
                  v0tag::IsdEdxLambda,
                  v0tag::IsdEdxAntiLambda,
                  v0tag::IsdEdxHypertriton,
                  v0tag::IsdEdxAntiHypertriton,
                  v0tag::IsFromCascade,
                  v0tag::IsFromTrackedCascade);

namespace kfcascdata
{
// declare in different namespace to 'overload' operator
DECLARE_SOA_COLUMN(MLambda, mLambda, float); //!
} // namespace kfcascdata

namespace cascdata
{
//______________________________________________________
// REGULAR COLUMNS FOR CASCINDICES
DECLARE_SOA_INDEX_COLUMN(V0, v0);                                                   //!
DECLARE_SOA_INDEX_COLUMN(Cascade, cascade);                                         //!
DECLARE_SOA_INDEX_COLUMN_FULL(PosTrack, posTrack, int, Tracks, "_Pos");             //!
DECLARE_SOA_INDEX_COLUMN_FULL(NegTrack, negTrack, int, Tracks, "_Neg");             //!
DECLARE_SOA_INDEX_COLUMN_FULL(Bachelor, bachelor, int, Tracks, "_Bach");            //!
DECLARE_SOA_INDEX_COLUMN_FULL(StrangeTrack, strangeTrack, int, Tracks, "_Strange"); //!
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                                     //!
// FOR DERIVED
DECLARE_SOA_INDEX_COLUMN_FULL(PosTrackExtra, posTrackExtra, int, DauTrackExtras, "_PosExtra");             //!
DECLARE_SOA_INDEX_COLUMN_FULL(NegTrackExtra, negTrackExtra, int, DauTrackExtras, "_NegExtra");             //!
DECLARE_SOA_INDEX_COLUMN_FULL(BachTrackExtra, bachTrackExtra, int, DauTrackExtras, "_BachExtra");          //!
DECLARE_SOA_INDEX_COLUMN_FULL(StrangeTrackExtra, strangeTrackExtra, int, DauTrackExtras, "_StrangeExtra"); //!
DECLARE_SOA_INDEX_COLUMN(StraCollision, straCollision);                                                    //!
DECLARE_SOA_INDEX_COLUMN(StraMCCollision, straMCCollision);                                                //!
DECLARE_SOA_INDEX_COLUMN(MotherMCPart, motherMCPart);                                                      //!

//______________________________________________________
// REGULAR COLUMNS FOR CASCCORES
// General cascade properties: position, momentum
DECLARE_SOA_COLUMN(Sign, sign, int);         //!
DECLARE_SOA_COLUMN(MXi, mXi, float);         //!
DECLARE_SOA_COLUMN(MOmega, mOmega, float);   //!
DECLARE_SOA_COLUMN(PxPos, pxpos, float);     //!
DECLARE_SOA_COLUMN(PyPos, pypos, float);     //!
DECLARE_SOA_COLUMN(PzPos, pzpos, float);     //!
DECLARE_SOA_COLUMN(PxNeg, pxneg, float);     //!
DECLARE_SOA_COLUMN(PyNeg, pyneg, float);     //!
DECLARE_SOA_COLUMN(PzNeg, pzneg, float);     //!
DECLARE_SOA_COLUMN(PxBach, pxbach, float);   //!
DECLARE_SOA_COLUMN(PyBach, pybach, float);   //!
DECLARE_SOA_COLUMN(PzBach, pzbach, float);   //!
DECLARE_SOA_COLUMN(Px, px, float);           //! cascade momentum X
DECLARE_SOA_COLUMN(Py, py, float);           //! cascade momentum Y
DECLARE_SOA_COLUMN(Pz, pz, float);           //! cascade momentum Z
DECLARE_SOA_COLUMN(X, x, float);             //!
DECLARE_SOA_COLUMN(Y, y, float);             //!
DECLARE_SOA_COLUMN(Z, z, float);             //!
DECLARE_SOA_COLUMN(Xlambda, xlambda, float); //!
DECLARE_SOA_COLUMN(Ylambda, ylambda, float); //!
DECLARE_SOA_COLUMN(Zlambda, zlambda, float); //!

// Saved from finding: DCAs
DECLARE_SOA_COLUMN(DCAV0Daughters, dcaV0daughters, float);     //!
DECLARE_SOA_COLUMN(DCACascDaughters, dcacascdaughters, float); //!
DECLARE_SOA_COLUMN(DCAPosToPV, dcapostopv, float);             //!
DECLARE_SOA_COLUMN(DCANegToPV, dcanegtopv, float);             //!
DECLARE_SOA_COLUMN(DCABachToPV, dcabachtopv, float);           //!
DECLARE_SOA_COLUMN(DCAXYCascToPV, dcaXYCascToPV, float);       //!
DECLARE_SOA_COLUMN(DCAZCascToPV, dcaZCascToPV, float);         //!

// Saved from finding: track position at minima
DECLARE_SOA_COLUMN(PosX, posX, float);   //! positive track X at min
DECLARE_SOA_COLUMN(NegX, negX, float);   //! negative track X at min
DECLARE_SOA_COLUMN(BachX, bachX, float); //! bachelor track X at min

//______________________________________________________
// REGULAR COLUMNS FOR CASCCOVS
// Saved from finding: covariance matrix of parent track (on request)
DECLARE_SOA_DYNAMIC_COLUMN(PositionCovMat, positionCovMat, //! for transparent handling
                           [](const float covMat[21]) -> std::vector<float> {
                            std::vector<float> posCovMat { covMat[0], covMat[1], covMat[2], covMat[3], covMat[4], covMat[5] };
                            return posCovMat; });
DECLARE_SOA_DYNAMIC_COLUMN(MomentumCovMat, momentumCovMat, //! for transparent handling
                           [](const float covMat[21]) -> std::vector<float> {
                            std::vector<float> momCovMat { covMat[9], covMat[13], covMat[14], covMat[18], covMat[19], covMat[20] };
                            return momCovMat; });
DECLARE_SOA_COLUMN(KFTrackCovMat, kfTrackCovMat, float[21]);                 //! covariance matrix elements for KF method (Cascade)
DECLARE_SOA_COLUMN(KFTrackCovMatV0, kfTrackCovMatV0, float[21]);             //! covariance matrix elements for KF method (V0)
DECLARE_SOA_COLUMN(KFTrackCovMatV0DauPos, kfTrackCovMatV0DauPos, float[21]); //! covariance matrix elements for KF method (V0 pos daughter)
DECLARE_SOA_COLUMN(KFTrackCovMatV0DauNeg, kfTrackCovMatV0DauNeg, float[21]); //! covariance matrix elements for KF method (V0 neg daughter)

// for CascCovs / TraCascCovs, meant to provide consistent interface everywhere
DECLARE_SOA_COLUMN(CovMat, covMat, float[21]); //! covariance matrix elements

//______________________________________________________
// REGULAR COLUMNS FOR CASCBBS
// General cascade properties: position, momentum
// Selection to avoid spurious invariant mass correlation
// bachelor-baryon cosine of pointing angle / DCA to PV
DECLARE_SOA_COLUMN(BachBaryonCosPA, bachBaryonCosPA, float);         //! avoid bach-baryon correlated inv mass structure in analysis
DECLARE_SOA_COLUMN(BachBaryonDCAxyToPV, bachBaryonDCAxyToPV, float); //! avoid bach-baryon correlated inv mass structure in analysis

//______________________________________________________
// REGULAR COLUMNS FOR KFCASCCORES
// General cascade properties: position, momentum
// Saved from KF particle fit for specic table
// note: separate chi2 is a consequence of fit -> conversion -> propagation -> fit logic
//       which, in turn, is necessary to do material corrections at the moment
//       this could be improved in the future!
DECLARE_SOA_COLUMN(KFV0Chi2, kfV0Chi2, float);           //!
DECLARE_SOA_COLUMN(KFCascadeChi2, kfCascadeChi2, float); //!
DECLARE_SOA_COLUMN(KFPxV0, kfpxv0, float);               //!
DECLARE_SOA_COLUMN(KFPyV0, kfpyv0, float);               //!
DECLARE_SOA_COLUMN(KFPzV0, kfpzv0, float);               //!
DECLARE_SOA_COLUMN(KFXPos, kfxpos, float);               //!
DECLARE_SOA_COLUMN(KFYPos, kfypos, float);               //!
DECLARE_SOA_COLUMN(KFZPos, kfzpos, float);               //!
DECLARE_SOA_COLUMN(KFXNeg, kfxneg, float);               //!
DECLARE_SOA_COLUMN(KFYNeg, kfyneg, float);               //!
DECLARE_SOA_COLUMN(KFZNeg, kfzneg, float);               //!

//______________________________________________________
// REGULAR COLUMNS FOR TRACASCCORES
// Saved from strangeness tracking
DECLARE_SOA_COLUMN(MatchingChi2, matchingChi2, float); //!
DECLARE_SOA_COLUMN(TopologyChi2, topologyChi2, float); //!
DECLARE_SOA_COLUMN(ItsClsSize, itsCluSize, float);     //!

//______________________________________________________
// REGULAR COLUMNS FOR CASCMCCORES
DECLARE_SOA_COLUMN(PDGCode, pdgCode, int);                      //! cascade PDG Code
DECLARE_SOA_COLUMN(PDGCodeMother, pdgCodeMother, int);          //! cascade mother PDG code (for feeddown)
DECLARE_SOA_COLUMN(PDGCodeV0, pdgCodeV0, int);                  //! cascade PDG Code
DECLARE_SOA_COLUMN(PDGCodePositive, pdgCodePositive, int);      //! V0 positive prong PDG code
DECLARE_SOA_COLUMN(PDGCodeNegative, pdgCodeNegative, int);      //! V0 negative prong PDG code
DECLARE_SOA_COLUMN(PDGCodeBachelor, pdgCodeBachelor, int);      //! cascade bachelor prong PDG code
DECLARE_SOA_COLUMN(IsPhysicalPrimary, isPhysicalPrimary, bool); //! is cascade physical primary
DECLARE_SOA_COLUMN(XMC, xMC, float);                            //! cascade decay position X (cm)
DECLARE_SOA_COLUMN(YMC, yMC, float);                            //! cascade decay position Y (cm)
DECLARE_SOA_COLUMN(ZMC, zMC, float);                            //! cascade decay position Z (cm)
DECLARE_SOA_COLUMN(XlambdaMC, xlambdaMC, float);                //! V0 decay position X (cm)
DECLARE_SOA_COLUMN(YlambdaMC, ylambdaMC, float);                //! V0 decay position Y (cm)
DECLARE_SOA_COLUMN(ZlambdaMC, zlambdaMC, float);                //! V0 decay position Z (cm)
DECLARE_SOA_COLUMN(PxPosMC, pxPosMC, float);                    //! V0 positive daughter px (GeV/c)
DECLARE_SOA_COLUMN(PyPosMC, pyPosMC, float);                    //! V0 positive daughter py (GeV/c)
DECLARE_SOA_COLUMN(PzPosMC, pzPosMC, float);                    //! V0 positive daughter pz (GeV/c)
DECLARE_SOA_COLUMN(PxNegMC, pxNegMC, float);                    //! V0 positive daughter px (GeV/c)
DECLARE_SOA_COLUMN(PyNegMC, pyNegMC, float);                    //! V0 positive daughter py (GeV/c)
DECLARE_SOA_COLUMN(PzNegMC, pzNegMC, float);                    //! V0 positive daughter pz (GeV/c)
DECLARE_SOA_COLUMN(PxBachMC, pxBachMC, float);                  //! cascade bachelor daughter px (GeV/c)
DECLARE_SOA_COLUMN(PyBachMC, pyBachMC, float);                  //! cascade bachelor daughter py (GeV/c)
DECLARE_SOA_COLUMN(PzBachMC, pzBachMC, float);                  //! cascade bachelor daughter pz (GeV/c)
DECLARE_SOA_COLUMN(PxMC, pxMC, float);                          //! cascade px (GeV/c)
DECLARE_SOA_COLUMN(PyMC, pyMC, float);                          //! cascade py (GeV/c)
DECLARE_SOA_COLUMN(PzMC, pzMC, float);                          //! cascade pz (GeV/c)

//______________________________________________________
// generated binned data
DECLARE_SOA_COLUMN(GeneratedXiMinus, generatedXiMinus, std::vector<uint32_t>);       //! XiMinus binned generated data
DECLARE_SOA_COLUMN(GeneratedXiPlus, generatedXiPlus, std::vector<uint32_t>);         //! XiPlus binned generated data
DECLARE_SOA_COLUMN(GeneratedOmegaMinus, generatedOmegaMinus, std::vector<uint32_t>); //! OmegaMinus binned generated data
DECLARE_SOA_COLUMN(GeneratedOmegaPlus, generatedOmegaPlus, std::vector<uint32_t>);   //! OmegaPlus binned generated data

//______________________________________________________
// DERIVED
// Length quantities
DECLARE_SOA_DYNAMIC_COLUMN(V0Radius, v0radius, //!
                           [](float xlambda, float ylambda) -> float { return RecoDecay::sqrtSumOfSquares(xlambda, ylambda); });
DECLARE_SOA_DYNAMIC_COLUMN(CascRadius, cascradius, //!
                           [](float x, float y) -> float { return RecoDecay::sqrtSumOfSquares(x, y); });

// CosPAs
DECLARE_SOA_DYNAMIC_COLUMN(V0CosPA, v0cosPA, //!
                           [](float Xlambda, float Ylambda, float Zlambda, float PxLambda, float PyLambda, float PzLambda, float pvX, float pvY, float pvZ) -> float { return RecoDecay::cpa(std::array{pvX, pvY, pvZ}, std::array{Xlambda, Ylambda, Zlambda}, std::array{PxLambda, PyLambda, PzLambda}); });
// DECLARE_SOA_DYNAMIC_COLUMN(CascCosPA, casccosPA, //!
//                            [](float X, float Y, float Z, float Px, float Py, float Pz, float pvX, float pvY, float pvZ) -> float { return RecoDecay::cpa(std::array{pvX, pvY, pvZ}, std::array{X, Y, Z}, std::array{Px, Py, Pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(CascCosPA, casccosPA, //!
                           [](float X, float Y, float Z, float PxBach, float PxPos, float PxNeg, float PyBach, float PyPos, float PyNeg, float PzBach, float PzPos, float PzNeg, float pvX, float pvY, float pvZ) -> float { return RecoDecay::cpa(std::array{pvX, pvY, pvZ}, std::array{X, Y, Z}, std::array{PxBach + PxPos + PxNeg, PyBach + PyPos + PyNeg, PzBach + PzPos + PzNeg}); });
DECLARE_SOA_DYNAMIC_COLUMN(DCAV0ToPV, dcav0topv, //!
                           [](float X, float Y, float Z, float Px, float Py, float Pz, float pvX, float pvY, float pvZ) -> float { return std::sqrt((std::pow((pvY - Y) * Pz - (pvZ - Z) * Py, 2) + std::pow((pvX - X) * Pz - (pvZ - Z) * Px, 2) + std::pow((pvX - X) * Py - (pvY - Y) * Px, 2)) / (Px * Px + Py * Py + Pz * Pz)); });

// Calculated on the fly with mass assumption + dynamic tables
DECLARE_SOA_DYNAMIC_COLUMN(MLambda, mLambda, //!
                           [](int charge, float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float { return RecoDecay::m(std::array{std::array{pxpos, pypos, pzpos}, std::array{pxneg, pyneg, pzneg}}, charge < 0 ? std::array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged} : std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton}); });
DECLARE_SOA_DYNAMIC_COLUMN(M, m, //! mass under a certain hypothesis (0:K0, 1:L, 2:Lbar, 3:gamma, 4:hyp, 5:ahyp)
                           [](float mXi, float mOmega, int value) -> float {
                             if (value == 0 || value == 1)
                               return mXi;
                             if (value == 2 || value == 3)
                               return mOmega;
                             return 0.0f;
                           });

DECLARE_SOA_DYNAMIC_COLUMN(YXi, yXi, //!
                           [](float Px, float Py, float Pz) -> float { return RecoDecay::y(std::array{Px, Py, Pz}, o2::constants::physics::MassXiMinus); });
DECLARE_SOA_DYNAMIC_COLUMN(YOmega, yOmega, //!
                           [](float Px, float Py, float Pz) -> float { return RecoDecay::y(std::array{Px, Py, Pz}, o2::constants::physics::MassOmegaMinus); });
DECLARE_SOA_DYNAMIC_COLUMN(Rapidity, rapidity, //! rapidity (0, 1: Xi; 2, 3: Omega)
                           [](float Px, float Py, float Pz, int value) -> float {
                             if (value == 0 || value == 1)
                               return RecoDecay::y(std::array{Px, Py, Pz}, o2::constants::physics::MassXiMinus);
                             if (value == 2 || value == 3)
                               return RecoDecay::y(std::array{Px, Py, Pz}, o2::constants::physics::MassOmegaMinus);
                             return 0.0f;
                           });

DECLARE_SOA_DYNAMIC_COLUMN(NegativePt, negativept, //! negative daughter pT
                           [](float pxneg, float pyneg) -> float { return RecoDecay::sqrtSumOfSquares(pxneg, pyneg); });
DECLARE_SOA_DYNAMIC_COLUMN(PositivePt, positivept, //! positive daughter pT
                           [](float pxpos, float pypos) -> float { return RecoDecay::sqrtSumOfSquares(pxpos, pypos); });
DECLARE_SOA_DYNAMIC_COLUMN(BachelorPt, bachelorpt, //! bachelor daughter pT
                           [](float pxpos, float pypos) -> float { return RecoDecay::sqrtSumOfSquares(pxpos, pypos); });
DECLARE_SOA_DYNAMIC_COLUMN(NegativeEta, negativeeta, //! negative daughter eta
                           [](float PxNeg, float PyNeg, float PzNeg) -> float { return RecoDecay::eta(std::array{PxNeg, PyNeg, PzNeg}); });
DECLARE_SOA_DYNAMIC_COLUMN(NegativePhi, negativephi, //! negative daughter phi
                           [](float PxNeg, float PyNeg) -> float { return RecoDecay::phi(PxNeg, PyNeg); });
DECLARE_SOA_DYNAMIC_COLUMN(PositiveEta, positiveeta, //! positive daughter eta
                           [](float PxPos, float PyPos, float PzPos) -> float { return RecoDecay::eta(std::array{PxPos, PyPos, PzPos}); });
DECLARE_SOA_DYNAMIC_COLUMN(PositivePhi, positivephi, //! positive daughter phi
                           [](float PxPos, float PyPos) -> float { return RecoDecay::phi(PxPos, PyPos); });
DECLARE_SOA_DYNAMIC_COLUMN(BachelorEta, bacheloreta, //! bachelor daughter eta
                           [](float PxPos, float PyPos, float PzPos) -> float { return RecoDecay::eta(std::array{PxPos, PyPos, PzPos}); });
DECLARE_SOA_DYNAMIC_COLUMN(BachelorPhi, bachelorphi, //! bachelor daughter phi
                           [](float PxPos, float PyPos) -> float { return RecoDecay::phi(PxPos, PyPos); });

DECLARE_SOA_DYNAMIC_COLUMN(RapidityMC, rapidityMC, //! rapidity (0, 1: Xi; 2, 3: Omega)
                           [](float PxMC, float PyMC, float PzMC, int value) -> float {
                             if (value == 0 || value == 1)
                               return RecoDecay::y(std::array{PxMC, PyMC, PzMC}, o2::constants::physics::MassXiMinus);
                             if (value == 2 || value == 3)
                               return RecoDecay::y(std::array{PxMC, PyMC, PzMC}, o2::constants::physics::MassOmegaMinus);
                             return 0.0f;
                           });

DECLARE_SOA_DYNAMIC_COLUMN(NegativePtMC, negativeptMC, //! negative daughter pT
                           [](float pxNegMC, float pyNegMC) -> float { return RecoDecay::sqrtSumOfSquares(pxNegMC, pyNegMC); });
DECLARE_SOA_DYNAMIC_COLUMN(PositivePtMC, positiveptMC, //! positive daughter pT
                           [](float pxPosMC, float pyPosMC) -> float { return RecoDecay::sqrtSumOfSquares(pxPosMC, pyPosMC); });
DECLARE_SOA_DYNAMIC_COLUMN(BachelorPtMC, bachelorptMC, //! bachelor daughter pT
                           [](float pxBachMC, float pyBachMC) -> float { return RecoDecay::sqrtSumOfSquares(pxBachMC, pyBachMC); });
DECLARE_SOA_DYNAMIC_COLUMN(PtMC, ptMC, //! cascade pT
                           [](float pxMC, float pyMC) -> float { return RecoDecay::sqrtSumOfSquares(pxMC, pyMC); });
} // namespace cascdata

//______________________________________________________
// EXPRESSION COLUMNS FOR TRACASCCORES
namespace cascdataext
{
DECLARE_SOA_EXPRESSION_COLUMN(PxLambda, pxlambda, //!
                              float, 1.f * aod::cascdata::pxpos + 1.f * aod::cascdata::pxneg);
DECLARE_SOA_EXPRESSION_COLUMN(PyLambda, pylambda, //!
                              float, 1.f * aod::cascdata::pypos + 1.f * aod::cascdata::pyneg);
DECLARE_SOA_EXPRESSION_COLUMN(PzLambda, pzlambda, //!
                              float, 1.f * aod::cascdata::pzpos + 1.f * aod::cascdata::pzneg);
DECLARE_SOA_EXPRESSION_COLUMN(Pt, pt, float, //! Transverse momentum in GeV/c
                              nsqrt(aod::cascdata::px* aod::cascdata::px +
                                    aod::cascdata::py * aod::cascdata::py));
DECLARE_SOA_EXPRESSION_COLUMN(P, p, float, //! Total momentum in GeV/c
                              nsqrt(aod::cascdata::px* aod::cascdata::px +
                                    aod::cascdata::py * aod::cascdata::py +
                                    aod::cascdata::pz * aod::cascdata::pz));
DECLARE_SOA_EXPRESSION_COLUMN(Phi, phi, float, //! Phi in the range [0, 2pi)
                              o2::constants::math::PI + natan2(-1.0f * aod::cascdata::py, -1.0f * aod::cascdata::px));
DECLARE_SOA_EXPRESSION_COLUMN(Eta, eta, float, //! Pseudorapidity, conditionally defined to avoid FPEs
                              ifnode((nsqrt(aod::cascdata::px * aod::cascdata::px +
                                            aod::cascdata::py * aod::cascdata::py +
                                            aod::cascdata::pz * aod::cascdata::pz) -
                                      aod::cascdata::pz) < static_cast<float>(1e-7),
                                     ifnode(aod::cascdata::pz < 0.f, -100.f, 100.f),
                                     0.5f * nlog((nsqrt(aod::cascdata::px * aod::cascdata::px +
                                                        aod::cascdata::py * aod::cascdata::py +
                                                        aod::cascdata::pz * aod::cascdata::pz) +
                                                  aod::cascdata::pz) /
                                                 (nsqrt(aod::cascdata::px * aod::cascdata::px +
                                                        aod::cascdata::py * aod::cascdata::py +
                                                        aod::cascdata::pz * aod::cascdata::pz) -
                                                  aod::cascdata::pz))));
} // namespace cascdataext

//______________________________________________________
// Cascade data model:
// --- standard, KF version and tracked version
// includes three variants to be able to account for
// cases in which multiple analyses are executed
// at the same time. Tables are kept as similar
// as possible for ease of use and comparison

DECLARE_SOA_TABLE(CascIndices, "AOD", "CascINDEX", //! index table when using AO2Ds
                  o2::soa::Index<>, cascdata::CascadeId, v0data::PosTrackId, v0data::NegTrackId, cascdata::BachelorId, cascdata::CollisionId, o2::soa::Marker<1>);
DECLARE_SOA_TABLE(KFCascIndices, "AOD", "KFCascINDEX", //! index table when using AO2Ds
                  o2::soa::Index<>, cascdata::CascadeId, v0data::PosTrackId, v0data::NegTrackId, cascdata::BachelorId, cascdata::CollisionId, o2::soa::Marker<2>);
DECLARE_SOA_TABLE(TraCascIndices, "AOD", "TraCascINDEX", //! index table when using AO2Ds
                  o2::soa::Index<>, cascdata::CascadeId, v0data::PosTrackId, v0data::NegTrackId, cascdata::BachelorId, cascdata::StrangeTrackId, cascdata::CollisionId);

DECLARE_SOA_TABLE(CascCollRefs, "AOD", "CASCCOLLREF", //! optional table to refer back to a collision
                  o2::soa::Index<>, cascdata::StraCollisionId, o2::soa::Marker<1>);
DECLARE_SOA_TABLE(KFCascCollRefs, "AOD", "KFCASCCOLLREF", //! optional table to refer back to a collision
                  o2::soa::Index<>, cascdata::StraCollisionId, o2::soa::Marker<2>);
DECLARE_SOA_TABLE(TraCascCollRefs, "AOD", "TRACASCCOLLREF", //! optional table to refer back to a collision
                  o2::soa::Index<>, cascdata::StraCollisionId, o2::soa::Marker<3>);

DECLARE_SOA_TABLE(CascExtras, "AOD", "CASCEXTRA", //! optional table to refer to custom track extras
                  o2::soa::Index<>, cascdata::PosTrackExtraId, cascdata::NegTrackExtraId,
                  cascdata::BachTrackExtraId, o2::soa::Marker<1>);
DECLARE_SOA_TABLE(StraTrackExtras, "AOD", "STRATRACKEXTRAS", //! optional table to refer to custom track extras
                  o2::soa::Index<>, cascdata::StrangeTrackExtraId);

// Track positions at minima (valid for regular cascade table), replayable
DECLARE_SOA_TABLE(CascTrackXs, "AOD", "CASCTRACKX", //! track X positions at minima when using AO2Ds
                  cascdata::PosX, cascdata::NegX, cascdata::BachX);

DECLARE_SOA_TABLE(StoredCascCores, "AOD", "CASCCORE", //! core information about decay, viable with AO2Ds or derived
                  cascdata::Sign, cascdata::MXi, cascdata::MOmega,
                  cascdata::X, cascdata::Y, cascdata::Z,
                  cascdata::Xlambda, cascdata::Ylambda, cascdata::Zlambda,
                  cascdata::PxPos, cascdata::PyPos, cascdata::PzPos,
                  cascdata::PxNeg, cascdata::PyNeg, cascdata::PzNeg,
                  cascdata::PxBach, cascdata::PyBach, cascdata::PzBach,
                  cascdata::Px, cascdata::Py, cascdata::Pz,
                  cascdata::DCAV0Daughters, cascdata::DCACascDaughters,
                  cascdata::DCAPosToPV, cascdata::DCANegToPV, cascdata::DCABachToPV, cascdata::DCAXYCascToPV, cascdata::DCAZCascToPV,

                  // Dynamic columns
                  cascdata::V0Radius<cascdata::Xlambda, cascdata::Ylambda>,
                  cascdata::CascRadius<cascdata::X, cascdata::Y>,
                  cascdata::V0CosPA<cascdata::Xlambda, cascdata::Ylambda, cascdata::Zlambda, cascdataext::PxLambda, cascdataext::PyLambda, cascdataext::PzLambda>,
                  cascdata::CascCosPA<cascdata::X, cascdata::Y, cascdata::Z, cascdata::PxBach, cascdata::PxPos, cascdata::PxNeg, cascdata::PyBach, cascdata::PyPos, cascdata::PyNeg, cascdata::PzBach, cascdata::PzPos, cascdata::PzNeg>,
                  cascdata::DCAV0ToPV<cascdata::Xlambda, cascdata::Ylambda, cascdata::Zlambda, cascdataext::PxLambda, cascdataext::PyLambda, cascdataext::PzLambda>,

                  // Invariant masses
                  cascdata::MLambda<cascdata::Sign, cascdata::PxPos, cascdata::PyPos, cascdata::PzPos, cascdata::PxNeg, cascdata::PyNeg, cascdata::PzNeg>,
                  cascdata::M<cascdata::MXi, cascdata::MOmega>,

                  // Longitudinal
                  cascdata::YXi<cascdata::Px, cascdata::Py, cascdata::Pz>,
                  cascdata::YOmega<cascdata::Px, cascdata::Py, cascdata::Pz>,
                  cascdata::Rapidity<cascdata::Px, cascdata::Py, cascdata::Pz>,

                  cascdata::NegativePt<cascdata::PxNeg, cascdata::PyNeg>,
                  cascdata::PositivePt<cascdata::PxPos, cascdata::PyPos>,
                  cascdata::BachelorPt<cascdata::PxBach, cascdata::PyBach>,
                  cascdata::NegativeEta<cascdata::PxNeg, cascdata::PyNeg, cascdata::PzNeg>,
                  cascdata::NegativePhi<cascdata::PxNeg, cascdata::PyNeg>,
                  cascdata::PositiveEta<cascdata::PxPos, cascdata::PyPos, cascdata::PzPos>,
                  cascdata::PositivePhi<cascdata::PxPos, cascdata::PyPos>,
                  cascdata::BachelorEta<cascdata::PxBach, cascdata::PyBach, cascdata::PzBach>,
                  cascdata::BachelorPhi<cascdata::PxBach, cascdata::PyBach>);

DECLARE_SOA_TABLE(StoredKFCascCores, "AOD", "KFCASCCORE", //!
                  cascdata::Sign, cascdata::MXi, cascdata::MOmega,
                  cascdata::X, cascdata::Y, cascdata::Z,
                  cascdata::Xlambda, cascdata::Ylambda, cascdata::Zlambda,
                  cascdata::KFXPos, cascdata::KFYPos, cascdata::KFZPos,
                  cascdata::KFXNeg, cascdata::KFYNeg, cascdata::KFZNeg,
                  cascdata::PxPos, cascdata::PyPos, cascdata::PzPos,
                  cascdata::PxNeg, cascdata::PyNeg, cascdata::PzNeg,
                  cascdata::PxBach, cascdata::PyBach, cascdata::PzBach,
                  cascdata::KFPxV0, cascdata::KFPyV0, cascdata::KFPzV0,
                  cascdata::Px, cascdata::Py, cascdata::Pz,
                  cascdata::DCAV0Daughters, cascdata::DCACascDaughters,
                  cascdata::DCAPosToPV, cascdata::DCANegToPV, cascdata::DCABachToPV, cascdata::DCAXYCascToPV, cascdata::DCAZCascToPV,

                  // KF particle fit specific
                  kfcascdata::MLambda, cascdata::KFV0Chi2, cascdata::KFCascadeChi2,

                  // Dynamic columns
                  cascdata::V0Radius<cascdata::Xlambda, cascdata::Ylambda>,
                  cascdata::CascRadius<cascdata::X, cascdata::Y>,
                  cascdata::V0CosPA<cascdata::Xlambda, cascdata::Ylambda, cascdata::Zlambda, cascdataext::PxLambda, cascdataext::PyLambda, cascdataext::PzLambda>,
                  cascdata::CascCosPA<cascdata::X, cascdata::Y, cascdata::Z, cascdata::PxBach, cascdata::PxPos, cascdata::PxNeg, cascdata::PyBach, cascdata::PyPos, cascdata::PyNeg, cascdata::PzBach, cascdata::PzPos, cascdata::PzNeg>,
                  cascdata::DCAV0ToPV<cascdata::Xlambda, cascdata::Ylambda, cascdata::Zlambda, cascdataext::PxLambda, cascdataext::PyLambda, cascdataext::PzLambda>,

                  // Invariant masses
                  cascdata::M<cascdata::MXi, cascdata::MOmega>,

                  // Longitudinal
                  cascdata::YXi<cascdata::Px, cascdata::Py, cascdata::Pz>,
                  cascdata::YOmega<cascdata::Px, cascdata::Py, cascdata::Pz>,

                  cascdata::NegativePt<cascdata::PxNeg, cascdata::PyNeg>,
                  cascdata::PositivePt<cascdata::PxPos, cascdata::PyPos>,
                  cascdata::BachelorPt<cascdata::PxBach, cascdata::PyBach>,
                  cascdata::NegativeEta<cascdata::PxNeg, cascdata::PyNeg, cascdata::PzNeg>,
                  cascdata::NegativePhi<cascdata::PxNeg, cascdata::PyNeg>,
                  cascdata::PositiveEta<cascdata::PxPos, cascdata::PyPos, cascdata::PzPos>,
                  cascdata::PositivePhi<cascdata::PxPos, cascdata::PyPos>,
                  cascdata::BachelorEta<cascdata::PxBach, cascdata::PyBach, cascdata::PzBach>,
                  cascdata::BachelorPhi<cascdata::PxBach, cascdata::PyBach>);

DECLARE_SOA_TABLE(StoredTraCascCores, "AOD", "TRACASCCORE", //!
                  cascdata::Sign, cascdata::MXi, cascdata::MOmega,
                  cascdata::X, cascdata::Y, cascdata::Z,
                  cascdata::Xlambda, cascdata::Ylambda, cascdata::Zlambda,
                  cascdata::PxPos, cascdata::PyPos, cascdata::PzPos,
                  cascdata::PxNeg, cascdata::PyNeg, cascdata::PzNeg,
                  cascdata::PxBach, cascdata::PyBach, cascdata::PzBach,
                  cascdata::Px, cascdata::Py, cascdata::Pz,
                  cascdata::DCAV0Daughters, cascdata::DCACascDaughters,
                  cascdata::DCAPosToPV, cascdata::DCANegToPV, cascdata::DCABachToPV, cascdata::DCAXYCascToPV, cascdata::DCAZCascToPV,

                  // Strangeness tracking specific
                  cascdata::MatchingChi2, cascdata::TopologyChi2, cascdata::ItsClsSize,

                  // Dynamic columns
                  cascdata::V0Radius<cascdata::Xlambda, cascdata::Ylambda>,
                  cascdata::CascRadius<cascdata::X, cascdata::Y>,
                  cascdata::V0CosPA<cascdata::Xlambda, cascdata::Ylambda, cascdata::Zlambda, cascdataext::PxLambda, cascdataext::PyLambda, cascdataext::PzLambda>,
                  cascdata::CascCosPA<cascdata::X, cascdata::Y, cascdata::Z, cascdata::PxBach, cascdata::PxPos, cascdata::PxNeg, cascdata::PyBach, cascdata::PyPos, cascdata::PyNeg, cascdata::PzBach, cascdata::PzPos, cascdata::PzNeg>,
                  cascdata::DCAV0ToPV<cascdata::Xlambda, cascdata::Ylambda, cascdata::Zlambda, cascdataext::PxLambda, cascdataext::PyLambda, cascdataext::PzLambda>,

                  // Invariant masses
                  cascdata::MLambda<cascdata::Sign, cascdata::PxPos, cascdata::PyPos, cascdata::PzPos, cascdata::PxNeg, cascdata::PyNeg, cascdata::PzNeg>,

                  // Longitudinal
                  cascdata::YXi<cascdata::Px, cascdata::Py, cascdata::Pz>,
                  cascdata::YOmega<cascdata::Px, cascdata::Py, cascdata::Pz>,

                  cascdata::NegativePt<cascdata::PxNeg, cascdata::PyNeg>,
                  cascdata::PositivePt<cascdata::PxPos, cascdata::PyPos>,
                  cascdata::BachelorPt<cascdata::PxBach, cascdata::PyBach>,
                  cascdata::NegativeEta<cascdata::PxNeg, cascdata::PyNeg, cascdata::PzNeg>,
                  cascdata::NegativePhi<cascdata::PxNeg, cascdata::PyNeg>,
                  cascdata::PositiveEta<cascdata::PxPos, cascdata::PyPos, cascdata::PzPos>,
                  cascdata::PositivePhi<cascdata::PxPos, cascdata::PyPos>,
                  cascdata::BachelorEta<cascdata::PxBach, cascdata::PyBach, cascdata::PzBach>,
                  cascdata::BachelorPhi<cascdata::PxBach, cascdata::PyBach>);

DECLARE_SOA_TABLE(CascMCCores, "AOD", "CASCMCCORE", //! bachelor-baryon correlation variables
                  cascdata::PDGCode, cascdata::PDGCodeMother, cascdata::PDGCodeV0, cascdata::IsPhysicalPrimary,
                  cascdata::PDGCodePositive, cascdata::PDGCodeNegative, cascdata::PDGCodeBachelor,
                  cascdata::XMC, cascdata::YMC, cascdata::ZMC,
                  cascdata::XlambdaMC, cascdata::YlambdaMC, cascdata::ZlambdaMC,
                  cascdata::PxPosMC, cascdata::PyPosMC, cascdata::PzPosMC,
                  cascdata::PxNegMC, cascdata::PyNegMC, cascdata::PzNegMC,
                  cascdata::PxBachMC, cascdata::PyBachMC, cascdata::PzBachMC,
                  cascdata::PxMC, cascdata::PyMC, cascdata::PzMC,
                  cascdata::RapidityMC<cascdata::PxMC, cascdata::PyMC, cascdata::PzMC>,
                  cascdata::NegativePtMC<cascdata::PxNegMC, cascdata::PyNegMC>,
                  cascdata::PositivePtMC<cascdata::PxPosMC, cascdata::PyPosMC>,
                  cascdata::BachelorPtMC<cascdata::PxBachMC, cascdata::PyBachMC>,
                  cascdata::PtMC<cascdata::PxMC, cascdata::PyMC>);

namespace cascdata
{
DECLARE_SOA_INDEX_COLUMN(CascMCCore, cascMCCore); //! Index to CascMCCore entry
}

DECLARE_SOA_TABLE(CascCoreMCLabels, "AOD", "CASCCOREMCLABEL", //! optional table to refer to CascMCCores if not joinable
                  o2::soa::Index<>, cascdata::CascMCCoreId);
DECLARE_SOA_TABLE(CascMCCollRefs, "AOD", "CASCMCCOLLREF", //! refers MC candidate back to proper MC Collision
                  o2::soa::Index<>, cascdata::StraMCCollisionId, o2::soa::Marker<3>);

DECLARE_SOA_TABLE(GeXiMinus, "AOD", "GeXiMinus", cascdata::GeneratedXiMinus);
DECLARE_SOA_TABLE(GeXiPlus, "AOD", "GeXiPlus", cascdata::GeneratedXiPlus);
DECLARE_SOA_TABLE(GeOmegaMinus, "AOD", "GeOmegaMinus", cascdata::GeneratedOmegaMinus);
DECLARE_SOA_TABLE(GeOmegaPlus, "AOD", "GeOmegaPlus", cascdata::GeneratedOmegaPlus);

DECLARE_SOA_TABLE(CascMCMothers, "AOD", "CASCMCMOTHER", //! optional table for MC mothers
                  o2::soa::Index<>, cascdata::MotherMCPartId);

DECLARE_SOA_TABLE(CascBBs, "AOD", "CASCBB", //! bachelor-baryon correlation variables
                  cascdata::BachBaryonCosPA, cascdata::BachBaryonDCAxyToPV)

DECLARE_SOA_TABLE(CascCovs, "AOD", "CASCCOVS", //!
                  cascdata::CovMat,
                  cascdata::PositionCovMat<cascdata::CovMat>,
                  cascdata::MomentumCovMat<cascdata::CovMat>,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(KFCascCovs, "AOD", "KFCASCCOVS", //!
                  cascdata::KFTrackCovMat, cascdata::KFTrackCovMatV0, cascdata::KFTrackCovMatV0DauPos, cascdata::KFTrackCovMatV0DauNeg);

DECLARE_SOA_TABLE(TraCascCovs, "AOD", "TRACASCCOVS", //!
                  cascdata::CovMat,
                  cascdata::PositionCovMat<cascdata::CovMat>,
                  cascdata::MomentumCovMat<cascdata::CovMat>,
                  o2::soa::Marker<2>);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(CascCores, StoredCascCores, "CascDATAEXT", //!
                                cascdataext::PxLambda, cascdataext::PyLambda, cascdataext::PzLambda, cascdataext::Pt, cascdataext::P, cascdataext::Eta, cascdataext::Phi);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(KFCascCores, StoredKFCascCores, "KFCascDATAEXT", //!
                                cascdataext::PxLambda, cascdataext::PyLambda, cascdataext::PzLambda,
                                cascdataext::Pt, cascdataext::P, cascdataext::Eta, cascdataext::Phi);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(TraCascCores, StoredTraCascCores, "TraCascDATAEXT", //!
                                cascdataext::PxLambda, cascdataext::PyLambda, cascdataext::PzLambda,
                                cascdataext::Pt, cascdataext::P, cascdataext::Eta, cascdataext::Phi);

namespace cascdata
{
// For cross-linking all cascade kinds
DECLARE_SOA_INDEX_COLUMN_FULL(TrackedCascade, trackedCascade, int, TraCascCores, "_Refs"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(KFCascade, kfCascade, int, KFCascCores, "_Refs");            //!
DECLARE_SOA_INDEX_COLUMN_FULL(StandardCascade, standardCascade, int, CascCores, "_Refs");  //!
} // namespace cascdata

// interlink different cascade types
DECLARE_SOA_TABLE(CascToTraRefs, "AOD", "CASCTOTRAREFS", //! standard -> tracked
                  o2::soa::Index<>, cascdata::TrackedCascadeId);
DECLARE_SOA_TABLE(CascToKFRefs, "AOD", "CASCTOKFREFS", //! standard -> KF
                  o2::soa::Index<>, cascdata::KFCascadeId);
DECLARE_SOA_TABLE(TraToCascRefs, "AOD", "TRATOCASCREFS", //! standard -> KF
                  o2::soa::Index<>, cascdata::StandardCascadeId, o2::soa::Marker<1>);
DECLARE_SOA_TABLE(KFToCascRefs, "AOD", "KFTOCASCREFS", //! standard -> KF
                  o2::soa::Index<>, cascdata::StandardCascadeId, o2::soa::Marker<2>);

using CascIndex = CascIndices::iterator;
using CascCore = CascCores::iterator;
using KFCascIndex = KFCascIndices::iterator;
using KFCascCore = KFCascCores::iterator;
using TraCascIndex = TraCascIndices::iterator;
using TraCascCore = TraCascCores::iterator;

using CascDatas = soa::Join<CascIndices, CascBBs, CascCores>;
using KFCascDatas = soa::Join<KFCascIndices, KFCascCores>;
using TraCascDatas = soa::Join<TraCascIndices, TraCascCores>;

using CascData = CascDatas::iterator;
using KFCascData = KFCascDatas::iterator;
using TraCascData = TraCascDatas::iterator;

using CascMCCore = CascMCCores::iterator;
using CascMCMother = CascMCMothers::iterator;
using CascMCDatas = soa::Join<CascMCCores, CascMCMothers>;
using CascMCData = CascMCDatas::iterator;

// For compatibility with previous table declarations
using CascDataFull = CascDatas;
using CascDataExt = CascDatas;

namespace cascdata
{
DECLARE_SOA_INDEX_COLUMN(CascData, cascData);       //! Index to CascData entry
DECLARE_SOA_INDEX_COLUMN(KFCascData, kfCascData);   //! Index to CascData entry
DECLARE_SOA_INDEX_COLUMN(TraCascData, traCascData); //! Index to CascData entry
} // namespace cascdata

DECLARE_SOA_TABLE(CascDataLink, "AOD", "CASCDATALINK", //! Joinable table with Cascades which links to CascData which is not produced for all entries
                  cascdata::CascDataId);
DECLARE_SOA_TABLE(KFCascDataLink, "AOD", "KFCASCDATALINK", //! Joinable table with Cascades which links to CascData which is not produced for all entries
                  cascdata::KFCascDataId);
DECLARE_SOA_TABLE(TraCascDataLink, "AOD", "TRACASCDATALINK", //! Joinable table with Cascades which links to CascData which is not produced for all entries
                  cascdata::TraCascDataId);

using CascadesLinked = soa::Join<Cascades, CascDataLink>;
using CascadeLinked = CascadesLinked::iterator;
using KFCascadesLinked = soa::Join<Cascades, KFCascDataLink>;
using KFCascadeLinked = KFCascadesLinked::iterator;
using TraCascadesLinked = soa::Join<Cascades, TraCascDataLink>;
using TraCascadeLinked = TraCascadesLinked::iterator;

namespace cascdata
{
DECLARE_SOA_INDEX_COLUMN(FindableV0, findableV0); //! V0 index
DECLARE_SOA_COLUMN(IsFound, isFound, bool);       //! is this FindableCascade actually in the Cascades table?
} // namespace cascdata

DECLARE_SOA_TABLE(FindableCascades, "AOD", "FINDABLECASCS", //! Run 3 cascade table
                  o2::soa::Index<>, cascade::CollisionId, cascdata::FindableV0Id, cascade::BachelorId, o2::soa::Marker<1>);

DECLARE_SOA_TABLE(CascFoundTags, "AOD", "CascFoundTag", //! found or not?
                  cascdata::IsFound);

using FindableCascadesLinked = soa::Join<FindableCascades, CascDataLink>;
using FindableCascadeLinked = FindableCascadesLinked::iterator;

namespace casctag
{
DECLARE_SOA_COLUMN(IsInteresting, isInteresting, bool); //! will this be built or not?

// MC association bools
DECLARE_SOA_COLUMN(IsTrueXiMinus, isTrueXiMinus, bool);         //! PDG checked correctly in MC
DECLARE_SOA_COLUMN(IsTrueXiPlus, isTrueXiPlus, bool);           //! PDG checked correctly in MC
DECLARE_SOA_COLUMN(IsTrueOmegaMinus, isTrueOmegaMinus, bool);   //! PDG checked correctly in MC
DECLARE_SOA_COLUMN(IsTrueOmegaPlus, isTrueOmegaPlus, bool);     //! PDG checked correctly in MC
DECLARE_SOA_COLUMN(IsPhysicalPrimary, isPhysicalPrimary, bool); //! physical primary

// dE/dx compatibility bools
DECLARE_SOA_COLUMN(IsdEdxXiMinus, isdEdxXiMinus, bool);       //! compatible with dE/dx hypotheses
DECLARE_SOA_COLUMN(IsdEdxXiPlus, isdEdxXiPlus, bool);         //! compatible with dE/dx hypotheses
DECLARE_SOA_COLUMN(IsdEdxOmegaMinus, isdEdxOmegaMinus, bool); //! compatible with dE/dx hypotheses
DECLARE_SOA_COLUMN(IsdEdxOmegaPlus, isdEdxOmegaPlus, bool);   //! compatible with dE/dx hypotheses
} // namespace casctag
DECLARE_SOA_TABLE(CascTags, "AOD", "CASCTAGS",
                  casctag::IsInteresting,
                  casctag::IsTrueXiMinus,
                  casctag::IsTrueXiPlus,
                  casctag::IsTrueOmegaMinus,
                  casctag::IsTrueOmegaPlus,
                  casctag::IsPhysicalPrimary,
                  casctag::IsdEdxXiMinus,
                  casctag::IsdEdxXiPlus,
                  casctag::IsdEdxOmegaMinus,
                  casctag::IsdEdxOmegaPlus);

// Definition of labels for V0s
namespace mcv0label
{
DECLARE_SOA_INDEX_COLUMN(McParticle, mcParticle);                                               //! MC particle for V0
DECLARE_SOA_INDEX_COLUMN_FULL(McMotherParticle, mcMotherParticle, int, McParticles, "_Mother"); //!
} // namespace mcv0label

DECLARE_SOA_TABLE(McV0Labels, "AOD", "MCV0LABEL", //! Table joinable with V0Data containing the MC labels
                  mcv0label::McParticleId, mcv0label::McMotherParticleId);
using McV0Label = McV0Labels::iterator;

// Definition of labels for V0s // Full table, joinable with V0 (CAUTION: NOT WITH V0DATA)
namespace mcfullv0label
{
DECLARE_SOA_INDEX_COLUMN(McParticle, mcParticle); //! MC particle for V0
} // namespace mcfullv0label

DECLARE_SOA_TABLE(McFullV0Labels, "AOD", "MCFULLV0LABEL", //! Table joinable with V0
                  mcfullv0label::McParticleId);
using McFullV0Label = McFullV0Labels::iterator;

// Definition of labels for cascades
namespace mccasclabel
{
DECLARE_SOA_INDEX_COLUMN(McParticle, mcParticle);                                               //! MC particle for Cascade
DECLARE_SOA_INDEX_COLUMN_FULL(McMotherParticle, mcMotherParticle, int, McParticles, "_Mother"); //!
DECLARE_SOA_COLUMN(IsBachBaryonCandidate, isBachBaryonCandidate, bool);                         //! will this be built or not?
} // namespace mccasclabel

DECLARE_SOA_TABLE(McCascLabels, "AOD", "MCCASCLABEL", //! Table joinable with CascData containing the MC labels
                  mccasclabel::McParticleId, mccasclabel::McMotherParticleId);
DECLARE_SOA_TABLE(McCascBBTags, "AOD", "MCCASCBBTAG", //! Table joinable with CascData containing yes / no for BB correlation
                  mccasclabel::IsBachBaryonCandidate);
using McCascLabel = McCascLabels::iterator;
using McCascBBTag = McCascBBTags::iterator;

// Definition of labels for KF cascades
namespace mckfcasclabel
{
DECLARE_SOA_INDEX_COLUMN(McParticle, mcParticle); //! MC particle for KF Cascade
} // namespace mckfcasclabel

DECLARE_SOA_TABLE(McKFCascLabels, "AOD", "MCKFCASCLABEL", //! Table joinable with KFCascData containing the MC labels
                  mckfcasclabel::McParticleId);
using McKFCascLabel = McKFCascLabels::iterator;

// Definition of labels for tracked cascades
namespace mctracasclabel
{
DECLARE_SOA_INDEX_COLUMN(McParticle, mcParticle); //! MC particle for V0
} // namespace mctracasclabel

DECLARE_SOA_TABLE(McTraCascLabels, "AOD", "MCTRACASCLABEL", //! Table joinable to cascdata containing the MC labels
                  mctracasclabel::McParticleId);
using McTraCascLabel = McTraCascLabels::iterator;

DECLARE_SOA_TABLE(TrackedCascadeColls, "AOD", "TRACASCCOLL", //! Table joinable with TrackedCascades containing collision ids
                  track::CollisionId, o2::soa::Marker<1>);
using TrackedCascadeColl = TrackedCascadeColls::iterator;
using AssignedTrackedCascades = soa::Join<aod::TrackedCascades, aod::TrackedCascadeColls>;
using AssignedTrackedCascade = AssignedTrackedCascades::iterator;

DECLARE_SOA_TABLE(TrackedV0Colls, "AOD", "TRAV0COLL", //! Table joinable with TrackedV0s containing collision ids
                  track::CollisionId, o2::soa::Marker<2>);
using TrackedV0Coll = TrackedV0Colls::iterator;
using AssignedTrackedV0s = soa::Join<aod::TrackedV0s, aod::TrackedV0Colls>;
using AssignedTrackedV0 = AssignedTrackedV0s::iterator;

DECLARE_SOA_TABLE(Tracked3BodyColls, "AOD", "TRA3BODYCOLL", //! Table joinable with Tracked3Bodys containing collision ids
                  track::CollisionId, o2::soa::Marker<3>);
using Tracked3BodyColl = Tracked3BodyColls::iterator;
using AssignedTracked3Bodys = soa::Join<aod::Tracked3Bodys, aod::Tracked3BodyColls>;
using AssignedTracked3Body = AssignedTracked3Bodys::iterator;
} // namespace o2::aod

//______________________________________________________
// Equivalency declarations
namespace o2::soa
{
DECLARE_EQUIVALENT_FOR_INDEX(aod::V0Indices, aod::V0Cores);
DECLARE_EQUIVALENT_FOR_INDEX(aod::V0TrackXs, aod::V0Cores);
DECLARE_EQUIVALENT_FOR_INDEX(aod::V0TrackXs, aod::V0Indices);
DECLARE_EQUIVALENT_FOR_INDEX(aod::CascIndices, aod::CascCores);
DECLARE_EQUIVALENT_FOR_INDEX(aod::CascIndices, aod::CascBBs);
DECLARE_EQUIVALENT_FOR_INDEX(aod::CascCores, aod::CascBBs);
DECLARE_EQUIVALENT_FOR_INDEX(aod::KFCascIndices, aod::KFCascCores);
DECLARE_EQUIVALENT_FOR_INDEX(aod::TraCascIndices, aod::TraCascCores);
} // namespace o2::soa

#endif // PWGLF_DATAMODEL_LFSTRANGENESSTABLES_H_
