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

#ifndef PWGDQ_DATAMODEL_REDUCEDINFOTABLES_H_
#define PWGDQ_DATAMODEL_REDUCEDINFOTABLES_H_

#include "PWGHF/Utils/utilsPid.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Qvectors.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "MathUtils/Utils.h"

#include <cmath>
#include <vector>

namespace o2::aod
{

namespace dqppfilter
{
DECLARE_SOA_COLUMN(EventFilter, eventFilter, uint64_t); //! Bit-field used for the high level event triggering
DECLARE_SOA_COLUMN(NewBcIndex, newBcIndex, uint64_t);   //! globalIndex of the new BC determined in filterPbPb
} // namespace dqppfilter

DECLARE_SOA_TABLE(DQEventFilter, "AOD", "EVENTFILTER", //! Store event-level decisions (DQ high level triggers)
                  dqppfilter::EventFilter);

DECLARE_SOA_TABLE(DQRapidityGapFilter, "AOD", "RAPIDITYGAPFILTER",
                  dqppfilter::EventFilter,
                  dqppfilter::NewBcIndex);

namespace reducedevent
{

// basic event information
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                                  //!
DECLARE_SOA_BITMAP_COLUMN(Tag, tag, 64);                                         //!  Bit-field for storing event information (e.g. high level info, cut decisions)
DECLARE_SOA_COLUMN(MCPosX, mcPosX, float);                                       //!  MC event position X
DECLARE_SOA_COLUMN(MCPosY, mcPosY, float);                                       //!  MC event position Y
DECLARE_SOA_COLUMN(MCPosZ, mcPosZ, float);                                       //!  MC event position Z
DECLARE_SOA_COLUMN(NTPCoccupContribLongA, nTPCoccupContribLongA, int);           //!  TPC pileup occupancy on A side (long time range)
DECLARE_SOA_COLUMN(NTPCoccupContribLongC, nTPCoccupContribLongC, int);           //!  TPC pileup occupancy on C side (long time range)
DECLARE_SOA_COLUMN(NTPCoccupMeanTimeLongA, nTPCoccupMeanTimeLongA, float);       //!  TPC pileup mean time on A side (long time range)
DECLARE_SOA_COLUMN(NTPCoccupMeanTimeLongC, nTPCoccupMeanTimeLongC, float);       //!  TPC pileup mean time on C side (long time range)
DECLARE_SOA_COLUMN(NTPCoccupMedianTimeLongA, nTPCoccupMedianTimeLongA, float);   //!  TPC pileup median time on A side (long time range)
DECLARE_SOA_COLUMN(NTPCoccupMedianTimeLongC, nTPCoccupMedianTimeLongC, float);   //!  TPC pileup median time on C side (long time range)
DECLARE_SOA_COLUMN(NTPCoccupContribShortA, nTPCoccupContribShortA, int);         //!  TPC pileup occupancy on A side (short time range)
DECLARE_SOA_COLUMN(NTPCoccupContribShortC, nTPCoccupContribShortC, int);         //!  TPC pileup occupancy on C side (short time range)
DECLARE_SOA_COLUMN(NTPCoccupMeanTimeShortA, nTPCoccupMeanTimeShortA, float);     //!  TPC pileup mean time on A side (short time range)
DECLARE_SOA_COLUMN(NTPCoccupMeanTimeShortC, nTPCoccupMeanTimeShortC, float);     //!  TPC pileup mean time on C side (short time range)
DECLARE_SOA_COLUMN(NTPCoccupMedianTimeShortA, nTPCoccupMedianTimeShortA, float); //!  TPC pileup median time on A side (short time range)
DECLARE_SOA_COLUMN(NTPCoccupMedianTimeShortC, nTPCoccupMedianTimeShortC, float); //!  TPC pileup median time on C side (short time range)
DECLARE_SOA_COLUMN(DCAzBimodalityCoefficient, dcazBimodalityCoefficient, float); //!  Bimodality coefficient of the DCAz distribution of the tracks in the event
DECLARE_SOA_COLUMN(DCAzBimodalityCoefficientBinned, dcazBimodalityCoefficientBinned, float);                 //!  Bimodality coefficient of the DCAz distribution of the tracks in the event, binned
DECLARE_SOA_COLUMN(DCAzBimodalityCoefficientBinnedTrimmed1, dcazBimodalityCoefficientBinnedTrimmed1, float); //!  Bimodality coefficient of the DCAz distribution of the tracks in the event, binned and trimmed 1
DECLARE_SOA_COLUMN(DCAzBimodalityCoefficientBinnedTrimmed2, dcazBimodalityCoefficientBinnedTrimmed2, float); //!  Bimodality coefficient of the DCAz distribution of the tracks in the event, binned and trimmed 2
DECLARE_SOA_COLUMN(DCAzBimodalityCoefficientBinnedTrimmed3, dcazBimodalityCoefficientBinnedTrimmed3, float); //!  Bimodality coefficient of the DCAz distribution of the tracks in the event, binned and trimmed 3
DECLARE_SOA_COLUMN(DCAzMean, dcazMean, float);                                                               //!  Mean of the DCAz distribution of the tracks in the event
DECLARE_SOA_COLUMN(DCAzMeanBinnedTrimmed1, dcazMeanBinnedTrimmed1, float);                                   //!  Mean of the DCAz distribution of the tracks in the event, binned and trimmed 1
DECLARE_SOA_COLUMN(DCAzMeanBinnedTrimmed2, dcazMeanBinnedTrimmed2, float);                                   //!  Mean of the DCAz distribution of the tracks in the event, binned and trimmed 2
DECLARE_SOA_COLUMN(DCAzMeanBinnedTrimmed3, dcazMeanBinnedTrimmed3, float);                                   //!  Mean of the DCAz distribution of the tracks in the event, binned and trimmed 3
DECLARE_SOA_COLUMN(DCAzRMS, dcazRMS, float);                                                                 //!  RMS of the DCAz distribution of the tracks in the event
DECLARE_SOA_COLUMN(DCAzRMSBinnedTrimmed1, dcazRMSBinnedTrimmed1, float);                                     //!  RMS of the DCAz distribution of the tracks in the event, binned and trimmed 1
DECLARE_SOA_COLUMN(DCAzRMSBinnedTrimmed2, dcazRMSBinnedTrimmed2, float);                                     //!  RMS of the DCAz distribution of the tracks in the event, binned and trimmed 2
DECLARE_SOA_COLUMN(DCAzRMSBinnedTrimmed3, dcazRMSBinnedTrimmed3, float);                                     //!  RMS of the DCAz distribution of the tracks in the event, binned and trimmed 3
DECLARE_SOA_COLUMN(DCAzSkewness, dcazSkewness, float);                                                       //!  Skewness of the DCAz distribution of the tracks in the event
DECLARE_SOA_COLUMN(DCAzKurtosis, dcazKurtosis, float);                                                       //!  Kurtosis of the DCAz distribution of the tracks in the event
DECLARE_SOA_COLUMN(DCAzFracAbove100um, dcazFracAbove100um, float);                                           //!  Fraction of tracks in the event with |DCAz| > 100um
DECLARE_SOA_COLUMN(DCAzFracAbove200um, dcazFracAbove200um, float);                                           //!  Fraction of tracks in the event with |DCAz| > 200um
DECLARE_SOA_COLUMN(DCAzFracAbove500um, dcazFracAbove500um, float);                                           //!  Fraction of tracks in the event with |DCAz| > 500um
DECLARE_SOA_COLUMN(DCAzFracAbove1mm, dcazFracAbove1mm, float);                                               //!  Fraction of tracks in the event with |DCAz| > 1mm
DECLARE_SOA_COLUMN(DCAzFracAbove2mm, dcazFracAbove2mm, float);                                               //!  Fraction of tracks in the event with |DCAz| > 2mm
DECLARE_SOA_COLUMN(DCAzFracAbove5mm, dcazFracAbove5mm, float);                                               //!  Fraction of tracks in the event with |DCAz| > 5mm
DECLARE_SOA_COLUMN(DCAzFracAbove10mm, dcazFracAbove10mm, float);                                             //!  Fraction of tracks in the event with |DCAz| > 10mm
DECLARE_SOA_COLUMN(DCAzNPeaks, dcazNPeaks, int);                                                             //!  Number of peaks in the DCAz distribution of the tracks in the event
DECLARE_SOA_COLUMN(DCAzNPeaksTrimmed1, dcazNPeaksTrimmed1, int);                                             //!  Number of peaks in the binned DCAz distribution (trimmed 1)
DECLARE_SOA_COLUMN(DCAzNPeaksTrimmed2, dcazNPeaksTrimmed2, int);                                             //!  Number of peaks in the binned DCAz distribution (trimmed 2)
DECLARE_SOA_COLUMN(DCAzNPeaksTrimmed3, dcazNPeaksTrimmed3, int);                                             //!  Number of peaks in the binned DCAz distribution (trimmed 3)

// Columns declared to guarantee the backward compatibility of the tables
DECLARE_SOA_COLUMN(QvecBPosRe, qvecBPosRe, float);
DECLARE_SOA_COLUMN(QvecBPosIm, qvecBPosIm, float);
DECLARE_SOA_COLUMN(QvecBNegRe, qvecBNegRe, float);
DECLARE_SOA_COLUMN(QvecBNegIm, qvecBNegIm, float);
DECLARE_SOA_COLUMN(QvecBAllRe, qvecBAllRe, float);
DECLARE_SOA_COLUMN(QvecBAllIm, qvecBAllIm, float);
DECLARE_SOA_COLUMN(NTrkBPos, nTrkBPos, int);
DECLARE_SOA_COLUMN(NTrkBNeg, nTrkBNeg, int);
DECLARE_SOA_COLUMN(NTrkBAll, nTrkBAll, int);

DECLARE_SOA_COLUMN(Q1ZNAX, q1znax, float);                 //!  Q-vector x component, evaluated with ZNA (harmonic 1 and power 1)
DECLARE_SOA_COLUMN(Q1ZNAY, q1znay, float);                 //!  Q-vector y component, evaluated with ZNA (harmonic 1 and power 1)
DECLARE_SOA_COLUMN(Q1ZNCX, q1zncx, float);                 //!  Q-vector x component, evaluated with ZNC (harmonic 1 and power 1)
DECLARE_SOA_COLUMN(Q1ZNCY, q1zncy, float);                 //!  Q-vector y component, evaluated with ZNC (harmonic 1 and power 1)
DECLARE_SOA_COLUMN(Q1X0A, q1x0a, float);                   //!  Q-vector x component, with event eta gap A (harmonic 1 and power 1)
DECLARE_SOA_COLUMN(Q1Y0A, q1y0a, float);                   //!  Q-vector y component, with event eta gap A (harmonic 1 and power 1)
DECLARE_SOA_COLUMN(Q1X0B, q1x0b, float);                   //!  Q-vector x component, with event eta gap B (harmonic 1 and power 1)
DECLARE_SOA_COLUMN(Q1Y0B, q1y0b, float);                   //!  Q-vector y component, with event eta gap B (harmonic 1 and power 1)
DECLARE_SOA_COLUMN(Q1X0C, q1x0c, float);                   //!  Q-vector x component, with event eta gap C (harmonic 1 and power 1)
DECLARE_SOA_COLUMN(Q1Y0C, q1y0c, float);                   //!  Q-vector y component, with event eta gap C (harmonic 1 and power 1)
DECLARE_SOA_COLUMN(Q2X0A, q2x0a, float);                   //!  Q-vector x component, with event eta gap A (harmonic 2 and power 1)
DECLARE_SOA_COLUMN(Q2Y0A, q2y0a, float);                   //!  Q-vector y component, with event eta gap A (harmonic 2 and power 1)
DECLARE_SOA_COLUMN(Q2X0B, q2x0b, float);                   //!  Q-vector x component, with event eta gap B (harmonic 2 and power 1)
DECLARE_SOA_COLUMN(Q2Y0B, q2y0b, float);                   //!  Q-vector y component, with event eta gap B (harmonic 2 and power 1)
DECLARE_SOA_COLUMN(Q2X0C, q2x0c, float);                   //!  Q-vector x component, with event eta gap C (harmonic 2 and power 1)
DECLARE_SOA_COLUMN(Q2Y0C, q2y0c, float);                   //!  Q-vector y component, with event eta gap C (harmonic 2 and power 1)
DECLARE_SOA_COLUMN(MultA, multa, float);                   //!  Event multiplicity eta gap A
DECLARE_SOA_COLUMN(MultB, multb, float);                   //!  Event multiplicity eta gap B
DECLARE_SOA_COLUMN(MultC, multc, float);                   //!  Event multiplicity eta gap C
DECLARE_SOA_COLUMN(Q3X0A, q3x0a, float);                   //!  Q-vector x component, with event eta gap A (harmonic 3 and power 1)
DECLARE_SOA_COLUMN(Q3Y0A, q3y0a, float);                   //!  Q-vector y component, with event eta gap A (harmonic 3 and power 1)
DECLARE_SOA_COLUMN(Q3X0B, q3x0b, float);                   //!  Q-vector x component, with event eta gap B (harmonic 3 and power 1)
DECLARE_SOA_COLUMN(Q3Y0B, q3y0b, float);                   //!  Q-vector y component, with event eta gap B (harmonic 3 and power 1)
DECLARE_SOA_COLUMN(Q3X0C, q3x0c, float);                   //!  Q-vector x component, with event eta gap C (harmonic 3 and power 1)
DECLARE_SOA_COLUMN(Q3Y0C, q3y0c, float);                   //!  Q-vector y component, with event eta gap C (harmonic 3 and power 1)
DECLARE_SOA_COLUMN(Q4X0A, q4x0a, float);                   //!  Q-vector x component, with event eta gap A (harmonic 4 and power 1)
DECLARE_SOA_COLUMN(Q4Y0A, q4y0a, float);                   //!  Q-vector y component, with event eta gap A (harmonic 4 and power 1)
DECLARE_SOA_COLUMN(Q4X0B, q4x0b, float);                   //!  Q-vector x component, with event eta gap B (harmonic 4 and power 1)
DECLARE_SOA_COLUMN(Q4Y0B, q4y0b, float);                   //!  Q-vector y component, with event eta gap B (harmonic 4 and power 1)
DECLARE_SOA_COLUMN(Q4X0C, q4x0c, float);                   //!  Q-vector x component, with event eta gap C (harmonic 4 and power 1)
DECLARE_SOA_COLUMN(Q4Y0C, q4y0c, float);                   //!  Q-vector y component, with event eta gap C (harmonic 4 and power 1)
DECLARE_SOA_COLUMN(Q42XA, q42xa, float);                   //!  Q-vector x component, with event eta gap A (harmonic 4 and power 2)
DECLARE_SOA_COLUMN(Q42YA, q42ya, float);                   //!  Q-vector y component, with event eta gap A (harmonic 4 and power 2)
DECLARE_SOA_COLUMN(Q23XA, q23xa, float);                   //!  Q-vector x component, with event eta gap A (harmonic 2 and power 3)
DECLARE_SOA_COLUMN(Q23YA, q23ya, float);                   //!  Q-vector y component, with event eta gap A (harmonic 2 and power 3)
DECLARE_SOA_COLUMN(S11A, s11a, float);                     //! Weighted multiplicity (p = 1, k = 1)
DECLARE_SOA_COLUMN(S12A, s12a, float);                     //! Weighted multiplicity (p = 1, k = 2)
DECLARE_SOA_COLUMN(S13A, s13a, float);                     //! Weighted multiplicity (p = 1, k = 3)
DECLARE_SOA_COLUMN(S31A, s31a, float);                     //! Weighted multiplicity (p = 3, k = 1)
DECLARE_SOA_COLUMN(CORR2REF, corr2ref, float);             //!  Ref Flow correlator <2>
DECLARE_SOA_COLUMN(CORR2REFetagap, corr2refetagap, float); //!  Ref Flow correlator <2>
DECLARE_SOA_COLUMN(CORR4REF, corr4ref, float);             //!  Ref Flow correlator <4>
DECLARE_SOA_COLUMN(M11REF, m11ref, float);                 //!  Weighted multiplicity of <<2>> for reference flow
DECLARE_SOA_COLUMN(M1111REF, m1111ref, float);             //!  Weighted multiplicity of <<4>> for reference flow
DECLARE_SOA_COLUMN(M11REFetagap, m11refetagap, float);     //!  Weighted multiplicity of <<2>>  etagap for reference flow
} // namespace reducedevent

DECLARE_SOA_TABLE_STAGED(ReducedEvents, "REDUCEDEVENT", //!   Main event information table
                         o2::soa::Index<>,
                         reducedevent::Tag, bc::RunNumber,
                         collision::PosX, collision::PosY, collision::PosZ, collision::NumContrib,
                         collision::CollisionTime, collision::CollisionTimeRes);

DECLARE_SOA_TABLE(ReducedEventsExtended_000, "AOD", "REEXTENDED", //!  Extended event information
                  bc::GlobalBC, evsel::Alias, evsel::Selection, timestamp::Timestamp, cent::CentRun2V0M,
                  mult::MultTPC, mult::MultFV0A, mult::MultFV0C, mult::MultFT0A, mult::MultFT0C,
                  mult::MultFDDA, mult::MultFDDC, mult::MultZNA, mult::MultZNC, mult::MultTracklets, mult::MultNTracksPV,
                  cent::CentFT0C);

DECLARE_SOA_TABLE_VERSIONED(ReducedEventsExtended_001, "AOD", "REEXTENDED", 1, //!  Extended event information
                            bc::GlobalBC, evsel::Alias, evsel::Selection, timestamp::Timestamp, cent::CentRun2V0M,
                            mult::MultTPC, mult::MultFV0A, mult::MultFV0C, mult::MultFT0A, mult::MultFT0C,
                            mult::MultFDDA, mult::MultFDDC, mult::MultZNA, mult::MultZNC, mult::MultTracklets, mult::MultNTracksPV,
                            cent::CentFT0C, cent::CentFT0A, cent::CentFT0M);

using ReducedEventsExtended = ReducedEventsExtended_001;

DECLARE_SOA_TABLE(ReducedEventsMultPV_000, "AOD", "REMULTPV", //!  Multiplicity information for primary vertex
                  mult::MultNTracksHasITS, mult::MultNTracksHasTPC, mult::MultNTracksHasTOF, mult::MultNTracksHasTRD,
                  mult::MultNTracksITSOnly, mult::MultNTracksTPCOnly, mult::MultNTracksITSTPC,
                  evsel::NumTracksInTimeRange);

DECLARE_SOA_TABLE_VERSIONED(ReducedEventsMultPV_001, "AOD", "REMULTPV", 1, //!  Multiplicity information for primary vertex
                            mult::MultNTracksHasITS, mult::MultNTracksHasTPC, mult::MultNTracksHasTOF, mult::MultNTracksHasTRD,
                            mult::MultNTracksITSOnly, mult::MultNTracksTPCOnly, mult::MultNTracksITSTPC,
                            mult::MultNTracksPVeta1, mult::MultNTracksPVetaHalf, evsel::NumTracksInTimeRange, evsel::SumAmpFT0CInTimeRange);

using ReducedEventsMultPV = ReducedEventsMultPV_001;

DECLARE_SOA_TABLE(ReducedEventsMultAll, "AOD", "REMULTALL", //!  Multiplicity information for all tracks in the event
                  mult::MultAllTracksTPCOnly, mult::MultAllTracksITSTPC,
                  reducedevent::NTPCoccupContribLongA, reducedevent::NTPCoccupContribLongC,
                  reducedevent::NTPCoccupMeanTimeLongA, reducedevent::NTPCoccupMeanTimeLongC,
                  reducedevent::NTPCoccupMedianTimeLongA, reducedevent::NTPCoccupMedianTimeLongC,
                  reducedevent::NTPCoccupContribShortA, reducedevent::NTPCoccupContribShortC,
                  reducedevent::NTPCoccupMeanTimeShortA, reducedevent::NTPCoccupMeanTimeShortC,
                  reducedevent::NTPCoccupMedianTimeShortA, reducedevent::NTPCoccupMedianTimeShortC);

DECLARE_SOA_TABLE(ReducedEventsVtxCov, "AOD", "REVTXCOV", //!    Event vertex covariance matrix
                  collision::CovXX, collision::CovXY, collision::CovXZ,
                  collision::CovYY, collision::CovYZ, collision::CovZZ, collision::Chi2);

DECLARE_SOA_TABLE(ReducedEventsQvector, "AOD", "REQVECTOR", //!    Event Q-vector information
                  reducedevent::Q1X0A, reducedevent::Q1Y0A, reducedevent::Q1X0B, reducedevent::Q1Y0B, reducedevent::Q1X0C, reducedevent::Q1Y0C,
                  reducedevent::Q2X0A, reducedevent::Q2Y0A, reducedevent::Q2X0B, reducedevent::Q2Y0B, reducedevent::Q2X0C, reducedevent::Q2Y0C,
                  reducedevent::MultA, reducedevent::MultB, reducedevent::MultC,
                  reducedevent::Q3X0A, reducedevent::Q3Y0A, reducedevent::Q3X0B, reducedevent::Q3Y0B, reducedevent::Q3X0C, reducedevent::Q3Y0C,
                  reducedevent::Q4X0A, reducedevent::Q4Y0A, reducedevent::Q4X0B, reducedevent::Q4Y0B, reducedevent::Q4X0C, reducedevent::Q4Y0C);

DECLARE_SOA_TABLE(ReducedEventsQvectorExtra, "AOD", "REQVECTOREXTRA", //!    Event Q-vector extra information
                  reducedevent::Q42XA, reducedevent::Q42YA, reducedevent::Q23XA, reducedevent::Q23YA,
                  reducedevent::S11A, reducedevent::S12A, reducedevent::S13A, reducedevent::S31A);

DECLARE_SOA_TABLE(ReducedEventsQvectorCentr, "AOD", "REQVECTORCTR", //!    Event Q-vector information from central framework
                  qvec::QvecFT0ARe, qvec::QvecFT0AIm, qvec::QvecFT0CRe, qvec::QvecFT0CIm, qvec::QvecFT0MRe, qvec::QvecFT0MIm, qvec::QvecFV0ARe, qvec::QvecFV0AIm, reducedevent::QvecBPosRe, reducedevent::QvecBPosIm, reducedevent::QvecBNegRe, reducedevent::QvecBNegIm,
                  qvec::SumAmplFT0A, qvec::SumAmplFT0C, qvec::SumAmplFT0M, qvec::SumAmplFV0A, reducedevent::NTrkBPos, reducedevent::NTrkBNeg);

DECLARE_SOA_TABLE(ReducedEventsQvectorCentrExtra, "AOD", "REQVECCTREXTA", //!    Event Q-vector information from central framework with TPC all
                  reducedevent::QvecBAllRe, reducedevent::QvecBAllIm, reducedevent::NTrkBAll);

DECLARE_SOA_TABLE(ReducedEventsRefFlow, "AOD", "REREFFLOW", //!    Event Ref Flow information
                  reducedevent::M11REF, reducedevent::M11REFetagap, reducedevent::M1111REF, reducedevent::CORR2REF, reducedevent::CORR2REFetagap, reducedevent::CORR4REF, cent::CentFT0C);

DECLARE_SOA_TABLE(ReducedEventsQvectorZN, "AOD", "REQVECTORZN", //!    Event Q-vector information from ZNs detectors
                  reducedevent::Q1ZNAX, reducedevent::Q1ZNAY, reducedevent::Q1ZNCX, reducedevent::Q1ZNCY);

DECLARE_SOA_TABLE(ReducedEventsInfo, "AOD", "REDUCEVENTINFO", //!   Main event index table
                  reducedevent::CollisionId);

DECLARE_SOA_TABLE(ReducedEventsMergingTable, "AOD", "REMERGE", //!   Collision merging quatities
                  reducedevent::DCAzBimodalityCoefficient, reducedevent::DCAzBimodalityCoefficientBinned,
                  reducedevent::DCAzBimodalityCoefficientBinnedTrimmed1, reducedevent::DCAzBimodalityCoefficientBinnedTrimmed2, reducedevent::DCAzBimodalityCoefficientBinnedTrimmed3,
                  reducedevent::DCAzMean, reducedevent::DCAzMeanBinnedTrimmed1, reducedevent::DCAzMeanBinnedTrimmed2, reducedevent::DCAzMeanBinnedTrimmed3,
                  reducedevent::DCAzRMS, reducedevent::DCAzRMSBinnedTrimmed1, reducedevent::DCAzRMSBinnedTrimmed2, reducedevent::DCAzRMSBinnedTrimmed3,
                  reducedevent::DCAzSkewness, reducedevent::DCAzKurtosis,
                  reducedevent::DCAzFracAbove100um, reducedevent::DCAzFracAbove200um, reducedevent::DCAzFracAbove500um,
                  reducedevent::DCAzFracAbove1mm, reducedevent::DCAzFracAbove2mm, reducedevent::DCAzFracAbove5mm, reducedevent::DCAzFracAbove10mm,
                  reducedevent::DCAzNPeaks, reducedevent::DCAzNPeaksTrimmed1, reducedevent::DCAzNPeaksTrimmed2, reducedevent::DCAzNPeaksTrimmed3);

// TODO and NOTE: This table is just an extension of the ReducedEvents table
//       There is no explicit accounting for MC events which were not reconstructed!!!
//       However, for analysis which will require these events, a special skimming process function
//           can be constructed and the same data model could be used

DECLARE_SOA_TABLE(ReducedMCEvents_000, "AOD", "REDUCEDMCEVENT", //!   Event level MC truth information
                  o2::soa::Index<>,
                  mccollision::GeneratorsID, reducedevent::MCPosX, reducedevent::MCPosY, reducedevent::MCPosZ,
                  mccollision::T, mccollision::Weight, mccollision::ImpactParameter);

DECLARE_SOA_TABLE_VERSIONED(ReducedMCEvents_001, "AOD", "REDUCEDMCEVENT", 1, //!   Event level MC truth information
                            o2::soa::Index<>,
                            mccollision::GeneratorsID, reducedevent::MCPosX, reducedevent::MCPosY, reducedevent::MCPosZ,
                            mccollision::T, mccollision::Weight, mccollision::ImpactParameter, cent::CentFT0C,
                            mult::MultMCNParticlesEta05, mult::MultMCNParticlesEta08, mult::MultMCNParticlesEta10);

using ReducedMCEvents = ReducedMCEvents_001;

using ReducedEvent = ReducedEvents::iterator;
using StoredReducedEvent = StoredReducedEvents::iterator;
using ReducedEventExtended = ReducedEventsExtended::iterator;
using ReducedEventVtxCov = ReducedEventsVtxCov::iterator;
using ReducedEventMultPV = ReducedEventsMultPV::iterator;
using ReducedEventMultAll = ReducedEventsMultAll::iterator;
using ReducedEventQvector = ReducedEventsQvector::iterator;
using ReducedEventQvectorExtra = ReducedEventsQvectorExtra::iterator;
using ReducedEventQvectorCentr = ReducedEventsQvectorCentr::iterator;
using ReducedEventQvectorCentrExtra = ReducedEventsQvectorCentrExtra::iterator;
using ReducedEventRefFlow = ReducedEventsRefFlow::iterator;
using ReducedEventQvectorZN = ReducedEventsQvectorZN::iterator;
using ReducedMCEvent = ReducedMCEvents::iterator;

namespace reducedeventlabel
{
DECLARE_SOA_INDEX_COLUMN(ReducedMCEvent, reducedMCevent); //! MC collision
DECLARE_SOA_COLUMN(McMask, mcMask, uint16_t);             //! Bit mask to indicate collision mismatches (bit ON means mismatch). Bit 15: indicates negative label
} // namespace reducedeventlabel

DECLARE_SOA_TABLE(ReducedMCEventLabels, "AOD", "REMCCOLLBL", //! Table joined to the ReducedEvents table containing the MC index
                  reducedeventlabel::ReducedMCEventId, reducedeventlabel::McMask);
using ReducedMCEventLabel = ReducedMCEventLabels::iterator;

namespace reducedzdc
{
DECLARE_SOA_COLUMN(EnergyCommonZNA, energyCommonZNA, float); //!
DECLARE_SOA_COLUMN(EnergyCommonZNC, energyCommonZNC, float); //!
DECLARE_SOA_COLUMN(EnergyCommonZPA, energyCommonZPA, float); //!
DECLARE_SOA_COLUMN(EnergyCommonZPC, energyCommonZPC, float); //!
DECLARE_SOA_COLUMN(EnergyZNA1, energyZNA1, float);           //!
DECLARE_SOA_COLUMN(EnergyZNA2, energyZNA2, float);           //!
DECLARE_SOA_COLUMN(EnergyZNA3, energyZNA3, float);           //!
DECLARE_SOA_COLUMN(EnergyZNA4, energyZNA4, float);           //!
DECLARE_SOA_COLUMN(EnergyZNC1, energyZNC1, float);           //!
DECLARE_SOA_COLUMN(EnergyZNC2, energyZNC2, float);           //!
DECLARE_SOA_COLUMN(EnergyZNC3, energyZNC3, float);           //!
DECLARE_SOA_COLUMN(EnergyZNC4, energyZNC4, float);           //!
DECLARE_SOA_COLUMN(TimeZNA, timeZNA, float);                 //!
DECLARE_SOA_COLUMN(TimeZNC, timeZNC, float);                 //!
DECLARE_SOA_COLUMN(TimeZPA, timeZPA, float);                 //!
DECLARE_SOA_COLUMN(TimeZPC, timeZPC, float);                 //!
} // namespace reducedzdc

DECLARE_SOA_TABLE(ReducedZdcs, "AOD", "REDUCEDZDC", //!   Event ZDC information
                  reducedzdc::EnergyCommonZNA, reducedzdc::EnergyCommonZNC,
                  reducedzdc::EnergyCommonZPA, reducedzdc::EnergyCommonZPC,
                  reducedzdc::TimeZNA, reducedzdc::TimeZNC,
                  reducedzdc::TimeZPA, reducedzdc::TimeZPC);

DECLARE_SOA_TABLE(ReducedZdcsExtra, "AOD", "REDUCEDZDCEXTRA", //!   Event ZDC extra information
                  reducedzdc::EnergyZNA1, reducedzdc::EnergyZNA2, reducedzdc::EnergyZNA3, reducedzdc::EnergyZNA4,
                  reducedzdc::EnergyZNC1, reducedzdc::EnergyZNC2, reducedzdc::EnergyZNC3, reducedzdc::EnergyZNC4);

using ReducedZdc = ReducedZdcs::iterator;
using ReducedZdcExtra = ReducedZdcsExtra::iterator;

namespace reducedfit
{
// FIT detector information (based on upchelpers::FITInfo structure)
DECLARE_SOA_COLUMN(AmplitudeFT0A, amplitudeFT0A, float);         //! FT0A total amplitude
DECLARE_SOA_COLUMN(AmplitudeFT0C, amplitudeFT0C, float);         //! FT0C total amplitude
DECLARE_SOA_COLUMN(TimeFT0A, timeFT0A, float);                   //! FT0A time
DECLARE_SOA_COLUMN(TimeFT0C, timeFT0C, float);                   //! FT0C time
DECLARE_SOA_COLUMN(TriggerMaskFT0, triggerMaskFT0, uint8_t);     //! FT0 trigger mask
DECLARE_SOA_COLUMN(NFiredChannelsFT0A, nFiredChannelsFT0A, int); //! Number of fired channels in FT0A
DECLARE_SOA_COLUMN(NFiredChannelsFT0C, nFiredChannelsFT0C, int); //! Number of fired channels in FT0C
DECLARE_SOA_COLUMN(AmplitudeFDDA, amplitudeFDDA, float);         //! FDDA total amplitude
DECLARE_SOA_COLUMN(AmplitudeFDDC, amplitudeFDDC, float);         //! FDDC total amplitude
DECLARE_SOA_COLUMN(TimeFDDA, timeFDDA, float);                   //! FDDA time
DECLARE_SOA_COLUMN(TimeFDDC, timeFDDC, float);                   //! FDDC time
DECLARE_SOA_COLUMN(TriggerMaskFDD, triggerMaskFDD, uint8_t);     //! FDD trigger mask
DECLARE_SOA_COLUMN(AmplitudeFV0A, amplitudeFV0A, float);         //! FV0A total amplitude
DECLARE_SOA_COLUMN(TimeFV0A, timeFV0A, float);                   //! FV0A time
DECLARE_SOA_COLUMN(TriggerMaskFV0A, triggerMaskFV0A, uint8_t);   //! FV0A trigger mask
DECLARE_SOA_COLUMN(NFiredChannelsFV0A, nFiredChannelsFV0A, int); //! Number of fired channels in FV0A
DECLARE_SOA_COLUMN(BBFT0Apf, bbFT0Apf, int32_t);                 //! Beam-beam flags for FT0A
DECLARE_SOA_COLUMN(BGFT0Apf, bgFT0Apf, int32_t);                 //! Beam-gas flags for FT0A
DECLARE_SOA_COLUMN(BBFT0Cpf, bbFT0Cpf, int32_t);                 //! Beam-beam flags for FT0C
DECLARE_SOA_COLUMN(BGFT0Cpf, bgFT0Cpf, int32_t);                 //! Beam-gas flags for FT0C
DECLARE_SOA_COLUMN(BBFV0Apf, bbFV0Apf, int32_t);                 //! Beam-beam flags for FV0A
DECLARE_SOA_COLUMN(BGFV0Apf, bgFV0Apf, int32_t);                 //! Beam-gas flags for FV0A
DECLARE_SOA_COLUMN(BBFDDApf, bbFDDApf, int32_t);                 //! Beam-beam flags for FDDA
DECLARE_SOA_COLUMN(BGFDDApf, bgFDDApf, int32_t);                 //! Beam-gas flags for FDDA
DECLARE_SOA_COLUMN(BBFDDCpf, bbFDDCpf, int32_t);                 //! Beam-beam flags for FDDC
DECLARE_SOA_COLUMN(BGFDDCpf, bgFDDCpf, int32_t);                 //! Beam-gas flags for FDDC
} // namespace reducedfit

DECLARE_SOA_TABLE(ReducedFITs, "AOD", "REDUCEDFIT", //! FIT detector information
                  reducedfit::AmplitudeFT0A, reducedfit::AmplitudeFT0C,
                  reducedfit::TimeFT0A, reducedfit::TimeFT0C,
                  reducedfit::TriggerMaskFT0,
                  reducedfit::NFiredChannelsFT0A, reducedfit::NFiredChannelsFT0C,
                  reducedfit::AmplitudeFDDA, reducedfit::AmplitudeFDDC,
                  reducedfit::TimeFDDA, reducedfit::TimeFDDC,
                  reducedfit::TriggerMaskFDD,
                  reducedfit::AmplitudeFV0A, reducedfit::TimeFV0A,
                  reducedfit::TriggerMaskFV0A,
                  reducedfit::NFiredChannelsFV0A,
                  reducedfit::BBFT0Apf, reducedfit::BGFT0Apf,
                  reducedfit::BBFT0Cpf, reducedfit::BGFT0Cpf,
                  reducedfit::BBFV0Apf, reducedfit::BGFV0Apf,
                  reducedfit::BBFDDApf, reducedfit::BGFDDApf,
                  reducedfit::BBFDDCpf, reducedfit::BGFDDCpf);

using ReducedFIT = ReducedFITs::iterator;

namespace reducedtrack
{
// basic track information
DECLARE_SOA_INDEX_COLUMN(ReducedEvent, reducedevent); //!
DECLARE_SOA_INDEX_COLUMN(Track, track);               //!
// ----  flags reserved for storing various information during filtering
DECLARE_SOA_BITMAP_COLUMN(FilteringFlags, filteringFlags, 64); //!
// -----------------------------------------------------
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
} // namespace reducedtrack

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
                  track::X, track::Alpha, track::IsWithinBeamPipe<track::X>,
                  track::Y, track::Z, track::Snp, track::Tgl, track::Signed1Pt,
                  track::TPCInnerParam, track::Flags, // tracking status flags
                  track::ITSClusterMap, track::ITSChi2NCl,
                  track::TPCNClsFindable, track::TPCNClsFindableMinusFound, track::TPCNClsFindableMinusCrossedRows,
                  track::TPCNClsShared, track::TPCChi2NCl,
                  track::TRDChi2, track::TRDPattern, track::TOFChi2, track::Length, reducedtrack::DcaXY, reducedtrack::DcaZ,
                  track::TrackTime, track::TrackTimeRes, track::TOFExpMom,
                  reducedtrack::DetectorMap,
                  track::IsPVContributor<track::Flags>,
                  track::TPCNClsFound<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  track::TPCNClsCrossedRows<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                  reducedtrack::HasITS<reducedtrack::DetectorMap>, reducedtrack::HasTRD<reducedtrack::DetectorMap>,
                  reducedtrack::HasTOF<reducedtrack::DetectorMap>, reducedtrack::HasTPC<reducedtrack::DetectorMap>);

// barrel covariance matrix  TODO: add all the elements required for secondary vertexing
DECLARE_SOA_TABLE(ReducedTracksBarrelCov, "AOD", "RTBARRELCOV", //!
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

// barrel collision information (joined with ReducedTracks) allowing to connect different tables (cross PWGs)
DECLARE_SOA_TABLE(ReducedTracksBarrelInfo, "AOD", "RTBARRELINFO",
                  reducedtrack::CollisionId, collision::PosX, collision::PosY, collision::PosZ, reducedtrack::TrackId);

using ReducedTrack = ReducedTracks::iterator;
using ReducedTrackBarrel = ReducedTracksBarrel::iterator;
using ReducedTrackBarrelCov = ReducedTracksBarrelCov::iterator;
using ReducedTrackBarrelPID = ReducedTracksBarrelPID::iterator;
using ReducedTrackBarrelInfo = ReducedTracksBarrelInfo::iterator;

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
DECLARE_SOA_TABLE(ReducedMCTracks, "AOD", "REDUCEDMCTRACK", //!  MC track information (on disk)
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
                  mcparticle::GetHepMCStatusCode<mcparticle::Flags, mcparticle::StatusCode>,
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

// MFT track quantities
namespace reducedmft
{
DECLARE_SOA_INDEX_COLUMN(ReducedEvent, reducedevent);        //!
DECLARE_SOA_COLUMN(FilteringFlags, filteringFlags, uint8_t); //!

DECLARE_SOA_COLUMN(Pt, pt, float);                                                        //!
DECLARE_SOA_COLUMN(Eta, eta, float);                                                      //!
DECLARE_SOA_COLUMN(Phi, phi, float);                                                      //!
DECLARE_SOA_COLUMN(Sign, sign, int);                                                      //!
DECLARE_SOA_COLUMN(FwdDcaX, fwdDcaX, float);                                              //!
DECLARE_SOA_COLUMN(FwdDcaY, fwdDcaY, float);                                              //!
DECLARE_SOA_COLUMN(MftClusterSizesAndTrackFlags, mftClusterSizesAndTrackFlags, uint64_t); //!
DECLARE_SOA_COLUMN(MftNClusters, mftNClusters, int);                                      //!
DECLARE_SOA_INDEX_COLUMN(ReducedMCTrack, reducedMCTrack);                                 //!
DECLARE_SOA_COLUMN(McMask, mcMask, uint16_t);                                             //!
DECLARE_SOA_COLUMN(McReducedFlags, mcReducedFlags, uint16_t);                             //!
} // namespace reducedmft

// MFT track kinematics
DECLARE_SOA_TABLE(ReducedMFTs, "AOD", "REDUCEDMFT", //!
                  o2::soa::Index<>, reducedmft::ReducedEventId, reducedmft::FilteringFlags,
                  reducedmft::Pt, reducedmft::Eta, reducedmft::Phi);

// MFT tracks extra info (cluster size, sign)
DECLARE_SOA_TABLE(ReducedMFTsExtra, "AOD", "RMFTEXTRA", //!
                  reducedmft::MftClusterSizesAndTrackFlags, reducedmft::Sign,
                  reducedmft::FwdDcaX, reducedmft::FwdDcaY, reducedmft::MftNClusters);

DECLARE_SOA_TABLE(ReducedMFTLabels, "AOD", "RTMFTLABELS", //!
                  reducedmft::ReducedMCTrackId, reducedmft::McMask, reducedmft::McReducedFlags);

// iterator
using ReducedMFT = ReducedMFTs::iterator;
using ReducedMFTExtra = ReducedMFTsExtra::iterator;
using ReducedMFTLabel = ReducedMFTLabels::iterator;

// muon quantities
namespace reducedmuon
{
DECLARE_SOA_INDEX_COLUMN(ReducedEvent, reducedevent);        //!
DECLARE_SOA_COLUMN(FilteringFlags, filteringFlags, uint8_t); //!
// the (pt,eta,phi,sign) will be computed in the skimming task //!
DECLARE_SOA_COLUMN(Pt, pt, float);                 //!
DECLARE_SOA_COLUMN(Eta, eta, float);               //!
DECLARE_SOA_COLUMN(Phi, phi, float);               //!
DECLARE_SOA_COLUMN(Sign, sign, int);               //!
DECLARE_SOA_COLUMN(FwdDcaX, fwdDcaX, float);       //!  Impact parameter in X of forward track to the primary vertex
DECLARE_SOA_COLUMN(FwdDcaY, fwdDcaY, float);       //!  Impact parameter in Y of forward track to the primary vertex
DECLARE_SOA_COLUMN(IsAmbiguous, isAmbiguous, int); //!
DECLARE_SOA_INDEX_COLUMN(Collision, collision);    //!
DECLARE_SOA_DYNAMIC_COLUMN(Px, px,                 //!
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
DECLARE_SOA_SELF_INDEX_COLUMN_FULL(MatchMCHTrack, matchMCHTrack, int, "ReducedMuons_MatchMCHTrack");
DECLARE_SOA_INDEX_COLUMN(ReducedMFT, matchMFTTrack); //!  matching index pointing to the ReducedMFTTrack table if filled
} // namespace reducedmuon

// Muon track kinematics
DECLARE_SOA_TABLE(ReducedMuons, "AOD", "REDUCEDMUON", //!
                  o2::soa::Index<>, reducedmuon::ReducedEventId,
                  reducedmuon::MatchMCHTrackId, reducedmuon::ReducedMFTId,
                  reducedmuon::FilteringFlags,
                  reducedmuon::Pt, reducedmuon::Eta, reducedmuon::Phi, reducedmuon::Sign, reducedmuon::IsAmbiguous,
                  reducedmuon::Px<reducedmuon::Pt, reducedmuon::Phi>,
                  reducedmuon::Py<reducedmuon::Pt, reducedmuon::Phi>,
                  reducedmuon::Pz<reducedmuon::Pt, reducedmuon::Eta>,
                  reducedmuon::P<reducedmuon::Pt, reducedmuon::Eta>);

// Muon track quality details
DECLARE_SOA_TABLE(ReducedMuonsExtra, "AOD", "RTMUONEXTRA", //!
                  fwdtrack::NClusters, fwdtrack::PDca, fwdtrack::RAtAbsorberEnd,
                  fwdtrack::Chi2, fwdtrack::Chi2MatchMCHMID, fwdtrack::Chi2MatchMCHMFT,
                  fwdtrack::MatchScoreMCHMFT,
                  fwdtrack::MCHBitMap, fwdtrack::MIDBitMap, fwdtrack::MIDBoards, fwdtrack::TrackType,
                  reducedmuon::FwdDcaX, reducedmuon::FwdDcaY,
                  fwdtrack::TrackTime, fwdtrack::TrackTimeRes);

// Muon covariance, TODO: the rest of the matrix should be added when needed
DECLARE_SOA_TABLE(ReducedMuonsCov, "AOD", "RTMUONCOV",
                  fwdtrack::X, fwdtrack::Y, fwdtrack::Z, reducedmuon::RawPhi, fwdtrack::Tgl, fwdtrack::Signed1Pt,
                  fwdtrack::CXX, fwdtrack::CXY, fwdtrack::CYY, fwdtrack::CPhiX, fwdtrack::CPhiY, fwdtrack::CPhiPhi,
                  fwdtrack::CTglX, fwdtrack::CTglY, fwdtrack::CTglPhi, fwdtrack::CTglTgl, fwdtrack::C1PtX,
                  fwdtrack::C1PtY, fwdtrack::C1PtPhi, fwdtrack::C1PtTgl, fwdtrack::C1Pt21Pt2);

// Muon collision information (joined with ReducedMuons) allowing to connect different tables (cross PWGs)
DECLARE_SOA_TABLE(ReducedMuonsInfo, "AOD", "RTMUONINFO",
                  reducedmuon::CollisionId, collision::PosX, collision::PosY, collision::PosZ);

// iterators
using ReducedMuon = ReducedMuons::iterator;
using ReducedMuonExtra = ReducedMuonsExtra::iterator;
using ReducedMuonCov = ReducedMuonsCov::iterator;
using ReducedMuonInfo = ReducedMuonsInfo::iterator;

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

namespace reducedtrack_association
{
DECLARE_SOA_INDEX_COLUMN(ReducedEvent, reducedevent); //! ReducedEvent index
DECLARE_SOA_INDEX_COLUMN(ReducedTrack, reducedtrack); //! ReducedTrack index
DECLARE_SOA_INDEX_COLUMN(ReducedMuon, reducedmuon);   //! ReducedMuon index
DECLARE_SOA_INDEX_COLUMN(ReducedMFT, reducedmft);     //! ReducedMFTTrack index
} // namespace reducedtrack_association

DECLARE_SOA_TABLE(ReducedTracksAssoc, "AOD", "RTASSOC", //! Table for reducedtrack-to-reducedcollision association
                  reducedtrack_association::ReducedEventId,
                  reducedtrack_association::ReducedTrackId);
DECLARE_SOA_TABLE(ReducedMuonsAssoc, "AOD", "RMASSOC", //! Table for reducedmuon-to-reducedcollision association
                  reducedtrack_association::ReducedEventId,
                  reducedtrack_association::ReducedMuonId);
DECLARE_SOA_TABLE(ReducedMFTAssoc, "AOD", "RMFTASSOC", //! Table for reducemft-to-reducedcollision association
                  reducedtrack_association::ReducedEventId,
                  reducedtrack_association::ReducedMFTId);

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
DECLARE_SOA_COLUMN(Chi21, chi21, float);                       //! Chi2 for Muon Tracks

DECLARE_SOA_COLUMN(Chi2MatchMCHMID2, chi2MatchMCHMID2, float); //! MCH-MID Match Chi2 for MUONStandalone tracks
DECLARE_SOA_COLUMN(Chi2MatchMCHMFT2, chi2MatchMCHMFT2, float); //! MCH-MFT Match Chi2 for GlobalMuonTracks
DECLARE_SOA_COLUMN(Chi22, chi22, float);                       //! Chi2 for Muon Tracks

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

DECLARE_SOA_COLUMN(IsAmbig1, isAmbig1, int); //!
DECLARE_SOA_COLUMN(IsAmbig2, isAmbig2, int); //!
DECLARE_SOA_COLUMN(IsCorrectAssoc1, isCorrectAssoc1, bool); //!
DECLARE_SOA_COLUMN(IsCorrectAssoc2, isCorrectAssoc2, bool); //!

DECLARE_SOA_COLUMN(FwdDcaX1, fwdDcaX1, float);               //! X component of forward DCA
DECLARE_SOA_COLUMN(FwdDcaY1, fwdDcaY1, float);               //! Y component of forward DCA
DECLARE_SOA_COLUMN(FwdDcaX2, fwdDcaX2, float);               //! X component of forward DCA
DECLARE_SOA_COLUMN(FwdDcaY2, fwdDcaY2, float);               //! Y component of forward DCA
DECLARE_SOA_COLUMN(ITSNCls1, itsNCls1, int);                 //! Number of ITS clusters
DECLARE_SOA_COLUMN(ITSClusterMap1, itsClusterMap1, uint8_t); //! ITS clusters map
DECLARE_SOA_COLUMN(ITSChi2NCl1, itsChi2NCl1, float);         //! ITS chi2/Ncls
DECLARE_SOA_COLUMN(TPCNClsFound1, tpcNClsFound1, float);     //! Number of TPC clusters found
DECLARE_SOA_COLUMN(TPCNClsCR1, tpcNClsCR1, float);           //! Number of TPC crossed rows
DECLARE_SOA_COLUMN(TPCChi2NCl1, tpcChi2NCl1, float);         //! TPC chi2/Ncls
DECLARE_SOA_COLUMN(DcaXY1, dcaXY1, float);                   //! DCA in XY plane
DECLARE_SOA_COLUMN(DcaZ1, dcaZ1, float);                     //! DCA in Z
DECLARE_SOA_COLUMN(TPCSignal1, tpcSignal1, float);           //! TPC dE/dx signal
DECLARE_SOA_COLUMN(TPCNSigmaEl1, tpcNSigmaEl1, float);       //! TPC nSigma electron
DECLARE_SOA_COLUMN(TPCNSigmaPi1, tpcNSigmaPi1, float);       //! TPC nSigma pion
DECLARE_SOA_COLUMN(TPCNSigmaPr1, tpcNSigmaPr1, float);       //! TPC nSigma proton
DECLARE_SOA_COLUMN(TOFBeta1, tofBeta1, float);               //! TOF beta
DECLARE_SOA_COLUMN(TOFNSigmaEl1, tofNSigmaEl1, float);       //! TOF nSigma electron
DECLARE_SOA_COLUMN(TOFNSigmaPi1, tofNSigmaPi1, float);       //! TOF nSigma pion
DECLARE_SOA_COLUMN(TOFNSigmaPr1, tofNSigmaPr1, float);       //! TOF nSigma proton
DECLARE_SOA_COLUMN(ITSNCls2, itsNCls2, int);                 //! Number of ITS clusters
DECLARE_SOA_COLUMN(ITSClusterMap2, itsClusterMap2, uint8_t); //! ITS clusters map
DECLARE_SOA_COLUMN(ITSChi2NCl2, itsChi2NCl2, float);         //! ITS chi2/Ncls
DECLARE_SOA_COLUMN(TPCNClsFound2, tpcNClsFound2, float);     //! Number of TPC clusters found
DECLARE_SOA_COLUMN(TPCNClsCR2, tpcNClsCR2, float);           //! Number of TPC crossed rows
DECLARE_SOA_COLUMN(TPCChi2NCl2, tpcChi2NCl2, float);         //! TPC chi2/Ncls
DECLARE_SOA_COLUMN(DcaXY2, dcaXY2, float);                   //! DCA in XY plane
DECLARE_SOA_COLUMN(DcaZ2, dcaZ2, float);                     //! DCA in Z
DECLARE_SOA_COLUMN(TPCSignal2, tpcSignal2, float);           //! TPC dE/dx signal
DECLARE_SOA_COLUMN(TPCNSigmaEl2, tpcNSigmaEl2, float);       //! TPC nSigma electron
DECLARE_SOA_COLUMN(TPCNSigmaPi2, tpcNSigmaPi2, float);       //! TPC nSigma pion
DECLARE_SOA_COLUMN(TPCNSigmaPr2, tpcNSigmaPr2, float);       //! TPC nSigma proton
DECLARE_SOA_COLUMN(TOFBeta2, tofBeta2, float);               //! TOF beta
DECLARE_SOA_COLUMN(TOFNSigmaEl2, tofNSigmaEl2, float);       //! TOF nSigma electron
DECLARE_SOA_COLUMN(TOFNSigmaPi2, tofNSigmaPi2, float);       //! TOF nSigma pion
DECLARE_SOA_COLUMN(TOFNSigmaPr2, tofNSigmaPr2, float);       //! TOF nSigma proton

DECLARE_SOA_COLUMN(DCAxyzTrk0KF, dcaxyztrk0KF, float); //! 3D DCA to primary vertex of the first track
DECLARE_SOA_COLUMN(DCAxyzTrk1KF, dcaxyztrk1KF, float); //! 3D DCA to primary vertex of the second track
DECLARE_SOA_COLUMN(DCAxyTrk0KF, dcaxytrk0KF, float);   //! 2D DCA to primary vertex of the first track
DECLARE_SOA_COLUMN(DCAxyTrk1KF, dcaxytrk1KF, float);   //! 2D DCA to primary vertex of the second track

DECLARE_SOA_COLUMN(DeviationTrk0KF, deviationTrk0KF, float);     //! 3D chi2 deviation to primary vertex of the first track
DECLARE_SOA_COLUMN(DeviationTrk1KF, deviationTrk1KF, float);     //! 3D chi2 deviation to primary vertex of the second track
DECLARE_SOA_COLUMN(DeviationxyTrk0KF, deviationxyTrk0KF, float); //! 2D chi2 deviation to primary vertex of the first track
DECLARE_SOA_COLUMN(DeviationxyTrk1KF, deviationxyTrk1KF, float); //! 2D chi2 deviation to primary vertex of the second track
} // namespace dilepton_track_index

// pair information
namespace reducedpair
{
DECLARE_SOA_INDEX_COLUMN(ReducedEvent, reducedevent);
DECLARE_SOA_INDEX_COLUMN_FULL(Index0, index0, int, ReducedTracks, "_0");                 //! Index to first prong
DECLARE_SOA_INDEX_COLUMN_FULL(Index1, index1, int, ReducedTracks, "_1");                 //! Index to second prong
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, Tracks, "_0");                        //! Index of first prong in Tracks table
DECLARE_SOA_INDEX_COLUMN_FULL(Prong1, prong1, int, Tracks, "_1");                        //! Index of second prong in Tracks table
DECLARE_SOA_BITMAP_COLUMN(EventSelection, evSelection, 8);                               //! Event selection bits (ambiguity, splitting candidate)
DECLARE_SOA_COLUMN(Mass, mass, float);                                                   //!
DECLARE_SOA_COLUMN(Pt, pt, float);                                                       //!
DECLARE_SOA_COLUMN(Eta, eta, float);                                                     //!
DECLARE_SOA_COLUMN(Phi, phi, float);                                                     //!
DECLARE_SOA_COLUMN(Sign, sign, int);                                                     //!
DECLARE_SOA_BITMAP_COLUMN(FilterMap, filterMap, 32);                                     //!
DECLARE_SOA_BITMAP_COLUMN(PairFilterMap, pairFilterMap, 32);                             //!
DECLARE_SOA_BITMAP_COLUMN(CommonFilterMap, commonFilterMap, 32);                         //!
DECLARE_SOA_COLUMN(McDecision, mcDecision, uint32_t);                                    //!
DECLARE_SOA_COLUMN(Tauz, tauz, float);                                                   //! Longitudinal pseudo-proper time of lepton pair (in ns)
DECLARE_SOA_COLUMN(TauzErr, tauzErr, float);                                             //! Error on longitudinal pseudo-proper time of lepton pair (in ns)
DECLARE_SOA_COLUMN(VertexPz, vertexPz, float);                                           //! Longitudinal projection of impulsion
DECLARE_SOA_COLUMN(SVertex, sVertex, float);                                             //! Secondary vertex of lepton pair
DECLARE_SOA_COLUMN(Tauxy, tauxy, float);                                                 //! Transverse pseudo-proper time of lepton pair (in ns)
DECLARE_SOA_COLUMN(TauxyErr, tauxyErr, float);                                           //! Error on transverse pseudo-proper time of lepton pair (in ns)
DECLARE_SOA_COLUMN(Lz, lz, float);                                                       //! Longitudinal projection of decay length
DECLARE_SOA_COLUMN(Lxy, lxy, float);                                                     //! Transverse projection of decay length
DECLARE_SOA_COLUMN(Chi2pca, chi2pca, float);                                             //! Chi2 for PCA of the dilepton
DECLARE_SOA_COLUMN(CosPointingAngle, cosPointingAngle, float);                           //! Cosine of the pointing angle
DECLARE_SOA_COLUMN(U2Q2, u2q2, float);                                                   //! Scalar product between unitary vector with event flow vector (harmonic 2)
DECLARE_SOA_COLUMN(U3Q3, u3q3, float);                                                   //! Scalar product between unitary vector with event flow vector (harmonic 3)
DECLARE_SOA_COLUMN(Cos2DeltaPhi, cos2deltaphi, float);                                   //! Cosinus term using event plane angle (harmonic 2)
DECLARE_SOA_COLUMN(Cos3DeltaPhi, cos3deltaphi, float);                                   //! Cosinus term using event plane angle (harmonic 3)
DECLARE_SOA_COLUMN(R2SP_AB, r2spab, float);                                              //! Event plane resolution for SP method n=2 (A,B) TPC-FT0A
DECLARE_SOA_COLUMN(R2SP_AC, r2spac, float);                                              //! Event plane resolution for SP method n=2 (A,C) TPC-FT0C
DECLARE_SOA_COLUMN(R2SP_BC, r2spbc, float);                                              //! Event plane resolution for SP method n=2 (B,C) FT0A-FT0C
DECLARE_SOA_COLUMN(R3SP, r3sp, float);                                                   //! Event plane resolution for SP method n=3
DECLARE_SOA_COLUMN(R2EP, r2ep, float);                                                   //! Event plane resolution for EP method n=2
DECLARE_SOA_COLUMN(R2EP_AB, r2epab, float);                                              //! Event plane resolution for EP method n=2 (A,B) TPC-FT0A
DECLARE_SOA_COLUMN(R2EP_AC, r2epac, float);                                              //! Event plane resolution for EP method n=2 (A,C) TPC-FT0C
DECLARE_SOA_COLUMN(R2EP_BC, r2epbc, float);                                              //! Event plane resolution for EP method n=2 (B,C) FT0A-FT0C
DECLARE_SOA_COLUMN(R3EP, r3ep, float);                                                   //! Event plane resolution for EP method n=3
DECLARE_SOA_COLUMN(CORR2POI, corr2poi, float);                                           //! POI FLOW CORRELATOR <2'>
DECLARE_SOA_COLUMN(CORR4POI, corr4poi, float);                                           //! POI FLOW CORRELATOR <4'>
DECLARE_SOA_COLUMN(M01POI, m01poi, float);                                               //! POI event weight for <2'>
DECLARE_SOA_COLUMN(M0111POI, m0111poi, float);                                           //! POI event weight for <4'>
DECLARE_SOA_COLUMN(MultDimuons, multdimuons, int);                                       //! Dimuon multiplicity
DECLARE_SOA_COLUMN(CentFT0C, centft0c, float);                                           //! Centrality information from FT0C
DECLARE_SOA_COLUMN(CollisionId, collisionId, int32_t);                                   //!
DECLARE_SOA_COLUMN(IsFirst, isfirst, int);                                               //! Flag for the first dilepton in the collision
DECLARE_SOA_COLUMN(DCAxyzBetweenTrksKF, dcaxyzbetweentrksKF, float);                     //! DCAxyz between the two tracks
DECLARE_SOA_COLUMN(DCAxyBetweenTrksKF, dcaxybetweentrksKF, float);                       //! DCAxy between the two tracks
DECLARE_SOA_COLUMN(MassKFGeo, massKFGeo, float);                                         //! Pair mass from KFParticle
DECLARE_SOA_COLUMN(CosPAKFGeo, cosPAKFGeo, float);                                       //! Cosine of the pointing angle from KFParticle
DECLARE_SOA_COLUMN(Chi2OverNDFKFGeo, chi2overndfKFGeo, float);                           //! Chi2 over NDF from KFParticle
DECLARE_SOA_COLUMN(DecayLengthKFGeo, decaylengthKFGeo, float);                           //! Decay length from KFParticle
DECLARE_SOA_COLUMN(DecayLengthOverErrKFGeo, decaylengthovererrKFGeo, float);             //! Decay length over error from KFParticle
DECLARE_SOA_COLUMN(DecayLengthXYKFGeo, decaylengthxyKFGeo, float);                       //! Decay length XY from KFParticle
DECLARE_SOA_COLUMN(DecayLengthXYOverErrKFGeo, decaylengthxyovererrKFGeo, float);         //! Decay length XY over error from KFParticle
DECLARE_SOA_COLUMN(PseudoproperDecayTimeKFGeo, pseudoproperdecaytimeKFGeo, float);       //! Pseudoproper decay time from KFParticle
DECLARE_SOA_COLUMN(PseudoproperDecayTimeErrKFGeo, pseudoproperdecaytimeErrKFGeo, float); //! Pseudoproper decay time error from KFParticle
DECLARE_SOA_COLUMN(MassKFGeoTop, massKFGeoTop, float);                                   //! Pair mass after topological constraint from KFParticle
DECLARE_SOA_COLUMN(Chi2OverNDFKFGeoTop, chi2overndfKFGeoTop, float);                     //! Chi2 over NDF after topological constraint from KFParticle
DECLARE_SOA_COLUMN(PairDCAxyz, pairDCAxyz, float);                                       //! Pair DCAxyz to PV from KFParticle
DECLARE_SOA_COLUMN(PairDCAxy, pairDCAxy, float);                                         //! Pair DCAxy to PV from KFParticle
DECLARE_SOA_COLUMN(DeviationPairKF, deviationPairKF, float);                             //! Pair chi2 deviation to PV from KFParticle
DECLARE_SOA_COLUMN(DeviationxyPairKF, deviationxyPairKF, float);                         //! Pair chi2 deviation to PV in XY from KFParticle
// DECLARE_SOA_INDEX_COLUMN(ReducedMuon, reducedmuon2); //!
DECLARE_SOA_COLUMN(CosThetaHE, costhetaHE, float);             //! Cosine in the helicity frame
DECLARE_SOA_COLUMN(PhiHE, phiHe, float);                       //! Phi in the helicity frame
DECLARE_SOA_COLUMN(PhiTildeHE, phiTildeHe, float);             //! Tilde Phi in the helicity frame
DECLARE_SOA_COLUMN(CosThetaCS, costhetaCS, float);             //! Cosine in the Collins-Soper frame
DECLARE_SOA_COLUMN(PhiCS, phiCS, float);                       //! Phi in the Collins-Soper frame
DECLARE_SOA_COLUMN(PhiTildeCS, phiTildeCS, float);             //! Tilde Phi in the Collins-Soper frame
DECLARE_SOA_COLUMN(CosThetaPP, costhetaPP, float);             //! Cosine in the Production Plane frame
DECLARE_SOA_COLUMN(PhiPP, phiPP, float);                       //! Phi in the Production Plane frame
DECLARE_SOA_COLUMN(PhiTildePP, phiTildePP, float);             //! Tilde Phi in the Production Plane frame
DECLARE_SOA_COLUMN(CosThetaRM, costhetaRM, float);             //! Cosine in the Random frame
DECLARE_SOA_COLUMN(CosThetaStarTPC, costhetaStarTPC, float);   //! global polarization, event plane reconstructed from TPC tracks
DECLARE_SOA_COLUMN(CosThetaStarFT0A, costhetaStarFT0A, float); //! global polarization, event plane reconstructed from FT0A tracks
DECLARE_SOA_COLUMN(CosThetaStarFT0C, costhetaStarFT0C, float); //! global polarization, event plane reconstructed from FT0C tracks
DECLARE_SOA_DYNAMIC_COLUMN(Px, px,                             //!
                           [](float pt, float phi) -> float { return pt * std::cos(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py, //!
                           [](float pt, float phi) -> float { return pt * std::sin(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz, //!
                           [](float pt, float eta) -> float { return pt * std::sinh(eta); });
DECLARE_SOA_DYNAMIC_COLUMN(P, p, //!
                           [](float pt, float eta) -> float { return pt * std::cosh(eta); });
DECLARE_SOA_DYNAMIC_COLUMN(Rap, rap, //!
                           [](float pt, float eta, float m) -> float { return std::log((std::sqrt(m * m + pt * pt * std::cosh(eta) * std::cosh(eta)) + pt * std::sinh(eta)) / std::sqrt(m * m + pt * pt)); });
DECLARE_SOA_DYNAMIC_COLUMN(Y, y, //!
                           [](float pt, float eta, float m) -> float { return std::log((std::sqrt(m * m + pt * pt * std::cosh(eta) * std::cosh(eta)) + pt * std::sinh(eta)) / std::sqrt(m * m + pt * pt)); });
} // namespace reducedpair

DECLARE_SOA_TABLE_STAGED(Dielectrons, "RTDIELECTRON", //!
                         o2::soa::Index<>, reducedpair::ReducedEventId,
                         reducedpair::Mass, reducedpair::Pt, reducedpair::Eta, reducedpair::Phi, reducedpair::Sign,
                         reducedpair::FilterMap, reducedpair::McDecision,
                         reducedpair::Rap<reducedpair::Pt, reducedpair::Eta, reducedpair::Mass>,
                         reducedpair::Y<reducedpair::Pt, reducedpair::Eta, reducedpair::Mass>,
                         reducedpair::Px<reducedpair::Pt, reducedpair::Phi>,
                         reducedpair::Py<reducedpair::Pt, reducedpair::Phi>,
                         reducedpair::Pz<reducedpair::Pt, reducedpair::Eta>,
                         reducedpair::P<reducedpair::Pt, reducedpair::Eta>);

DECLARE_SOA_TABLE(Dimuons, "AOD", "RTDIMUON", //!
                  o2::soa::Index<>, reducedpair::ReducedEventId,
                  reducedpair::Mass, reducedpair::Pt, reducedpair::Eta, reducedpair::Phi, reducedpair::Sign,
                  reducedpair::FilterMap, reducedpair::McDecision,
                  reducedpair::Px<reducedpair::Pt, reducedpair::Phi>,
                  reducedpair::Py<reducedpair::Pt, reducedpair::Phi>,
                  reducedpair::Pz<reducedpair::Pt, reducedpair::Eta>,
                  reducedpair::P<reducedpair::Pt, reducedpair::Eta>,
                  reducedpair::Rap<reducedpair::Pt, reducedpair::Eta, reducedpair::Mass>,
                  reducedpair::Y<reducedpair::Pt, reducedpair::Eta, reducedpair::Mass>);

DECLARE_SOA_TABLE(ElectronMuons, "AOD", "RTELECTRONMUON", //!
                  o2::soa::Index<>, reducedpair::ReducedEventId,
                  reducedpair::Mass, reducedpair::Pt, reducedpair::Eta, reducedpair::Phi, reducedpair::Sign,
                  reducedpair::FilterMap, reducedpair::McDecision,
                  reducedpair::Px<reducedpair::Pt, reducedpair::Phi>,
                  reducedpair::Py<reducedpair::Pt, reducedpair::Phi>,
                  reducedpair::Pz<reducedpair::Pt, reducedpair::Eta>,
                  reducedpair::P<reducedpair::Pt, reducedpair::Eta>,
                  reducedpair::Rap<reducedpair::Pt, reducedpair::Eta, reducedpair::Mass>,
                  reducedpair::Y<reducedpair::Pt, reducedpair::Eta, reducedpair::Mass>);

DECLARE_SOA_TABLE(DielectronsExtra, "AOD", "RTDIELEEXTRA", //!
                  reducedpair::Index0Id, reducedpair::Index1Id,
                  reducedpair::Tauz,
                  reducedpair::Lz,
                  reducedpair::Lxy);

DECLARE_SOA_TABLE(DielectronsInfo, "AOD", "RTDIELINFO",
                  reducedpair::CollisionId, reducedpair::Prong0Id, reducedpair::Prong1Id);

DECLARE_SOA_TABLE(DimuonsExtra, "AOD", "RTDIMUEXTRA", //!
                  dilepton_track_index::Index0Id, dilepton_track_index::Index1Id,
                  reducedpair::Tauz,
                  reducedpair::Lz,
                  reducedpair::Lxy);

DECLARE_SOA_TABLE(DileptonsFlow, "AOD", "RTDILEPTONFLOW", //!
                  reducedpair::CollisionId,
                  reducedpair::Mass,
                  reducedpair::CentFT0C,
                  reducedpair::Pt, reducedpair::Eta, reducedpair::Phi, reducedpair::Sign,
                  reducedpair::IsFirst,
                  reducedpair::U2Q2, reducedpair::R2SP_AB, reducedpair::R2SP_AC, reducedpair::R2SP_BC,
                  reducedpair::U3Q3, reducedpair::R3SP,
                  reducedpair::Cos2DeltaPhi, reducedpair::R2EP_AB, reducedpair::R2EP_AC, reducedpair::R2EP_BC,
                  reducedpair::Cos3DeltaPhi, reducedpair::R3EP,
                  reducedpair::CORR2POI, reducedpair::CORR4POI, reducedpair::M01POI, reducedpair::M0111POI,
                  reducedevent::CORR2REF, reducedevent::CORR4REF, reducedevent::M11REF, reducedevent::M1111REF,
                  reducedpair::MultDimuons, reducedevent::MultA);

// Dilepton collision information (joined with DileptonsExtra) allowing to connect different tables (cross PWGs)
DECLARE_SOA_TABLE(DileptonsInfo, "AOD", "RTDILEPTONINFO",
                  reducedpair::CollisionId, collision::PosX, collision::PosY, collision::PosZ);

DECLARE_SOA_TABLE_STAGED(DielectronsAll, "RTDIELECTRONALL", //!
                         reducedpair::Mass,
                         reducedpair::Pt, reducedpair::Eta, reducedpair::Phi, reducedpair::Sign,
                         reducedpair::FilterMap,
                         reducedpair::McDecision,
                         dilepton_track_index::Pt1, dilepton_track_index::Eta1, dilepton_track_index::Phi1, dilepton_track_index::ITSClusterMap1, dilepton_track_index::ITSChi2NCl1, dilepton_track_index::TPCNClsCR1, dilepton_track_index::TPCNClsFound1, dilepton_track_index::TPCChi2NCl1, dilepton_track_index::DcaXY1, dilepton_track_index::DcaZ1, dilepton_track_index::TPCSignal1, dilepton_track_index::TPCNSigmaEl1, dilepton_track_index::TPCNSigmaPi1, dilepton_track_index::TPCNSigmaPr1, dilepton_track_index::TOFBeta1, dilepton_track_index::TOFNSigmaEl1, dilepton_track_index::TOFNSigmaPi1, dilepton_track_index::TOFNSigmaPr1,
                         dilepton_track_index::Pt2, dilepton_track_index::Eta2, dilepton_track_index::Phi2, dilepton_track_index::ITSClusterMap2, dilepton_track_index::ITSChi2NCl2, dilepton_track_index::TPCNClsCR2, dilepton_track_index::TPCNClsFound2, dilepton_track_index::TPCChi2NCl2, dilepton_track_index::DcaXY2, dilepton_track_index::DcaZ2, dilepton_track_index::TPCSignal2, dilepton_track_index::TPCNSigmaEl2, dilepton_track_index::TPCNSigmaPi2, dilepton_track_index::TPCNSigmaPr2, dilepton_track_index::TOFBeta2, dilepton_track_index::TOFNSigmaEl2, dilepton_track_index::TOFNSigmaPi2, dilepton_track_index::TOFNSigmaPr2,
                         dilepton_track_index::DCAxyzTrk0KF, dilepton_track_index::DCAxyzTrk1KF, reducedpair::DCAxyzBetweenTrksKF, dilepton_track_index::DCAxyTrk0KF, dilepton_track_index::DCAxyTrk1KF, reducedpair::DCAxyBetweenTrksKF,
                         dilepton_track_index::DeviationTrk0KF, dilepton_track_index::DeviationTrk1KF, dilepton_track_index::DeviationxyTrk0KF, dilepton_track_index::DeviationxyTrk1KF,
                         reducedpair::MassKFGeo, reducedpair::Chi2OverNDFKFGeo, reducedpair::DecayLengthKFGeo, reducedpair::DecayLengthOverErrKFGeo, reducedpair::DecayLengthXYKFGeo, reducedpair::DecayLengthXYOverErrKFGeo, reducedpair::PseudoproperDecayTimeKFGeo, reducedpair::PseudoproperDecayTimeErrKFGeo, reducedpair::CosPAKFGeo, reducedpair::PairDCAxyz, reducedpair::PairDCAxy,
                         reducedpair::DeviationPairKF, reducedpair::DeviationxyPairKF,
                         reducedpair::MassKFGeoTop, reducedpair::Chi2OverNDFKFGeoTop,
                         reducedpair::Tauz, reducedpair::Tauxy,
                         reducedpair::Lz,
                         reducedpair::Lxy);

DECLARE_SOA_TABLE(DimuonsAll, "AOD", "RTDIMUONALL", //!
                  collision::PosX, collision::PosY, collision::PosZ, collision::NumContrib,
                  evsel::Selection, reducedpair::EventSelection,
                  reducedevent::MCPosX, reducedevent::MCPosY, reducedevent::MCPosZ,
                  reducedpair::Mass,
                  reducedpair::McDecision,
                  reducedpair::Pt, reducedpair::Eta, reducedpair::Phi, reducedpair::Sign, reducedpair::Chi2pca,
                  reducedpair::Tauz, reducedpair::TauzErr,
                  reducedpair::Tauxy, reducedpair::TauxyErr,
                  reducedpair::CosPointingAngle,
                  dilepton_track_index::Pt1, dilepton_track_index::Eta1, dilepton_track_index::Phi1, dilepton_track_index::Sign1,
                  dilepton_track_index::Pt2, dilepton_track_index::Eta2, dilepton_track_index::Phi2, dilepton_track_index::Sign2,
                  dilepton_track_index::FwdDcaX1, dilepton_track_index::FwdDcaY1, dilepton_track_index::FwdDcaX2, dilepton_track_index::FwdDcaY2,
                  dilepton_track_index::McMask1, dilepton_track_index::McMask2,
                  dilepton_track_index::Chi2MatchMCHMID1, dilepton_track_index::Chi2MatchMCHMID2,
                  dilepton_track_index::Chi2MatchMCHMFT1, dilepton_track_index::Chi2MatchMCHMFT2,
                  dilepton_track_index::Chi21, dilepton_track_index::Chi22,
                  dilepton_track_index::PtMC1, dilepton_track_index::EtaMC1, dilepton_track_index::PhiMC1, dilepton_track_index::EMC1,
                  dilepton_track_index::PtMC2, dilepton_track_index::EtaMC2, dilepton_track_index::PhiMC2, dilepton_track_index::EMC2,
                  dilepton_track_index::Vx1, dilepton_track_index::Vy1, dilepton_track_index::Vz1, dilepton_track_index::Vt1,
                  dilepton_track_index::Vx2, dilepton_track_index::Vy2, dilepton_track_index::Vz2, dilepton_track_index::Vt2,
                  dilepton_track_index::IsAmbig1, dilepton_track_index::IsAmbig2,
                  dilepton_track_index::IsCorrectAssoc1, dilepton_track_index::IsCorrectAssoc2,
                  reducedpair::U2Q2,
                  reducedpair::U3Q3,
                  reducedpair::R2EP_AB,
                  reducedpair::R2SP_AB,
                  reducedpair::CentFT0C,
                  reducedpair::Cos2DeltaPhi,
                  reducedpair::Cos3DeltaPhi,
                  reducedpair::CORR2POI,
                  reducedpair::CORR4POI,
                  reducedpair::M01POI,
                  reducedpair::M0111POI,
                  reducedpair::MultDimuons,
                  reducedpair::VertexPz,
                  reducedpair::SVertex);

DECLARE_SOA_TABLE(DileptonsMiniTree, "AOD", "RTDILEPTMTREE", //!
                  reducedpair::Mass, reducedpair::Pt, reducedpair::Eta, reducedpair::CentFT0C, reducedpair::Cos2DeltaPhi,
                  dilepton_track_index::Pt1, dilepton_track_index::Eta1, dilepton_track_index::Phi1,
                  dilepton_track_index::Pt2, dilepton_track_index::Eta2, dilepton_track_index::Phi2);

DECLARE_SOA_TABLE(DileptonsMiniTreeGen, "AOD", "RTDILMTREEGEN", //!
                  reducedpair::McDecision, mccollision::ImpactParameter,
                  dilepton_track_index::PtMC1, dilepton_track_index::EtaMC1, dilepton_track_index::PhiMC1,
                  dilepton_track_index::PtMC2, dilepton_track_index::EtaMC2, dilepton_track_index::PhiMC2);

DECLARE_SOA_TABLE(DileptonsMiniTreeRec, "AOD", "RTDILMTREEREC", //!
                  reducedpair::McDecision, reducedpair::Mass, reducedpair::Pt, reducedpair::Eta, reducedpair::Phi, reducedpair::CentFT0C,
                  dilepton_track_index::PtMC1, dilepton_track_index::EtaMC1, dilepton_track_index::PhiMC1,
                  dilepton_track_index::PtMC2, dilepton_track_index::EtaMC2, dilepton_track_index::PhiMC2,
                  dilepton_track_index::Pt1, dilepton_track_index::Eta1, dilepton_track_index::Phi1,
                  dilepton_track_index::Pt2, dilepton_track_index::Eta2, dilepton_track_index::Phi2);

DECLARE_SOA_TABLE(DileptonsPolarization, "AOD", "RTDILPOLAR", //!
                  reducedpair::CosThetaHE, reducedpair::PhiHE, reducedpair::PhiTildeHE,
                  reducedpair::CosThetaCS, reducedpair::PhiCS, reducedpair::PhiTildeCS,
                  reducedpair::CosThetaPP, reducedpair::PhiPP, reducedpair::PhiTildePP,
                  reducedpair::CosThetaRM,
                  reducedpair::CosThetaStarTPC, reducedpair::CosThetaStarFT0A, reducedpair::CosThetaStarFT0C);

using Dielectron = Dielectrons::iterator;
using StoredDielectron = StoredDielectrons::iterator;
using Dimuon = Dimuons::iterator;
using DielectronExtra = DielectronsExtra::iterator;
using DielectronInfo = DielectronsInfo::iterator;
using DimuonExtra = DimuonsExtra::iterator;
using DileptonFlow = DileptonsFlow::iterator;
using DileptonInfo = DileptonsInfo::iterator;
using DielectronAll = DielectronsAll::iterator;
using DimuonAll = DimuonsAll::iterator;
using DileptonMiniTree = DileptonsMiniTree::iterator;
using DileptonMiniTreeGen = DileptonsMiniTreeGen::iterator;
using DileptonMiniTreeRec = DileptonsMiniTreeRec::iterator;
using DileptonPolarization = DileptonsPolarization::iterator;

// Tables for using analysis-dilepton-track with analysis-asymmetric-pairing
DECLARE_SOA_TABLE(Ditracks, "AOD", "RTDITRACK", //!
                  o2::soa::Index<>, reducedpair::ReducedEventId,
                  reducedpair::Mass, reducedpair::Pt, reducedpair::Eta, reducedpair::Phi, reducedpair::Sign,
                  reducedpair::FilterMap, reducedpair::PairFilterMap, reducedpair::CommonFilterMap,
                  reducedpair::Rap<reducedpair::Pt, reducedpair::Eta, reducedpair::Mass>,
                  reducedpair::Y<reducedpair::Pt, reducedpair::Eta, reducedpair::Mass>,
                  reducedpair::Px<reducedpair::Pt, reducedpair::Phi>,
                  reducedpair::Py<reducedpair::Pt, reducedpair::Phi>,
                  reducedpair::Pz<reducedpair::Pt, reducedpair::Eta>,
                  reducedpair::P<reducedpair::Pt, reducedpair::Eta>);

DECLARE_SOA_TABLE(DitracksExtra, "AOD", "RTDITRKEXTRA", //!
                  reducedpair::Index0Id, reducedpair::Index1Id,
                  reducedpair::Tauz,
                  reducedpair::Lz,
                  reducedpair::Lxy,
                  o2::soa::Marker<1>);

// mft PID reduced data model
namespace fwdpid
{
DECLARE_SOA_COLUMN(Pt, pt, float);                    //!
DECLARE_SOA_COLUMN(Eta, eta, float);                  //!
DECLARE_SOA_COLUMN(Phi, phi, float);                  //!
DECLARE_SOA_COLUMN(Sign, sign, int);                  //!
DECLARE_SOA_COLUMN(McDecision, mcDecision, uint32_t); //!
} // namespace fwdpid

DECLARE_SOA_TABLE(FwdPidsAll, "AOD", "RTFWDPIDALL", //!
                  fwdtrack::TrackType, collision::PosX, collision::PosY, collision::PosZ, collision::NumContrib,
                  fwdpid::Pt, fwdpid::Eta, fwdpid::Phi, fwdpid::Sign,
                  reducedmft::MftClusterSizesAndTrackFlags,
                  reducedmft::FwdDcaX, reducedmft::FwdDcaY, fwdtrack::Chi2MatchMCHMID, fwdtrack::Chi2MatchMCHMFT, fwdpid::McDecision);

using FwdPidAll = FwdPidsAll::iterator;

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

// candidate information
namespace dileptonTrackTrackCandidate
{
// infotmation about the dilepton-track-track
DECLARE_SOA_COLUMN(Mass, mass, float);                                       //!
DECLARE_SOA_COLUMN(Pt, pt, float);                                           //!
DECLARE_SOA_COLUMN(Eta, eta, float);                                         //!
DECLARE_SOA_COLUMN(Phi, phi, float);                                         //!
DECLARE_SOA_COLUMN(Rap, rap, float);                                         //!
DECLARE_SOA_COLUMN(DeltaQ, deltaQ, float);                                   //!
DECLARE_SOA_COLUMN(R1, r1, float);                                           //! distance between the dilepton and the track1 in theta-phi plane
DECLARE_SOA_COLUMN(R2, r2, float);                                           //! distance between the dilepton and the track2 in theta-phi plane
DECLARE_SOA_COLUMN(R, r, float);                                             //!
DECLARE_SOA_COLUMN(DileptonMass, dileptonMass, float);                       //!
DECLARE_SOA_COLUMN(DileptonPt, dileptonPt, float);                           //!
DECLARE_SOA_COLUMN(DileptonEta, dileptonEta, float);                         //!
DECLARE_SOA_COLUMN(DileptonPhi, dileptonPhi, float);                         //!
DECLARE_SOA_COLUMN(DileptonSign, dileptonSign, int);                         //!
DECLARE_SOA_COLUMN(DileptonTPCnSigmaEl1, dileptonTPCnSigmaEl1, float);       //!
DECLARE_SOA_COLUMN(DileptonTPCnSigmaPi1, dileptonTPCnSigmaPi1, float);       //!
DECLARE_SOA_COLUMN(DileptonTPCnSigmaPr1, dileptonTPCnSigmaPr1, float);       //!
DECLARE_SOA_COLUMN(DileptonTPCnCls1, dileptonTPCnCls1, float);               //!
DECLARE_SOA_COLUMN(DileptonTPCnSigmaEl2, dileptonTPCnSigmaEl2, float);       //!
DECLARE_SOA_COLUMN(DileptonTPCnSigmaPi2, dileptonTPCnSigmaPi2, float);       //!
DECLARE_SOA_COLUMN(DileptonTPCnSigmaPr2, dileptonTPCnSigmaPr2, float);       //!
DECLARE_SOA_COLUMN(DileptonTPCnCls2, dileptonTPCnCls2, float);               //!
DECLARE_SOA_COLUMN(DiTracksMass, diTracksMass, float);                       //!
DECLARE_SOA_COLUMN(DiTracksPt, diTracksPt, float);                           //!
DECLARE_SOA_COLUMN(TrackPt1, trackPt1, float);                               //!
DECLARE_SOA_COLUMN(TrackPt2, trackPt2, float);                               //!
DECLARE_SOA_COLUMN(TrackEta1, trackEta1, float);                             //!
DECLARE_SOA_COLUMN(TrackEta2, trackEta2, float);                             //!
DECLARE_SOA_COLUMN(TrackPhi1, trackPhi1, float);                             //!
DECLARE_SOA_COLUMN(TrackPhi2, trackPhi2, float);                             //!
DECLARE_SOA_COLUMN(TrackSign1, trackSign1, int);                             //!
DECLARE_SOA_COLUMN(TrackSign2, trackSign2, int);                             //!
DECLARE_SOA_COLUMN(TrackTPCNSigmaPi1, trackTPCNSigmaPi1, float);             //!
DECLARE_SOA_COLUMN(TrackTPCNSigmaPi2, trackTPCNSigmaPi2, float);             //!
DECLARE_SOA_COLUMN(TrackTPCNSigmaKa1, trackTPCNSigmaKa1, float);             //!
DECLARE_SOA_COLUMN(TrackTPCNSigmaKa2, trackTPCNSigmaKa2, float);             //!
DECLARE_SOA_COLUMN(TrackTPCNSigmaPr1, trackTPCNSigmaPr1, float);             //!
DECLARE_SOA_COLUMN(TrackTPCNSigmaPr2, trackTPCNSigmaPr2, float);             //!
DECLARE_SOA_COLUMN(TrackTPCNCls1, trackTPCNCls1, float);                     //!
DECLARE_SOA_COLUMN(TrackTPCNCls2, trackTPCNCls2, float);                     //!
DECLARE_SOA_COLUMN(KFMass, kfMass, float);                                   //!
DECLARE_SOA_COLUMN(VertexingProcCode, vertexingProcCode, float);             //!
DECLARE_SOA_COLUMN(VertexingChi2PCA, vertexingChi2PCA, float);               //!
DECLARE_SOA_COLUMN(CosPointingAngle, cosPointingAngle, float);               //!
DECLARE_SOA_COLUMN(KFDCAxyzBetweenProngs, kfDCAxyzBetweenProngs, float);     //!
DECLARE_SOA_COLUMN(KFChi2OverNDFGeo, kfChi2OverNDFGeo, float);               //!
DECLARE_SOA_COLUMN(VertexingLz, vertexingLz, float);                         //!
DECLARE_SOA_COLUMN(VertexingLxy, vertexingLxy, float);                       //!
DECLARE_SOA_COLUMN(VertexingLxyz, vertexingLxyz, float);                     //!
DECLARE_SOA_COLUMN(VertexingTauz, vertexingTauz, float);                     //!
DECLARE_SOA_COLUMN(VertexingTauxy, vertexingTauxy, float);                   //!
DECLARE_SOA_COLUMN(VertexingLzErr, vertexingLzErr, float);                   //!
DECLARE_SOA_COLUMN(VertexingLxyzErr, vertexingLxyzErr, float);               //!
DECLARE_SOA_COLUMN(VertexingTauzErr, vertexingTauzErr, float);               //!
DECLARE_SOA_COLUMN(VertexingLzProjected, vertexingLzProjected, float);       //!
DECLARE_SOA_COLUMN(VertexingLxyProjected, vertexingLxyProjected, float);     //!
DECLARE_SOA_COLUMN(VertexingLxyzProjected, vertexingLxyzProjected, float);   //!
DECLARE_SOA_COLUMN(VertexingTauzProjected, vertexingTauzProjected, float);   //!
DECLARE_SOA_COLUMN(VertexingTauxyProjected, vertexingTauxyProjected, float); //!
} // namespace dileptonTrackTrackCandidate

DECLARE_SOA_TABLE(DileptonTrackTrackCandidates, "AOD", "RTDQUADPLET", //!
                  dileptonTrackTrackCandidate::Mass,
                  dileptonTrackTrackCandidate::Pt,
                  dileptonTrackTrackCandidate::Eta,
                  dileptonTrackTrackCandidate::Phi,
                  dileptonTrackTrackCandidate::Rap,
                  dileptonTrackTrackCandidate::DeltaQ,
                  dileptonTrackTrackCandidate::R1,
                  dileptonTrackTrackCandidate::R2,
                  dileptonTrackTrackCandidate::R,
                  dileptonTrackTrackCandidate::DileptonMass,
                  dileptonTrackTrackCandidate::DileptonPt,
                  dileptonTrackTrackCandidate::DileptonEta,
                  dileptonTrackTrackCandidate::DileptonPhi,
                  dileptonTrackTrackCandidate::DileptonSign,
                  dileptonTrackTrackCandidate::DileptonTPCnSigmaEl1,
                  dileptonTrackTrackCandidate::DileptonTPCnSigmaPi1,
                  dileptonTrackTrackCandidate::DileptonTPCnSigmaPr1,
                  dileptonTrackTrackCandidate::DileptonTPCnCls1,
                  dileptonTrackTrackCandidate::DileptonTPCnSigmaEl2,
                  dileptonTrackTrackCandidate::DileptonTPCnSigmaPi2,
                  dileptonTrackTrackCandidate::DileptonTPCnSigmaPr2,
                  dileptonTrackTrackCandidate::DileptonTPCnCls2,
                  dileptonTrackTrackCandidate::DiTracksMass,
                  dileptonTrackTrackCandidate::DiTracksPt,
                  dileptonTrackTrackCandidate::TrackPt1,
                  dileptonTrackTrackCandidate::TrackPt2,
                  dileptonTrackTrackCandidate::TrackEta1,
                  dileptonTrackTrackCandidate::TrackEta2,
                  dileptonTrackTrackCandidate::TrackPhi1,
                  dileptonTrackTrackCandidate::TrackPhi2,
                  dileptonTrackTrackCandidate::TrackSign1,
                  dileptonTrackTrackCandidate::TrackSign2,
                  dileptonTrackTrackCandidate::TrackTPCNSigmaPi1,
                  dileptonTrackTrackCandidate::TrackTPCNSigmaPi2,
                  dileptonTrackTrackCandidate::TrackTPCNSigmaKa1,
                  dileptonTrackTrackCandidate::TrackTPCNSigmaKa2,
                  dileptonTrackTrackCandidate::TrackTPCNSigmaPr1,
                  dileptonTrackTrackCandidate::TrackTPCNSigmaPr2,
                  dileptonTrackTrackCandidate::TrackTPCNCls1,
                  dileptonTrackTrackCandidate::TrackTPCNCls2,
                  dileptonTrackTrackCandidate::KFMass,
                  dileptonTrackTrackCandidate::VertexingProcCode,
                  dileptonTrackTrackCandidate::VertexingChi2PCA,
                  dileptonTrackTrackCandidate::CosPointingAngle,
                  dileptonTrackTrackCandidate::KFDCAxyzBetweenProngs,
                  dileptonTrackTrackCandidate::KFChi2OverNDFGeo,
                  dileptonTrackTrackCandidate::VertexingLz,
                  dileptonTrackTrackCandidate::VertexingLxy,
                  dileptonTrackTrackCandidate::VertexingLxyz,
                  dileptonTrackTrackCandidate::VertexingTauz,
                  dileptonTrackTrackCandidate::VertexingTauxy,
                  dileptonTrackTrackCandidate::VertexingLzErr,
                  dileptonTrackTrackCandidate::VertexingLxyzErr,
                  dileptonTrackTrackCandidate::VertexingTauzErr,
                  dileptonTrackTrackCandidate::VertexingLzProjected,
                  dileptonTrackTrackCandidate::VertexingLxyProjected,
                  dileptonTrackTrackCandidate::VertexingLxyzProjected,
                  dileptonTrackTrackCandidate::VertexingTauzProjected,
                  dileptonTrackTrackCandidate::VertexingTauxyProjected);

using DileptonTrackTrackCandidate = DileptonTrackTrackCandidates::iterator;

namespace v0bits
{
DECLARE_SOA_COLUMN(PIDBit, pidbit, uint8_t); //!
} // namespace v0bits

// bit information for particle species.
DECLARE_SOA_TABLE(V0Bits, "AOD", "V0BITS", //!
                  v0bits::PIDBit);

// iterators
using V0Bit = V0Bits::iterator;

namespace v0mapID
{
DECLARE_SOA_COLUMN(V0AddID, v0addid, int8_t); //!
} // namespace v0mapID

DECLARE_SOA_TABLE(V0MapID, "AOD", "V0MAPID", //!
                  v0mapID::V0AddID);

namespace cascmapID
{
DECLARE_SOA_COLUMN(CascAddID, cascaddid, int8_t); //!
} // namespace cascmapID

DECLARE_SOA_TABLE(CascMapID, "AOD", "CASCMAPID", //!
                  cascmapID::CascAddID);

namespace DalBits
{
DECLARE_SOA_COLUMN(DALITZBits, dalitzBits, uint8_t); //!
} // namespace DalBits

// bit information for particle species.
DECLARE_SOA_TABLE(DalitzBits, "AOD", "DALITZBITS", DalBits::DALITZBits);

DECLARE_SOA_TABLE(RedJpDmColls, "AOD", "REDJPDMCOLL", //!
                  o2::soa::Index<>,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  collision::NumContrib);

namespace jpsidmescorr
{
DECLARE_SOA_INDEX_COLUMN(RedJpDmColl, redJpDmColl);                                    //!
DECLARE_SOA_COLUMN(MassDmes, massDmes, float);                                         //!
DECLARE_SOA_COLUMN(MassD0, massD0, float);                                             //!
DECLARE_SOA_COLUMN(MassD0bar, massD0bar, float);                                       //!
DECLARE_SOA_COLUMN(Px, px, float);                                                     //!
DECLARE_SOA_COLUMN(Py, py, float);                                                     //!
DECLARE_SOA_COLUMN(Pz, pz, float);                                                     //!
DECLARE_SOA_COLUMN(DecVtxX, decVtxX, float);                                           //!
DECLARE_SOA_COLUMN(DecVtxY, decVtxY, float);                                           //!
DECLARE_SOA_COLUMN(DecVtxZ, decVtxZ, float);                                           //!
DECLARE_SOA_COLUMN(BdtBkgMassHypo0, bdtBkgMassHypo0, float);                           //!
DECLARE_SOA_COLUMN(BdtPromptMassHypo0, bdtPromptMassHypo0, float);                     //!
DECLARE_SOA_COLUMN(BdtNonpromptMassHypo0, bdtNonpromptMassHypo0, float);               //!
DECLARE_SOA_COLUMN(BdtBkg, bdtBkg, float);                                             //!
DECLARE_SOA_COLUMN(BdtPrompt, bdtPrompt, float);                                       //!
DECLARE_SOA_COLUMN(BdtNonprompt, bdtNonprompt, float);                                 //!
DECLARE_SOA_COLUMN(BdtBkgMassHypo1, bdtBkgMassHypo1, float);                           //!
DECLARE_SOA_COLUMN(BdtPromptMassHypo1, bdtPromptMassHypo1, float);                     //!
DECLARE_SOA_COLUMN(BdtNonpromptMassHypo1, bdtNonpromptMassHypo1, float);               //!
DECLARE_SOA_COLUMN(NumColls, numColls, uint64_t);                                      //!
DECLARE_SOA_COLUMN(PtDmes, ptDmes, float);                                             //!
DECLARE_SOA_COLUMN(PtJpsi, ptJpsi, float);                                             //!
DECLARE_SOA_COLUMN(RapDmes, rapDmes, float);                                           //!
DECLARE_SOA_COLUMN(RapJpsi, rapJpsi, float);                                           //!
DECLARE_SOA_COLUMN(PhiDmes, phiDmes, float);                                           //!
DECLARE_SOA_COLUMN(PhiJpsi, phiJpsi, float);                                           //!
DECLARE_SOA_COLUMN(DeltaY, deltaY, float);                                             //!
DECLARE_SOA_COLUMN(DeltaPhi, deltaPhi, float);                                         //!
DECLARE_SOA_COLUMN(NumItsClsDmesProng0, numItsClsDmesProng0, int);                     //!
DECLARE_SOA_COLUMN(NumItsClsDmesProng1, numItsClsDmesProng1, int);                     //!
DECLARE_SOA_COLUMN(NumTpcCrossedRowsDmesProng0, numTpcCrossedRowsDmesProng0, int);     //!
DECLARE_SOA_COLUMN(NumTpcCrossedRowsDmesProng1, numTpcCrossedRowsDmesProng1, int);     //!
DECLARE_SOA_COLUMN(EtaDmesProng0, etaDmesProng0, float);                               //!
DECLARE_SOA_COLUMN(EtaDmesProng1, etaDmesProng1, float);                               //!
DECLARE_SOA_COLUMN(PtDmesProng0, ptDmesProng0, float);                                 //!
DECLARE_SOA_COLUMN(PtDmesProng1, ptDmesProng1, float);                                 //!
DECLARE_SOA_COLUMN(MinNumItsClsDmesProng, minNumItsClsDmesProng, int);                 //!
DECLARE_SOA_COLUMN(MinNumTpcCrossedRowsDmesProng, minNumTpcCrossedRowsDmesProng, int); //!
DECLARE_SOA_COLUMN(MinAbsEtaDmesProng, minAbsEtaDmesProng, float);                     //!
DECLARE_SOA_COLUMN(MinPtDmesProng, minPtDmesProng, float);                             //!
DECLARE_SOA_COLUMN(NumSigmaTpcPiProng0, numSigmaTpcPiProng0, float);                   //!
DECLARE_SOA_COLUMN(NumSigmaTpcPiProng1, numSigmaTpcPiProng1, float);                   //!
DECLARE_SOA_COLUMN(NumSigmaTofPiProng0, numSigmaTofPiProng0, float);                   //!
DECLARE_SOA_COLUMN(NumSigmaTofPiProng1, numSigmaTofPiProng1, float);                   //!
DECLARE_SOA_COLUMN(NumSigmaTpcKaProng0, numSigmaTpcKaProng0, float);                   //!
DECLARE_SOA_COLUMN(NumSigmaTpcKaProng1, numSigmaTpcKaProng1, float);                   //!
DECLARE_SOA_COLUMN(NumSigmaTofKaProng0, numSigmaTofKaProng0, float);                   //!
DECLARE_SOA_COLUMN(NumSigmaTofKaProng1, numSigmaTofKaProng1, float);                   //!
DECLARE_SOA_DYNAMIC_COLUMN(NumSigmaTpcTofPiProng0, numSigmaTpcTofPiProng0,             //!
                           [](float tpcNSigmaPi0, float tofNSigmaPi0) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaPi0, tofNSigmaPi0); });
DECLARE_SOA_DYNAMIC_COLUMN(NumSigmaTpcTofPiProng1, numSigmaTpcTofPiProng1, //!
                           [](float tpcNSigmaPi1, float tofNSigmaPi1) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaPi1, tofNSigmaPi1); });
DECLARE_SOA_DYNAMIC_COLUMN(NumSigmaTpcTofKaProng0, numSigmaTpcTofKaProng0, //!
                           [](float tpcNSigmaKa0, float tofNSigmaKa0) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaKa0, tofNSigmaKa0); });
DECLARE_SOA_DYNAMIC_COLUMN(NumSigmaTpcTofKaProng1, numSigmaTpcTofKaProng1, //!
                           [](float tpcNSigmaKa1, float tofNSigmaKa1) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaKa1, tofNSigmaKa1); });
} // namespace jpsidmescorr

DECLARE_SOA_TABLE(RedJpDmDileptons, "AOD", "REDJPDMDILEPTON", //!
                  o2::soa::Index<>,
                  jpsidmescorr::RedJpDmCollId,
                  jpsidmescorr::Px,
                  jpsidmescorr::Py,
                  jpsidmescorr::Pz,
                  reducedpair::Mass,
                  reducedpair::Sign,
                  reducedpair::McDecision,
                  reducedpair::Tauz,
                  reducedpair::Lz,
                  reducedpair::Lxy);

DECLARE_SOA_TABLE(RedJpDmColCounts, "AOD", "REDJPDMCOLCOUNT", //!
                  jpsidmescorr::NumColls);

DECLARE_SOA_TABLE(RedJpDmDmesons, "AOD", "REDJPDMDMESON", //!
                  o2::soa::Index<>,
                  jpsidmescorr::RedJpDmCollId,
                  jpsidmescorr::Px,
                  jpsidmescorr::Py,
                  jpsidmescorr::Pz,
                  jpsidmescorr::DecVtxX,
                  jpsidmescorr::DecVtxY,
                  jpsidmescorr::DecVtxZ,
                  reducedpair::Sign,
                  reducedpair::McDecision);

DECLARE_SOA_TABLE(RedJpDmDmDau0s, "AOD", "REDJPDMDMDAU0", //!
                  jpsidmescorr::PtDmesProng0,
                  jpsidmescorr::EtaDmesProng0,
                  jpsidmescorr::NumItsClsDmesProng0,
                  jpsidmescorr::NumTpcCrossedRowsDmesProng0,
                  jpsidmescorr::NumSigmaTpcPiProng0,
                  jpsidmescorr::NumSigmaTofPiProng0,
                  jpsidmescorr::NumSigmaTpcKaProng0,
                  jpsidmescorr::NumSigmaTofKaProng0,
                  jpsidmescorr::NumSigmaTpcTofPiProng0<jpsidmescorr::NumSigmaTpcPiProng0, jpsidmescorr::NumSigmaTofPiProng0>,
                  jpsidmescorr::NumSigmaTpcTofKaProng0<jpsidmescorr::NumSigmaTpcKaProng0, jpsidmescorr::NumSigmaTofKaProng0>);

DECLARE_SOA_TABLE(RedJpDmDmDau1s, "AOD", "REDJPDMDMDAU1", //!
                  jpsidmescorr::PtDmesProng1,
                  jpsidmescorr::EtaDmesProng1,
                  jpsidmescorr::NumItsClsDmesProng1,
                  jpsidmescorr::NumTpcCrossedRowsDmesProng1,
                  jpsidmescorr::NumSigmaTpcPiProng1,
                  jpsidmescorr::NumSigmaTofPiProng1,
                  jpsidmescorr::NumSigmaTpcKaProng1,
                  jpsidmescorr::NumSigmaTofKaProng1,
                  jpsidmescorr::NumSigmaTpcTofPiProng1<jpsidmescorr::NumSigmaTpcPiProng1, jpsidmescorr::NumSigmaTofPiProng1>,
                  jpsidmescorr::NumSigmaTpcTofKaProng1<jpsidmescorr::NumSigmaTpcKaProng1, jpsidmescorr::NumSigmaTofKaProng1>);

DECLARE_SOA_TABLE(RedJpDmD0Masss, "AOD", "REDJPDMD0MASS", //!
                  jpsidmescorr::MassD0,
                  jpsidmescorr::MassD0bar);

DECLARE_SOA_TABLE(RedJpDmDmesBdts, "AOD", "REDJPDMDMESBDT", //!
                  jpsidmescorr::BdtBkgMassHypo0,
                  jpsidmescorr::BdtPromptMassHypo0,
                  jpsidmescorr::BdtNonpromptMassHypo0,
                  jpsidmescorr::BdtBkgMassHypo1,
                  jpsidmescorr::BdtPromptMassHypo1,
                  jpsidmescorr::BdtNonpromptMassHypo1);

DECLARE_SOA_TABLE(RedDleptDmesAll, "AOD", "RTDILPTDMESALL", //!
                  reducedpair::Mass,
                  jpsidmescorr::MassDmes,
                  jpsidmescorr::PtJpsi,
                  jpsidmescorr::PtDmes,
                  jpsidmescorr::RapJpsi,
                  jpsidmescorr::RapDmes,
                  jpsidmescorr::PhiJpsi,
                  jpsidmescorr::PhiDmes,
                  jpsidmescorr::DeltaY,
                  jpsidmescorr::DeltaPhi,
                  jpsidmescorr::BdtBkg,
                  jpsidmescorr::BdtPrompt,
                  jpsidmescorr::BdtNonprompt,
                  jpsidmescorr::MinPtDmesProng,
                  jpsidmescorr::MinAbsEtaDmesProng,
                  jpsidmescorr::MinNumItsClsDmesProng,
                  jpsidmescorr::MinNumTpcCrossedRowsDmesProng);

namespace muondca
{
DECLARE_SOA_COLUMN(pDCA, pdca, float); //!
DECLARE_SOA_COLUMN(DCA, dca, float);   //!
DECLARE_SOA_COLUMN(DCAx, dcax, float); //!
DECLARE_SOA_COLUMN(DCAy, dcay, float); //!
DECLARE_SOA_COLUMN(Rabs, rabs, float); //!
DECLARE_SOA_COLUMN(Px, px, float);     //!
DECLARE_SOA_COLUMN(Py, py, float);     //!
DECLARE_SOA_COLUMN(Pz, pz, float);     //!
} // namespace muondca

DECLARE_SOA_TABLE(ReducedMuonsDca, "AOD", "RTMUONDCA",
                  muondca::pDCA,
                  muondca::DCA,
                  muondca::DCAx,
                  muondca::DCAy,
                  muondca::Rabs,
                  reducedmuon::Pt,
                  reducedmuon::Eta, reducedmuon::Phi,
                  reducedmuon::Sign, reducedmuon::IsAmbiguous,
                  muondca::Px,
                  muondca::Py,
                  muondca::Pz);

using ReducedMuonDca = ReducedMuonsDca::iterator;

//______________________________________________________
namespace generatedquarkoniamc
{
//______________________________________________________
// Binned content for generated particles: derived data
DECLARE_SOA_COLUMN(GeneratedEtaC1S, generatedEtaC1S, std::vector<uint32_t>); //! Eta(1S) binned generated data
DECLARE_SOA_COLUMN(GeneratedJPsi, generatedJPsi, std::vector<uint32_t>);     //! J/Psi binned generated data
DECLARE_SOA_COLUMN(GeneratedChiC0, generatedChiC0, std::vector<uint32_t>);   //! ChiC0(1P) binned generated data
DECLARE_SOA_COLUMN(GeneratedChiC1, generatedChiC1, std::vector<uint32_t>);   //! ChiC1(1P) binned generated data
DECLARE_SOA_COLUMN(GeneratedHC, generatedHC, std::vector<uint32_t>);         //! hC binned generated data
DECLARE_SOA_COLUMN(GeneratedChiC2, generatedChiC2, std::vector<uint32_t>);   //! ChiC2(1P) binned generated data
DECLARE_SOA_COLUMN(GeneratedEtaC2S, generatedEtaC2S, std::vector<uint32_t>); //! EtaC(2S) binned generated data
DECLARE_SOA_COLUMN(GeneratedPsi2S, generatedPsi2S, std::vector<uint32_t>);   //! Psi(2S) binned generated data
} // namespace generatedquarkoniamc

DECLARE_SOA_TABLE(GeEtaC1S, "AOD", "GEETAC1S", generatedquarkoniamc::GeneratedEtaC1S);
DECLARE_SOA_TABLE(GeJPsi, "AOD", "GEJPSI", generatedquarkoniamc::GeneratedJPsi);
DECLARE_SOA_TABLE(GeChiC0, "AOD", "GECHIC0", generatedquarkoniamc::GeneratedChiC0);
DECLARE_SOA_TABLE(GeChiC1, "AOD", "GECHIC1", generatedquarkoniamc::GeneratedChiC1);
DECLARE_SOA_TABLE(GeHC, "AOD", "GEHC", generatedquarkoniamc::GeneratedHC);
DECLARE_SOA_TABLE(GeChiC2, "AOD", "GECHIC2", generatedquarkoniamc::GeneratedChiC2);
DECLARE_SOA_TABLE(GeEtaC2S, "AOD", "GEETAC2S", generatedquarkoniamc::GeneratedEtaC2S);
DECLARE_SOA_TABLE(GePsi2S, "AOD", "GEPSI2S", generatedquarkoniamc::GeneratedPsi2S);
} // namespace o2::aod

#endif // PWGDQ_DATAMODEL_REDUCEDINFOTABLES_H_
