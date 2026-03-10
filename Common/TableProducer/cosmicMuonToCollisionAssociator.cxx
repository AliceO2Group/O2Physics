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
/// \file   cosmicMuonToCollisionAssociator.cxx
/// \brief  Cosmic Muon tracking for Alice Data
/// \author Rahul Verma (rahul.verma@iitb.ac.in) :: Marian I Ivanov (mivanov@cern.ch)

#include "TableHelper.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/OccupancyTables.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Tools/TrackTuner.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsFT0/Digit.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

#include "TPDGCode.h"
#include <TRandom.h>

#include <chrono>
#include <random>
#include <typeinfo>

namespace o2
{
namespace aod
{
using Tracks = aod::Tracks;

namespace csmu
{
DECLARE_SOA_INDEX_COLUMN_FULL(UTrack, uTrack, int64_t, Tracks, "_UTrk"); //! Positive track
DECLARE_SOA_INDEX_COLUMN_FULL(LTrack, lTrack, int64_t, Tracks, "_LTrk"); //! Negative track

// Bunch Crossing and Time Frame information.
DECLARE_SOA_COLUMN(UTrkGlobalBcStart, uTrkGlobalBcStart, int64_t);
DECLARE_SOA_COLUMN(UTrkGlobalBcEnd, uTrkGlobalBcEnd, int64_t);
DECLARE_SOA_COLUMN(UTimestamp, uTimestamp, int64_t);
DECLARE_SOA_COLUMN(UTfId, uTfId, int64_t);
DECLARE_SOA_COLUMN(UBCinTF, uBCinTF, int64_t);

DECLARE_SOA_COLUMN(LTrkGlobalBcStart, lTrkGlobalBcStart, int64_t);
DECLARE_SOA_COLUMN(LTrkGlobalBcEnd, lTrkGlobalBcEnd, int64_t);
DECLARE_SOA_COLUMN(LTimestamp, lTimestamp, int64_t);
DECLARE_SOA_COLUMN(LTfId, lTfId, int64_t);
DECLARE_SOA_COLUMN(LBCinTF, lBCinTF, int64_t);

// aod::Tracks table columns
DECLARE_SOA_COLUMN(UTrkX, uTrkX, float);
DECLARE_SOA_COLUMN(UTrkAlpha, uTrkAlpha, float);
DECLARE_SOA_COLUMN(UTrkY, uTrkY, float);
DECLARE_SOA_COLUMN(UTrkZ, uTrkZ, float);
DECLARE_SOA_COLUMN(UTrkSnp, uTrkSnp, float);
DECLARE_SOA_COLUMN(UTrkTgl, uTrkTgl, float);
DECLARE_SOA_COLUMN(UTrkSigned1Pt, uTrkSigned1Pt, float);
DECLARE_SOA_COLUMN(UTrkPt, uTrkPt, float);
DECLARE_SOA_COLUMN(UTrkP, uTrkP, float);
DECLARE_SOA_COLUMN(UTrkEta, uTrkEta, float);
DECLARE_SOA_COLUMN(UTrkPhi, uTrkPhi, float);

DECLARE_SOA_COLUMN(LTrkX, lTrkX, float);
DECLARE_SOA_COLUMN(LTrkAlpha, lTrkAlpha, float);
DECLARE_SOA_COLUMN(LTrkY, lTrkY, float);
DECLARE_SOA_COLUMN(LTrkZ, lTrkZ, float);
DECLARE_SOA_COLUMN(LTrkSnp, lTrkSnp, float);
DECLARE_SOA_COLUMN(LTrkTgl, lTrkTgl, float);
DECLARE_SOA_COLUMN(LTrkSigned1Pt, lTrkSigned1Pt, float);
DECLARE_SOA_COLUMN(LTrkPt, lTrkPt, float);
DECLARE_SOA_COLUMN(LTrkP, lTrkP, float);
DECLARE_SOA_COLUMN(LTrkEta, lTrkEta, float);
DECLARE_SOA_COLUMN(LTrkPhi, lTrkPhi, float);

// aod::TracksExtra table columns
DECLARE_SOA_COLUMN(UTrkTPCInnerParam, uTrkTPCInnerParam, float);
DECLARE_SOA_COLUMN(UTrkFlags, uTrkFlags, uint32_t);

DECLARE_SOA_COLUMN(UTrkITSClusterSizes, uTrkITSClusterSizes, uint32_t); // aod::TracksExtra_002 additinal column

DECLARE_SOA_COLUMN(UTrkTPCNClsFindable, uTrkTPCNClsFindable, uint8_t);
DECLARE_SOA_COLUMN(UTrkTPCNClsFindableMinusFound, uTrkTPCNClsFindableMinusFound, int8_t);

DECLARE_SOA_COLUMN(UTrkTPCNClsFindableMinusPID, uTrkTPCNClsFindableMinusPID, int8_t); // aod::TracksExtra_002 additinal column

DECLARE_SOA_COLUMN(UTrkTPCNClsFindableMinusCrossedRows, uTrkTPCNClsFindableMinusCrossedRows, int8_t);
DECLARE_SOA_COLUMN(UTrkTPCNClsShared, uTrkTPCNClsShared, uint8_t);

DECLARE_SOA_COLUMN(UTrkTPCDeltaTFwd, uTrkTPCDeltaTFwd, float); // aod::TracksExtra_002 additinal column
DECLARE_SOA_COLUMN(UTrkTPCDeltaTBwd, uTrkTPCDeltaTBwd, float); // aod::TracksExtra_002 additinal column

DECLARE_SOA_COLUMN(UTrkTRDPattern, uTrkTRDPattern, uint8_t);
DECLARE_SOA_COLUMN(UTrkITSChi2NCl, uTrkITSChi2NCl, float);
DECLARE_SOA_COLUMN(UTrkTPCChi2NCl, uTrkTPCChi2NCl, float);
DECLARE_SOA_COLUMN(UTrkTRDChi2, uTrkTRDChi2, float);
DECLARE_SOA_COLUMN(UTrkTOFChi2, uTrkTOFChi2, float);
DECLARE_SOA_COLUMN(UTrkTPCSignal, uTrkTPCSignal, float);
DECLARE_SOA_COLUMN(UTrkTRDSignal, uTrkTRDSignal, float);
DECLARE_SOA_COLUMN(UTrkLength, uTrkLength, float);
DECLARE_SOA_COLUMN(UTrkTOFExpMom, uTrkTOFExpMom, float);
DECLARE_SOA_COLUMN(UTrkPIDForTracking, uTrkPIDForTracking, uint32_t);
DECLARE_SOA_COLUMN(UTrkIsPVContributor, uTrkIsPVContributor, bool);
DECLARE_SOA_COLUMN(UTrkHasITS, uTrkHasITS, bool);
DECLARE_SOA_COLUMN(UTrkHasTPC, uTrkHasTPC, bool);
DECLARE_SOA_COLUMN(UTrkHasTRD, uTrkHasTRD, bool);
DECLARE_SOA_COLUMN(UTrkHasTOF, uTrkHasTOF, bool);
DECLARE_SOA_COLUMN(UTrkTPCNClsFound, uTrkTPCNClsFound, int16_t);

DECLARE_SOA_COLUMN(UTrkTPCNClsPID, uTrkTPCNClsPID, int16_t); // aod::TracksExtra_002 additinal column

DECLARE_SOA_COLUMN(UTrkTPCNClsCrossedRows, uTrkTPCNClsCrossedRows, int16_t);

DECLARE_SOA_COLUMN(UTrkITSClusterMap, uTrkITSClusterMap, uint8_t); // Version updated in 002
DECLARE_SOA_COLUMN(UTrkITSNCls, uTrkITSNCls, uint8_t);             // Version updated in 002

DECLARE_SOA_COLUMN(UTrkITSNClsInnerBarrel, uTrkITSNClsInnerBarrel, uint8_t); // Version updated in 002
DECLARE_SOA_COLUMN(UTrkITSClsSizeInLayer, uTrkITSClsSizeInLayer, uint8_t);   // aod::TracksExtra_002 additinal column
DECLARE_SOA_COLUMN(UTrkIsITSAfterburner, uTrkIsITSAfterburner, bool);        // aod::TracksExtra_002 additinal column

DECLARE_SOA_COLUMN(UTrkTOFExpTimeEl, uTrkTOFExpTimeEl, float);
DECLARE_SOA_COLUMN(UTrkTOFExpTimeMu, uTrkTOFExpTimeMu, float);
DECLARE_SOA_COLUMN(UTrkTOFExpTimePi, uTrkTOFExpTimePi, float);
DECLARE_SOA_COLUMN(UTrkTOFExpTimeKa, uTrkTOFExpTimeKa, float);
DECLARE_SOA_COLUMN(UTrkTOFExpTimePr, uTrkTOFExpTimePr, float);
DECLARE_SOA_COLUMN(UTrkTOFExpTimeDe, uTrkTOFExpTimeDe, float);
DECLARE_SOA_COLUMN(UTrkTOFExpTimeTr, uTrkTOFExpTimeTr, float);
DECLARE_SOA_COLUMN(UTrkTOFExpTimeHe, uTrkTOFExpTimeHe, float);
DECLARE_SOA_COLUMN(UTrkTOFExpTimeAl, uTrkTOFExpTimeAl, float);

DECLARE_SOA_COLUMN(UTrkTPCCrossedRowsOverFindableCls, uTrkTPCCrossedRowsOverFindableCls, float);
DECLARE_SOA_COLUMN(UTrkTPCFoundOverFindableCls, uTrkTPCFoundOverFindableCls, float);
DECLARE_SOA_COLUMN(UTrkTPCFractionSharedCls, uTrkTPCFractionSharedCls, float);
DECLARE_SOA_COLUMN(UTrkTrackEtaEmcal, uTrkTrackEtaEmcal, float);
DECLARE_SOA_COLUMN(UTrkTrackPhiEmcal, uTrkTrackPhiEmcal, float);
DECLARE_SOA_COLUMN(UTrkTrackTime, uTrkTrackTime, float);
DECLARE_SOA_COLUMN(UTrkTrackTimeRes, uTrkTrackTimeRes, float);
DECLARE_SOA_COLUMN(UTrkDetectorMap, uTrkDetectorMap, uint8_t);

DECLARE_SOA_COLUMN(LTrkTPCInnerParam, lTrkTPCInnerParam, float);
DECLARE_SOA_COLUMN(LTrkFlags, lTrkFlags, uint32_t);

DECLARE_SOA_COLUMN(LTrkITSClusterSizes, lTrkITSClusterSizes, uint32_t); // New in _002

DECLARE_SOA_COLUMN(LTrkTPCNClsFindable, lTrkTPCNClsFindable, uint8_t);
DECLARE_SOA_COLUMN(LTrkTPCNClsFindableMinusFound, lTrkTPCNClsFindableMinusFound, int8_t);

DECLARE_SOA_COLUMN(LTrkTPCNClsFindableMinusPID, lTrkTPCNClsFindableMinusPID, int8_t); // New in _002

DECLARE_SOA_COLUMN(LTrkTPCNClsFindableMinusCrossedRows, lTrkTPCNClsFindableMinusCrossedRows, int8_t);
DECLARE_SOA_COLUMN(LTrkTPCNClsShared, lTrkTPCNClsShared, uint8_t);

DECLARE_SOA_COLUMN(LTrkTPCDeltaTFwd, lTrkTPCDeltaTFwd, float); // New in _002
DECLARE_SOA_COLUMN(LTrkTPCDeltaTBwd, lTrkTPCDeltaTBwd, float); // New in _002

DECLARE_SOA_COLUMN(LTrkTRDPattern, lTrkTRDPattern, uint8_t);
DECLARE_SOA_COLUMN(LTrkITSChi2NCl, lTrkITSChi2NCl, float);
DECLARE_SOA_COLUMN(LTrkTPCChi2NCl, lTrkTPCChi2NCl, float);
DECLARE_SOA_COLUMN(LTrkTRDChi2, lTrkTRDChi2, float);
DECLARE_SOA_COLUMN(LTrkTOFChi2, lTrkTOFChi2, float);
DECLARE_SOA_COLUMN(LTrkTPCSignal, lTrkTPCSignal, float);
DECLARE_SOA_COLUMN(LTrkTRDSignal, lTrkTRDSignal, float);
DECLARE_SOA_COLUMN(LTrkLength, lTrkLength, float);
DECLARE_SOA_COLUMN(LTrkTOFExpMom, lTrkTOFExpMom, float);
DECLARE_SOA_COLUMN(LTrkPIDForTracking, lTrkPIDForTracking, uint32_t);
DECLARE_SOA_COLUMN(LTrkIsPVContributor, lTrkIsPVContributor, bool);
DECLARE_SOA_COLUMN(LTrkHasITS, lTrkHasITS, bool);
DECLARE_SOA_COLUMN(LTrkHasTPC, lTrkHasTPC, bool);
DECLARE_SOA_COLUMN(LTrkHasTRD, lTrkHasTRD, bool);
DECLARE_SOA_COLUMN(LTrkHasTOF, lTrkHasTOF, bool);
DECLARE_SOA_COLUMN(LTrkTPCNClsFound, lTrkTPCNClsFound, int16_t);

DECLARE_SOA_COLUMN(LTrkTPCNClsPID, lTrkTPCNClsPID, int16_t); // New in _002

DECLARE_SOA_COLUMN(LTrkTPCNClsCrossedRows, lTrkTPCNClsCrossedRows, int16_t);

DECLARE_SOA_COLUMN(LTrkITSClusterMap, lTrkITSClusterMap, uint8_t);           // Version updated in _002
DECLARE_SOA_COLUMN(LTrkITSNCls, lTrkITSNCls, uint8_t);                       // Version updated in _002
DECLARE_SOA_COLUMN(LTrkITSNClsInnerBarrel, lTrkITSNClsInnerBarrel, uint8_t); // Version updated in _002
DECLARE_SOA_COLUMN(LTrkITSClsSizeInLayer, lTrkITSClsSizeInLayer, uint8_t);   // New in _002
DECLARE_SOA_COLUMN(LTrkIsITSAfterburner, lTrkIsITSAfterburner, bool);        // New in _002

DECLARE_SOA_COLUMN(LTrkTOFExpTimeEl, lTrkTOFExpTimeEl, float);
DECLARE_SOA_COLUMN(LTrkTOFExpTimeMu, lTrkTOFExpTimeMu, float);
DECLARE_SOA_COLUMN(LTrkTOFExpTimePi, lTrkTOFExpTimePi, float);
DECLARE_SOA_COLUMN(LTrkTOFExpTimeKa, lTrkTOFExpTimeKa, float);
DECLARE_SOA_COLUMN(LTrkTOFExpTimePr, lTrkTOFExpTimePr, float);
DECLARE_SOA_COLUMN(LTrkTOFExpTimeDe, lTrkTOFExpTimeDe, float);
DECLARE_SOA_COLUMN(LTrkTOFExpTimeTr, lTrkTOFExpTimeTr, float);
DECLARE_SOA_COLUMN(LTrkTOFExpTimeHe, lTrkTOFExpTimeHe, float);
DECLARE_SOA_COLUMN(LTrkTOFExpTimeAl, lTrkTOFExpTimeAl, float);

DECLARE_SOA_COLUMN(LTrkTPCCrossedRowsOverFindableCls, lTrkTPCCrossedRowsOverFindableCls, float);
DECLARE_SOA_COLUMN(LTrkTPCFoundOverFindableCls, lTrkTPCFoundOverFindableCls, float);
DECLARE_SOA_COLUMN(LTrkTPCFractionSharedCls, lTrkTPCFractionSharedCls, float);
DECLARE_SOA_COLUMN(LTrkTrackEtaEmcal, lTrkTrackEtaEmcal, float);
DECLARE_SOA_COLUMN(LTrkTrackPhiEmcal, lTrkTrackPhiEmcal, float);
DECLARE_SOA_COLUMN(LTrkTrackTime, lTrkTrackTime, float);
DECLARE_SOA_COLUMN(LTrkTrackTimeRes, lTrkTrackTimeRes, float);
DECLARE_SOA_COLUMN(LTrkDetectorMap, lTrkDetectorMap, uint8_t);

// Marian and Rahul's addition calculated DCA Column
DECLARE_SOA_COLUMN(UTrkDcaCalc, uTrkDcaCalc, float);
DECLARE_SOA_COLUMN(LTrkDcaCalc, lTrkDcaCalc, float);

// aod::TracksDCA table columns
DECLARE_SOA_COLUMN(UTrkDcaXY, uTrkDcaXY, float);
DECLARE_SOA_COLUMN(UTrkDcaZ, uTrkDcaZ, float);

DECLARE_SOA_COLUMN(LTrkDcaXY, lTrkDcaXY, float);
DECLARE_SOA_COLUMN(LTrkDcaZ, lTrkDcaZ, float);

// aod::TracksDCACov table columns
DECLARE_SOA_COLUMN(UTrkSigmaDcaXY2, uTrkSigmaDcaXY2, float);
DECLARE_SOA_COLUMN(UTrkSigmaDcaZ2, uTrkSigmaDcaZ2, float);

DECLARE_SOA_COLUMN(LTrkSigmaDcaXY2, lTrkSigmaDcaXY2, float);
DECLARE_SOA_COLUMN(LTrkSigmaDcaZ2, lTrkSigmaDcaZ2, float);

// aod::pidTOFFlags
DECLARE_SOA_COLUMN(UTrkGoodTOFMatch, uTrkGoodTOFMatch, bool);
DECLARE_SOA_COLUMN(LTrkGoodTOFMatch, lTrkGoodTOFMatch, bool);

// aod::pidEvTimeFlags
DECLARE_SOA_COLUMN(UTrkTOFFlags, uTrkTOFFlags, uint8_t);
DECLARE_SOA_COLUMN(UTrkIsEvTimeDefined, uTrkIsEvTimeDefined, bool);
DECLARE_SOA_COLUMN(UTrkIsEvTimeTOF, uTrkIsEvTimeTOF, bool);
DECLARE_SOA_COLUMN(UTrkIsEvTimeT0AC, uTrkIsEvTimeT0AC, bool);
DECLARE_SOA_COLUMN(UTrkIsEvTimeTOFT0AC, uTrkIsEvTimeTOFT0AC, bool);

DECLARE_SOA_COLUMN(LTrkTOFFlags, lTrkTOFFlags, uint8_t);
DECLARE_SOA_COLUMN(LTrkIsEvTimeDefined, lTrkIsEvTimeDefined, bool);
DECLARE_SOA_COLUMN(LTrkIsEvTimeTOF, lTrkIsEvTimeTOF, bool);
DECLARE_SOA_COLUMN(LTrkIsEvTimeT0AC, lTrkIsEvTimeT0AC, bool);
DECLARE_SOA_COLUMN(LTrkIsEvTimeTOFT0AC, lTrkIsEvTimeTOFT0AC, bool);

// aod::TOFSignal
DECLARE_SOA_COLUMN(UTrkTOFSignal, uTrkTOFSignal, float);
DECLARE_SOA_COLUMN(UTrkEventCollisionTime, uTrkEventCollisionTime, float);

DECLARE_SOA_COLUMN(LTrkTOFSignal, lTrkTOFSignal, float);
DECLARE_SOA_COLUMN(LTrkEventCollisionTime, lTrkEventCollisionTime, float);

// aod::TOFEvTime
DECLARE_SOA_COLUMN(UTrkTOFEvTime, uTrkTOFEvTime, float);
DECLARE_SOA_COLUMN(UTrkTOFEvTimeErr, uTrkTOFEvTimeErr, float);

DECLARE_SOA_COLUMN(LTrkTOFEvTime, lTrkTOFEvTime, float);
DECLARE_SOA_COLUMN(LTrkTOFEvTimeErr, lTrkTOFEvTimeErr, float);

// aod::TrackQA version 000 columns
DECLARE_SOA_COLUMN(UTrkTPCTime0, uTrkTPCTime0, float);
DECLARE_SOA_COLUMN(UTrkTPCdEdxNorm, uTrkTPCdEdxNorm, float); // New in _003
DECLARE_SOA_COLUMN(UTrkTPCDCAR, uTrkTPCDCAR, int16_t);
DECLARE_SOA_COLUMN(UTrkTPCDCAZ, uTrkTPCDCAZ, int16_t);
DECLARE_SOA_COLUMN(UTrkTPCClusterByteMask, uTrkTPCClusterByteMask, uint8_t);
DECLARE_SOA_COLUMN(UTrkTPCdEdxMax0R, uTrkTPCdEdxMax0R, uint8_t);
DECLARE_SOA_COLUMN(UTrkTPCdEdxMax1R, uTrkTPCdEdxMax1R, uint8_t);
DECLARE_SOA_COLUMN(UTrkTPCdEdxMax2R, uTrkTPCdEdxMax2R, uint8_t);
DECLARE_SOA_COLUMN(UTrkTPCdEdxMax3R, uTrkTPCdEdxMax3R, uint8_t);
DECLARE_SOA_COLUMN(UTrkTPCdEdxTot0R, uTrkTPCdEdxTot0R, uint8_t);
DECLARE_SOA_COLUMN(UTrkTPCdEdxTot1R, uTrkTPCdEdxTot1R, uint8_t);
DECLARE_SOA_COLUMN(UTrkTPCdEdxTot2R, uTrkTPCdEdxTot2R, uint8_t);
DECLARE_SOA_COLUMN(UTrkTPCdEdxTot3R, uTrkTPCdEdxTot3R, uint8_t);

// aod::TrackQA version 000 columns for LTrk
DECLARE_SOA_COLUMN(LTrkTPCTime0, lTrkTPCTime0, float);
DECLARE_SOA_COLUMN(LTrkTPCdEdxNorm, lTrkTPCdEdxNorm, float); // New in _003
DECLARE_SOA_COLUMN(LTrkTPCDCAR, lTrkTPCDCAR, int16_t);
DECLARE_SOA_COLUMN(LTrkTPCDCAZ, lTrkTPCDCAZ, int16_t);
DECLARE_SOA_COLUMN(LTrkTPCClusterByteMask, lTrkTPCClusterByteMask, uint8_t);
DECLARE_SOA_COLUMN(LTrkTPCdEdxMax0R, lTrkTPCdEdxMax0R, uint8_t);
DECLARE_SOA_COLUMN(LTrkTPCdEdxMax1R, lTrkTPCdEdxMax1R, uint8_t);
DECLARE_SOA_COLUMN(LTrkTPCdEdxMax2R, lTrkTPCdEdxMax2R, uint8_t);
DECLARE_SOA_COLUMN(LTrkTPCdEdxMax3R, lTrkTPCdEdxMax3R, uint8_t);
DECLARE_SOA_COLUMN(LTrkTPCdEdxTot0R, lTrkTPCdEdxTot0R, uint8_t);
DECLARE_SOA_COLUMN(LTrkTPCdEdxTot1R, lTrkTPCdEdxTot1R, uint8_t);
DECLARE_SOA_COLUMN(LTrkTPCdEdxTot2R, lTrkTPCdEdxTot2R, uint8_t);
DECLARE_SOA_COLUMN(LTrkTPCdEdxTot3R, lTrkTPCdEdxTot3R, uint8_t);

// aod::TrackQA version 003 columns
DECLARE_SOA_COLUMN(UTrkDeltaRefContParamY, uTrkDeltaRefContParamY, int8_t);
DECLARE_SOA_COLUMN(UTrkDeltaRefContParamZ, uTrkDeltaRefContParamZ, int8_t);
DECLARE_SOA_COLUMN(UTrkDeltaRefContParamSnp, uTrkDeltaRefContParamSnp, int8_t);
DECLARE_SOA_COLUMN(UTrkDeltaRefContParamTgl, uTrkDeltaRefContParamTgl, int8_t);
DECLARE_SOA_COLUMN(UTrkDeltaRefContParamQ2Pt, uTrkDeltaRefContParamQ2Pt, int8_t);
DECLARE_SOA_COLUMN(UTrkDeltaRefGloParamY, uTrkDeltaRefGloParamY, int8_t);
DECLARE_SOA_COLUMN(UTrkDeltaRefGloParamZ, uTrkDeltaRefGloParamZ, int8_t);
DECLARE_SOA_COLUMN(UTrkDeltaRefGloParamSnp, uTrkDeltaRefGloParamSnp, int8_t);
DECLARE_SOA_COLUMN(UTrkDeltaRefGloParamTgl, uTrkDeltaRefGloParamTgl, int8_t);
DECLARE_SOA_COLUMN(UTrkDeltaRefGloParamQ2Pt, uTrkDeltaRefGloParamQ2Pt, int8_t);
DECLARE_SOA_COLUMN(UTrkDeltaTOFdX, uTrkDeltaTOFdX, int8_t);
DECLARE_SOA_COLUMN(UTrkDeltaTOFdZ, uTrkDeltaTOFdZ, int8_t);

// aod::TrackQA version 003 columns for LTrk
DECLARE_SOA_COLUMN(LTrkDeltaRefContParamY, lTrkDeltaRefContParamY, int8_t);
DECLARE_SOA_COLUMN(LTrkDeltaRefContParamZ, lTrkDeltaRefContParamZ, int8_t);
DECLARE_SOA_COLUMN(LTrkDeltaRefContParamSnp, lTrkDeltaRefContParamSnp, int8_t);
DECLARE_SOA_COLUMN(LTrkDeltaRefContParamTgl, lTrkDeltaRefContParamTgl, int8_t);
DECLARE_SOA_COLUMN(LTrkDeltaRefContParamQ2Pt, lTrkDeltaRefContParamQ2Pt, int8_t);
DECLARE_SOA_COLUMN(LTrkDeltaRefGloParamY, lTrkDeltaRefGloParamY, int8_t);
DECLARE_SOA_COLUMN(LTrkDeltaRefGloParamZ, lTrkDeltaRefGloParamZ, int8_t);
DECLARE_SOA_COLUMN(LTrkDeltaRefGloParamSnp, lTrkDeltaRefGloParamSnp, int8_t);
DECLARE_SOA_COLUMN(LTrkDeltaRefGloParamTgl, lTrkDeltaRefGloParamTgl, int8_t);
DECLARE_SOA_COLUMN(LTrkDeltaRefGloParamQ2Pt, lTrkDeltaRefGloParamQ2Pt, int8_t);
DECLARE_SOA_COLUMN(LTrkDeltaTOFdX, lTrkDeltaTOFdX, int8_t);
DECLARE_SOA_COLUMN(LTrkDeltaTOFdZ, lTrkDeltaTOFdZ, int8_t);

// Occupancy realted information
DECLARE_SOA_COLUMN(UTrkTmoPrimUnfm80, uTrkTmoPrimUnfm80, float);
DECLARE_SOA_COLUMN(UTrkTmoRobustT0V0PrimUnfm80, uTrkTmoRobustT0V0PrimUnfm80, float);
DECLARE_SOA_COLUMN(UTrkTwmoPrimUnfm80, uTrkTwmoPrimUnfm80, float);
DECLARE_SOA_COLUMN(UTrkTwmoRobustT0V0PrimUnfm80, uTrkTwmoRobustT0V0PrimUnfm80, float);

DECLARE_SOA_COLUMN(LTrkTmoPrimUnfm80, lTrkTmoPrimUnfm80, float);
DECLARE_SOA_COLUMN(LTrkTmoRobustT0V0PrimUnfm80, lTrkTmoRobustT0V0PrimUnfm80, float);
DECLARE_SOA_COLUMN(LTrkTwmoPrimUnfm80, lTrkTwmoPrimUnfm80, float);
DECLARE_SOA_COLUMN(LTrkTwmoRobustT0V0PrimUnfm80, lTrkTwmoRobustT0V0PrimUnfm80, float);

// Additional Columns for Track information
DECLARE_SOA_COLUMN(UTrkHasColl, uTrkHasColl, bool);
DECLARE_SOA_COLUMN(UTrkVtxZ, uTrkVtxZ, float);

DECLARE_SOA_COLUMN(LTrkHasColl, lTrkHasColl, bool);
DECLARE_SOA_COLUMN(LTrkVtxZ, lTrkVtxZ, float);

// Cosmic Pair Information Columns
DECLARE_SOA_COLUMN(BFieldKG, bFieldKG, int8_t);

DECLARE_SOA_COLUMN(PtCut15, ptCut15, bool);
DECLARE_SOA_COLUMN(IsOKSum6, isOKSum6, bool);
DECLARE_SOA_COLUMN(IsOK6, isOK6, bool);
DECLARE_SOA_COLUMN(TOK0, tOK0, bool);
DECLARE_SOA_COLUMN(IsSelected, isSelected, bool);
DECLARE_SOA_COLUMN(DP0Corr, dP0Corr, float);
DECLARE_SOA_COLUMN(IsOK3Corr, isOK3Corr, bool);
DECLARE_SOA_COLUMN(IsOK3Clean, isOK3Clean, bool);

DECLARE_SOA_COLUMN(HasCosmicPair, hasCosmicPair, bool);
DECLARE_SOA_COLUMN(HasIsSelectedCosmic, hasIsSelectedCosmic, bool);
DECLARE_SOA_COLUMN(CosmicPairIndexList, cosmicPairIndexList, std::vector<int>);
}; // namespace csmu

// CosmicPairs00 Table: contains the comic selection flag information
DECLARE_SOA_TABLE(CosmicPairs00, "AOD", "COSMICPAIRS00", o2::soa::Index<>,
                  o2::aod::csmu::BFieldKG,
                  o2::aod::csmu::PtCut15,
                  o2::aod::csmu::IsOKSum6,
                  o2::aod::csmu::IsOK6,
                  o2::aod::csmu::TOK0,
                  o2::aod::csmu::IsSelected,

                  o2::aod::csmu::DP0Corr,
                  o2::aod::csmu::IsOK3Corr,
                  o2::aod::csmu::IsOK3Clean);

// CosmicPairs01 Table: contains the BC/TF info
DECLARE_SOA_TABLE(CosmicPairs01, "AOD", "COSMICPAIRS01",
                  o2::aod::csmu::UTrackId,
                  o2::aod::csmu::UTrkGlobalBcStart,
                  o2::aod::csmu::UTrkGlobalBcEnd,
                  // o2::aod::csmu::UTimestamp,
                  o2::aod::csmu::UTfId,
                  o2::aod::csmu::UBCinTF,

                  o2::aod::csmu::LTrackId,
                  o2::aod::csmu::LTrkGlobalBcStart,
                  o2::aod::csmu::LTrkGlobalBcEnd,
                  // o2::aod::csmu::LTimestamp,
                  o2::aod::csmu::LTfId,
                  o2::aod::csmu::LBCinTF);

// CosmicPairs02 Table: contains the track kinematics
DECLARE_SOA_TABLE(CosmicPairs02, "AOD", "COSMICPAIRS02",
                  o2::aod::csmu::UTrkX,
                  o2::aod::csmu::UTrkAlpha,
                  o2::aod::csmu::UTrkY,
                  o2::aod::csmu::UTrkZ,
                  o2::aod::csmu::UTrkSnp,
                  o2::aod::csmu::UTrkTgl,
                  o2::aod::csmu::UTrkSigned1Pt,
                  o2::aod::csmu::UTrkPt,
                  o2::aod::csmu::UTrkP,
                  o2::aod::csmu::UTrkEta,
                  o2::aod::csmu::UTrkPhi,

                  o2::aod::csmu::LTrkX,
                  o2::aod::csmu::LTrkAlpha,
                  o2::aod::csmu::LTrkY,
                  o2::aod::csmu::LTrkZ,
                  o2::aod::csmu::LTrkSnp,
                  o2::aod::csmu::LTrkTgl,
                  o2::aod::csmu::LTrkSigned1Pt,
                  o2::aod::csmu::LTrkPt,
                  o2::aod::csmu::LTrkP,
                  o2::aod::csmu::LTrkEta,
                  o2::aod::csmu::LTrkPhi);

// CosmicPairs03 Table: Track DCA information
DECLARE_SOA_TABLE(CosmicPairs03, "AOD", "COSMICPAIRS03",
                  o2::aod::csmu::UTrkDcaCalc,
                  o2::aod::csmu::UTrkDcaXY,
                  o2::aod::csmu::UTrkDcaZ,
                  o2::aod::csmu::UTrkSigmaDcaXY2,
                  o2::aod::csmu::UTrkSigmaDcaZ2,

                  o2::aod::csmu::LTrkDcaCalc,
                  o2::aod::csmu::LTrkDcaXY,
                  o2::aod::csmu::LTrkDcaZ,
                  o2::aod::csmu::LTrkSigmaDcaXY2,
                  o2::aod::csmu::LTrkSigmaDcaZ2);

DECLARE_SOA_TABLE(CosmicPairs04, "AOD", "COSMICPAIRS04",
                  o2::aod::csmu::UTrkHasColl,
                  o2::aod::csmu::UTrkVtxZ,
                  o2::aod::csmu::LTrkHasColl,
                  o2::aod::csmu::LTrkVtxZ);

// CosmicPairs05 Table: full UTrk and LTrk extra track information
DECLARE_SOA_TABLE(CosmicPairs05, "AOD", "COSMICPAIRS05",
                  o2::aod::csmu::UTrkFlags,
                  o2::aod::csmu::UTrkLength,
                  o2::aod::csmu::UTrkPIDForTracking,
                  o2::aod::csmu::UTrkIsPVContributor,
                  o2::aod::csmu::UTrkTrackTime,
                  o2::aod::csmu::UTrkTrackTimeRes,
                  o2::aod::csmu::UTrkDetectorMap,

                  o2::aod::csmu::LTrkFlags,
                  o2::aod::csmu::LTrkLength,
                  o2::aod::csmu::LTrkPIDForTracking,
                  o2::aod::csmu::LTrkIsPVContributor,
                  o2::aod::csmu::LTrkTrackTime,
                  o2::aod::csmu::LTrkTrackTimeRes,
                  o2::aod::csmu::LTrkDetectorMap);

// CosmicPairs06 Table: Contains Track Extra TPC information
DECLARE_SOA_TABLE(CosmicPairs06, "AOD", "COSMICPAIRS06",
                  o2::aod::csmu::UTrkTPCInnerParam,
                  o2::aod::csmu::UTrkTPCNClsFindable,
                  o2::aod::csmu::UTrkTPCNClsFindableMinusFound,
                  o2::aod::csmu::UTrkTPCNClsFindableMinusPID,
                  o2::aod::csmu::UTrkTPCNClsFindableMinusCrossedRows,
                  o2::aod::csmu::UTrkTPCNClsShared,
                  o2::aod::csmu::UTrkTPCDeltaTFwd,
                  o2::aod::csmu::UTrkTPCDeltaTBwd,
                  o2::aod::csmu::UTrkTPCChi2NCl,
                  o2::aod::csmu::UTrkTPCSignal,
                  o2::aod::csmu::UTrkHasTPC,
                  o2::aod::csmu::UTrkTPCNClsFound,
                  o2::aod::csmu::UTrkTPCNClsPID,
                  o2::aod::csmu::UTrkTPCNClsCrossedRows,
                  o2::aod::csmu::UTrkTPCCrossedRowsOverFindableCls,
                  o2::aod::csmu::UTrkTPCFoundOverFindableCls,
                  o2::aod::csmu::UTrkTPCFractionSharedCls,

                  o2::aod::csmu::LTrkTPCInnerParam,
                  o2::aod::csmu::LTrkTPCNClsFindable,
                  o2::aod::csmu::LTrkTPCNClsFindableMinusFound,
                  o2::aod::csmu::LTrkTPCNClsFindableMinusPID,
                  o2::aod::csmu::LTrkTPCNClsFindableMinusCrossedRows,
                  o2::aod::csmu::LTrkTPCNClsShared,
                  o2::aod::csmu::LTrkTPCDeltaTFwd,
                  o2::aod::csmu::LTrkTPCDeltaTBwd,
                  o2::aod::csmu::LTrkTPCChi2NCl,
                  o2::aod::csmu::LTrkTPCSignal,
                  o2::aod::csmu::LTrkHasTPC,
                  o2::aod::csmu::LTrkTPCNClsFound,
                  o2::aod::csmu::LTrkTPCNClsPID,
                  o2::aod::csmu::LTrkTPCNClsCrossedRows,
                  o2::aod::csmu::LTrkTPCCrossedRowsOverFindableCls,
                  o2::aod::csmu::LTrkTPCFoundOverFindableCls,
                  o2::aod::csmu::LTrkTPCFractionSharedCls);

// CosmicPairs07 Table: Contains Track Extra ITS information
DECLARE_SOA_TABLE(CosmicPairs07, "AOD", "COSMICPAIRS07",
                  o2::aod::csmu::UTrkITSClusterSizes,
                  o2::aod::csmu::UTrkITSChi2NCl,
                  o2::aod::csmu::UTrkHasITS,
                  o2::aod::csmu::UTrkITSClusterMap,
                  o2::aod::csmu::UTrkITSNCls,
                  o2::aod::csmu::UTrkITSNClsInnerBarrel,
                  // o2::aod::csmu::UTrkITSClsSizeInLayer,
                  o2::aod::csmu::UTrkIsITSAfterburner,

                  o2::aod::csmu::LTrkITSClusterSizes,
                  o2::aod::csmu::LTrkITSChi2NCl,
                  o2::aod::csmu::LTrkHasITS,
                  o2::aod::csmu::LTrkITSClusterMap,
                  o2::aod::csmu::LTrkITSNCls,
                  o2::aod::csmu::LTrkITSNClsInnerBarrel,
                  // o2::aod::csmu::LTrkITSClsSizeInLayer,
                  o2::aod::csmu::LTrkIsITSAfterburner);

// CosmicPairs08 Table: Contains Track Extra TRD information
DECLARE_SOA_TABLE(CosmicPairs08, "AOD", "COSMICPAIRS08",
                  o2::aod::csmu::UTrkTRDPattern,
                  o2::aod::csmu::UTrkTRDChi2,
                  o2::aod::csmu::UTrkTRDSignal,
                  o2::aod::csmu::UTrkHasTRD,

                  o2::aod::csmu::LTrkTRDPattern,
                  o2::aod::csmu::LTrkTRDChi2,
                  o2::aod::csmu::LTrkTRDSignal,
                  o2::aod::csmu::LTrkHasTRD);

// CosmicPairs09 Table: Contains Track Extra TOF information
DECLARE_SOA_TABLE(CosmicPairs09, "AOD", "COSMICPAIRS09",
                  o2::aod::csmu::UTrkTOFChi2,
                  o2::aod::csmu::UTrkTOFExpMom,
                  o2::aod::csmu::UTrkHasTOF,
                  o2::aod::csmu::UTrkTOFExpTimeEl,
                  o2::aod::csmu::UTrkTOFExpTimeMu,
                  o2::aod::csmu::UTrkTOFExpTimePi,
                  o2::aod::csmu::UTrkTOFExpTimeKa,
                  o2::aod::csmu::UTrkTOFExpTimePr,
                  o2::aod::csmu::UTrkTOFExpTimeDe,
                  o2::aod::csmu::UTrkTOFExpTimeTr,
                  o2::aod::csmu::UTrkTOFExpTimeHe,
                  o2::aod::csmu::UTrkTOFExpTimeAl,

                  o2::aod::csmu::LTrkTOFChi2,
                  o2::aod::csmu::LTrkTOFExpMom,
                  o2::aod::csmu::LTrkHasTOF,
                  o2::aod::csmu::LTrkTOFExpTimeEl,
                  o2::aod::csmu::LTrkTOFExpTimeMu,
                  o2::aod::csmu::LTrkTOFExpTimePi,
                  o2::aod::csmu::LTrkTOFExpTimeKa,
                  o2::aod::csmu::LTrkTOFExpTimePr,
                  o2::aod::csmu::LTrkTOFExpTimeDe,
                  o2::aod::csmu::LTrkTOFExpTimeTr,
                  o2::aod::csmu::LTrkTOFExpTimeHe,
                  o2::aod::csmu::LTrkTOFExpTimeAl);

// CosmicPairs10 Table: Contains Track Extra EMCal information
DECLARE_SOA_TABLE(CosmicPairs10, "AOD", "COSMICPAIRS10",
                  o2::aod::csmu::UTrkTrackEtaEmcal,
                  o2::aod::csmu::UTrkTrackPhiEmcal,

                  o2::aod::csmu::LTrkTrackEtaEmcal,
                  o2::aod::csmu::LTrkTrackPhiEmcal);

// aod::pidTOFFlags && aod::TOFSignal information
DECLARE_SOA_TABLE(CosmicPairs11, "AOD", "COSMICPAIRS11",
                  o2::aod::csmu::UTrkGoodTOFMatch,
                  o2::aod::csmu::UTrkTOFSignal,
                  // o2::aod::csmu::UTrkEventCollisionTime,

                  o2::aod::csmu::LTrkGoodTOFMatch,
                  o2::aod::csmu::LTrkTOFSignal //,
                                               // o2::aod::csmu::LTrkEventCollisionTime
);

DECLARE_SOA_TABLE(CosmicPairs12, "AOD", "COSMICPAIRS12",
                  o2::aod::csmu::UTrkTOFFlags,
                  o2::aod::csmu::UTrkIsEvTimeDefined,
                  o2::aod::csmu::UTrkIsEvTimeTOF,
                  o2::aod::csmu::UTrkIsEvTimeT0AC,
                  o2::aod::csmu::UTrkIsEvTimeTOFT0AC,

                  o2::aod::csmu::LTrkTOFFlags,
                  o2::aod::csmu::LTrkIsEvTimeDefined,
                  o2::aod::csmu::LTrkIsEvTimeTOF,
                  o2::aod::csmu::LTrkIsEvTimeT0AC,
                  o2::aod::csmu::LTrkIsEvTimeTOFT0AC);

DECLARE_SOA_TABLE(CosmicPairs13, "AOD", "COSMICPAIRS13",
                  o2::aod::csmu::UTrkTOFEvTime,
                  o2::aod::csmu::UTrkTOFEvTimeErr,

                  o2::aod::csmu::LTrkTOFEvTime,
                  o2::aod::csmu::LTrkTOFEvTimeErr);

DECLARE_SOA_TABLE(CosmicPairs14, "AOD", "COSMICPAIRS14",
                  o2::aod::csmu::UTrkTPCTime0,
                  o2::aod::csmu::UTrkTPCdEdxNorm,
                  o2::aod::csmu::UTrkTPCDCAR,
                  o2::aod::csmu::UTrkTPCDCAZ,
                  o2::aod::csmu::UTrkTPCClusterByteMask,
                  o2::aod::csmu::UTrkTPCdEdxMax0R,
                  o2::aod::csmu::UTrkTPCdEdxMax1R,
                  o2::aod::csmu::UTrkTPCdEdxMax2R,
                  o2::aod::csmu::UTrkTPCdEdxMax3R,
                  o2::aod::csmu::UTrkTPCdEdxTot0R,
                  o2::aod::csmu::UTrkTPCdEdxTot1R,
                  o2::aod::csmu::UTrkTPCdEdxTot2R,
                  o2::aod::csmu::UTrkTPCdEdxTot3R,

                  o2::aod::csmu::LTrkTPCTime0,
                  o2::aod::csmu::LTrkTPCdEdxNorm,
                  o2::aod::csmu::LTrkTPCDCAR,
                  o2::aod::csmu::LTrkTPCDCAZ,
                  o2::aod::csmu::LTrkTPCClusterByteMask,
                  o2::aod::csmu::LTrkTPCdEdxMax0R,
                  o2::aod::csmu::LTrkTPCdEdxMax1R,
                  o2::aod::csmu::LTrkTPCdEdxMax2R,
                  o2::aod::csmu::LTrkTPCdEdxMax3R,
                  o2::aod::csmu::LTrkTPCdEdxTot0R,
                  o2::aod::csmu::LTrkTPCdEdxTot1R,
                  o2::aod::csmu::LTrkTPCdEdxTot2R,
                  o2::aod::csmu::LTrkTPCdEdxTot3R);

DECLARE_SOA_TABLE(CosmicPairs15, "AOD", "COSMICPAIRS15",
                  o2::aod::csmu::UTrkDeltaRefContParamY,
                  o2::aod::csmu::UTrkDeltaRefContParamZ,
                  o2::aod::csmu::UTrkDeltaRefContParamSnp,
                  o2::aod::csmu::UTrkDeltaRefContParamTgl,
                  o2::aod::csmu::UTrkDeltaRefContParamQ2Pt,
                  o2::aod::csmu::UTrkDeltaRefGloParamY,
                  o2::aod::csmu::UTrkDeltaRefGloParamZ,
                  o2::aod::csmu::UTrkDeltaRefGloParamSnp,
                  o2::aod::csmu::UTrkDeltaRefGloParamTgl,
                  o2::aod::csmu::UTrkDeltaRefGloParamQ2Pt,
                  o2::aod::csmu::UTrkDeltaTOFdX,
                  o2::aod::csmu::UTrkDeltaTOFdZ,

                  o2::aod::csmu::LTrkDeltaRefContParamY,
                  o2::aod::csmu::LTrkDeltaRefContParamZ,
                  o2::aod::csmu::LTrkDeltaRefContParamSnp,
                  o2::aod::csmu::LTrkDeltaRefContParamTgl,
                  o2::aod::csmu::LTrkDeltaRefContParamQ2Pt,
                  o2::aod::csmu::LTrkDeltaRefGloParamY,
                  o2::aod::csmu::LTrkDeltaRefGloParamZ,
                  o2::aod::csmu::LTrkDeltaRefGloParamSnp,
                  o2::aod::csmu::LTrkDeltaRefGloParamTgl,
                  o2::aod::csmu::LTrkDeltaRefGloParamQ2Pt,
                  o2::aod::csmu::LTrkDeltaTOFdX,
                  o2::aod::csmu::LTrkDeltaTOFdZ);

DECLARE_SOA_TABLE(CosmicPairs16, "AOD", "COSMICPAIRS16",
                  o2::aod::csmu::UTrkTmoPrimUnfm80,
                  o2::aod::csmu::UTrkTmoRobustT0V0PrimUnfm80,
                  o2::aod::csmu::UTrkTwmoPrimUnfm80,
                  o2::aod::csmu::UTrkTwmoRobustT0V0PrimUnfm80,

                  o2::aod::csmu::LTrkTmoPrimUnfm80,
                  o2::aod::csmu::LTrkTmoRobustT0V0PrimUnfm80,
                  o2::aod::csmu::LTrkTwmoPrimUnfm80,
                  o2::aod::csmu::LTrkTwmoRobustT0V0PrimUnfm80);

DECLARE_SOA_TABLE(CollCosmicFlags, "AOD", "COLLCOSMICFLAGS", o2::soa::Index<>,
                  o2::aod::csmu::HasCosmicPair, o2::aod::csmu::HasIsSelectedCosmic);
DECLARE_SOA_TABLE(CollToCosmic, "AOD", "COLLTOCOSMIC", o2::aod::csmu::CosmicPairIndexList);

} // namespace aod
} // namespace o2

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace o2::constants::math;

template <typename T>
void printTime(T Start, std::string String)
{
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = duration_cast<std::chrono::microseconds>(stop - Start);
  LOG(info) << String << float(duration.count()) / float(1000000) << " seconds"; //<<endl;
}

int32_t nBCsPerOrbit = o2::constants::lhc::LHCMaxBunches;

struct CosmicMuonToCollisionAssociator {

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // declare production of tables
  Produces<aod::CosmicPairs00> genCosmicPairs00;
  Produces<aod::CosmicPairs01> genCosmicPairs01;
  Produces<aod::CosmicPairs02> genCosmicPairs02;
  Produces<aod::CosmicPairs03> genCosmicPairs03;
  Produces<aod::CosmicPairs04> genCosmicPairs04;
  Produces<aod::CosmicPairs05> genCosmicPairs05;
  Produces<aod::CosmicPairs06> genCosmicPairs06;
  Produces<aod::CosmicPairs07> genCosmicPairs07;
  Produces<aod::CosmicPairs08> genCosmicPairs08;
  Produces<aod::CosmicPairs09> genCosmicPairs09;
  Produces<aod::CosmicPairs10> genCosmicPairs10;
  Produces<aod::CosmicPairs11> genCosmicPairs11;
  Produces<aod::CosmicPairs12> genCosmicPairs12;
  Produces<aod::CosmicPairs13> genCosmicPairs13;
  Produces<aod::CosmicPairs14> genCosmicPairs14;
  Produces<aod::CosmicPairs15> genCosmicPairs15;
  Produces<aod::CosmicPairs16> genCosmicPairs16;

  Produces<aod::CollCosmicFlags> genCollisionCosmicFlags;
  Produces<aod::CollToCosmic> genCollToCosmic;

  // Histogram registry;
  HistogramRegistry histReg{"histReg", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  // Configurables

  Configurable<int> customOrbitOffset{"customOrbitOffset", 0, "customOrbitOffset for MC"};
  Configurable<int> nBCinTF{"nBCinTF", 114048, "nBCinTF"};

  struct : ConfigurableGroup {
    Configurable<float> csmMuThrPt{"csmMuThrPt", 0.5, "csmMuThrPt"};
    Configurable<float> csmMuThrDCAMin{"csmMuThrDCAMin", 0.5, "csmMuThrDCAMin"};
    Configurable<float> csmMuThrDCAMax{"csmMuThrDCAMax", 120, "csmMuThrDCAMax"};
    Configurable<int> csmMuTpcNClsFound{"csmMuTpcNClsFound", 100, "csmMuTpcNClsFound"};

    Configurable<float> csmMuSumPtPair{"csmMuSumPtPair", 2.0, "csmMuSumPtPair"};
    Configurable<float> csmMuSumQPtPair{"csmMuSumQPtPair", 0.2, "csmMuSumQPtPair"};
    Configurable<float> csmMuSumTglPair{"csmMuSumTglPair", 0.01, "csmMuSumTglPair"};
    Configurable<float> csmMuSumDcaXY{"csmMuSumDcaXY", 6.0, "csmMuSumDcaXY"};
    Configurable<float> csmMuDiffDcaXY{"csmMuDiffDcaXY", 0.6, "csmMuDiffDcaXY"};
    Configurable<float> csmMuDiffAlphaPair{"csmMuDiffAlphaPair", 0.1, "csmMuDiffAlphaPair"};
    Configurable<int> deltaBC{"deltaBC", 10000, "deltaBC"};
    Configurable<int> deltaBCForColl{"deltaBCForColl", 4000, "deltaBCForColl"};

    Configurable<float> csmMuDiffQPtPair{"csmMuDiffQPtPair", 0.2, "csmMuDiffQPtPair"};
    Configurable<float> csmMuDiffTglPair{"csmMuDiffTglPair", 0.2, "csmMuDiffTglPair"};
    Configurable<float> csmMuDiffDcaCalcPair{"csmMuDiffDcaCalcPair", 5, "csmMuDiffDcaCalcPair"};

    Configurable<bool> produceTable00{"produceTable00", true, "Produce CosmicPairs00: pair selection flags"};
    Configurable<bool> produceTable01{"produceTable01", true, "Produce CosmicPairs01: track IDs and BC/TF info"};
    Configurable<bool> produceTable02{"produceTable02", true, "Produce CosmicPairs02: track kinematics"};
    Configurable<bool> produceTable03{"produceTable03", true, "Produce CosmicPairs03: track DCA info"};
    Configurable<bool> produceTable04{"produceTable04", true, "Produce CosmicPairs04: track HasColl and VtxZ"};
    Configurable<bool> produceTable05{"produceTable05", false, "Produce CosmicPairs05: track flags, length, PID, time"};
    Configurable<bool> produceTable06{"produceTable06", false, "Produce CosmicPairs06: TPC cluster and signal info"};
    Configurable<bool> produceTable07{"produceTable07", false, "Produce CosmicPairs07: ITS cluster info"};
    Configurable<bool> produceTable08{"produceTable08", false, "Produce CosmicPairs08: TRD pattern, chi2, signal"};
    Configurable<bool> produceTable09{"produceTable09", false, "Produce CosmicPairs09: TOF expected times"};
    Configurable<bool> produceTable10{"produceTable10", false, "Produce CosmicPairs10: EMCal eta and phi"};
    Configurable<bool> produceTable11{"produceTable11", false, "Produce CosmicPairs11: TOF signal and good match"};
    Configurable<bool> produceTable12{"produceTable12", false, "Produce CosmicPairs12: TOF event time flags"};
    Configurable<bool> produceTable13{"produceTable13", false, "Produce CosmicPairs13: TOF event time and error"};
    Configurable<bool> produceTable14{"produceTable14", false, "Produce CosmicPairs14: TrackQA TPC dEdx and DCA"};
    Configurable<bool> produceTable15{"produceTable15", false, "Produce CosmicPairs15: TrackQA delta ref parameters"};
    Configurable<bool> produceTable16{"produceTable16", false, "Produce CosmicPairs16: track mean occupancies"};

    Configurable<bool> csmMuDoBkgDownSampling{"csmMuDoBkgDownSampling", true, "csmMuDoBkgDownSampling"};
    Configurable<float> csmMuBkgDownSampling{"csmMuBkgDownSampling", 0.01, "csmMuBkgDownSampling"};

    Configurable<bool> csmMuCheckSumDcaXY{"csmMuCheckSumDcaXY", true, "csmMuCheckSumDcaXY"};
    Configurable<bool> csmMuCheckDiffDcaXY{"csmMuCheckDiffDcaXY", true, "csmMuCheckDiffDcaXY"};
    Configurable<bool> csmMuCheckDiffAlphaPair{"csmMuCheckDiffAlphaPair", true, "csmMuCheckDiffAlphaPair"};
  } cfgCM;

  struct : ConfigurableGroup {
    Configurable<bool> printDebugMessages{"printDebugMessages", true, "printDebugMessages"};
    Configurable<bool> resetHistograms{"resetHistograms", false, "resetHistograms"};
  } cfgDebug;

  Configurable<std::vector<double>> countBins{"countBins",
                                              {// 0 to 10 (step 1)
                                               0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                                               // 10 to 100 (step 10)
                                               20, 30, 40, 50, 60, 70, 80, 90, 100,
                                               // 100 to 1,000 (step 100)
                                               200, 300, 400, 500, 600, 700, 800, 900, 1000,
                                               // 1,000 to 10,000 (step 1,000)
                                               2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000,
                                               // 10,000 to 100,000 (step 10,000)
                                               20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000,
                                               // 100,000 to 1,000,000 (step 100,000)
                                               200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000,
                                               // 1,000,000 to 10,000,000 (step 1,000,000)
                                               2000000, 3000000, 4000000, 5000000, 6000000, 7000000, 8000000, 9000000, 10000000,
                                               // 10,000,000 to 100,000,000 (step 10,000,000)
                                               20000000, 30000000, 40000000, 50000000, 60000000, 70000000, 80000000, 90000000, 100000000,
                                               // 100,000,000 to 1,000,000,000 (step 100,000,000)
                                               200000000, 300000000, 400000000, 500000000, 600000000, 700000000, 800000000, 900000000, 1000000000},
                                              "Particle Count Bins"};

  enum CosmicPairRejectionType {
    kPairPassed = 0,
    kFailSumPtPair,
    kFailSumQPtPair,
    kFailSumTglPair,
    kFailSumDcaXY,
    kFailDiffDcaXY,
    kFailDiffAlphaPair
  };

  static constexpr std::string_view CMRejectionTag[]{
    "kPairPassed",
    "kFailSumPtPair",
    "kFailSumQPtPair",
    "kFailSumTglPair",
    "kFailSumDcaXY",
    "kFailDiffDcaXY",
    "kFailDiffAlphaPair"};

  void init(InitContext const&)
  {
    const AxisSpec axisCount{countBins, "Counts"};

    histReg.add("hTiming_nDFsProcessed", "hTiming_nDFsProcessed", kTH1F, {{1, 0, 1}});
    histReg.add("hTiming_nTFsPerDataFrame", "hTiming_nTFsPerDataFrame", kTH1F, {{5000, 0, 5000}});
    histReg.add("hTiming_nTFsProcessed", "hTiming_nTFsProcessed", kTH1F, {{1, 0, 1}});

    histReg.add("hCosmicPairs_nCosmicsPerDF", "hCosmicPairs_nCosmicsPerDF", kTH1F, {{5000, 0, 5000}});
    histReg.add("hCosmicPairs_nIsSelCosmicsPerDF", "hCosmicPairs_nIsSelCosmicsPerDF", kTH1F, {{5000, 0, 5000}});

    histReg.add("hRatio", "Ratio;Ratio;Entries", kTH1F, {{1000, 0, 1}});
    histReg.add("hRatioVtx8", "Ratio Vtx<8;Ratio;Entries", kTH1F, {{1000, 0, 1}});
    histReg.add("hRatioVtx10", "Ratio Vtx<10;Ratio;Entries", kTH1F, {{1000, 0, 1}});

    histReg.add("hRatioIsSel", "Ratio IsSelected;Ratio;Entries", kTH1F, {{1000, 0, 1}});
    histReg.add("hRatioIsSelVtx8", "Ratio IsSelected Vtx<8;Ratio;Entries", kTH1F, {{1000, 0, 1}});
    histReg.add("hRatioIsSelVtx10", "Ratio IsSelected Vtx<10;Ratio;Entries", kTH1F, {{1000, 0, 1}});

    histReg.add("hCosmicFlag", "Cosmic Flag;Count;Entries", kTH1F, {axisCount});
    histReg.add("hCosmicFlagVtx8", "Cosmic Flag Vtx<8;Count;Entries", kTH1F, {axisCount});
    histReg.add("hCosmicFlagVtx10", "Cosmic Flag Vtx<10;Count;Entries", kTH1F, {axisCount});

    histReg.add("hCosmicFlagIsSel", "IsSelected Cosmic Flag;Count;Entries", kTH1F, {axisCount});
    histReg.add("hCosmicFlagIsSelVtx8", "IsSelected Cosmic Flag Vtx<8;Count;Entries", kTH1F, {axisCount});
    histReg.add("hCosmicFlagIsSelVtx10", "IsSelected Cosmic Flag Vtx<10;Count;Entries", kTH1F, {axisCount});

    histReg.add("hColls", "Total Colls;Count;Entries", kTH1F, {axisCount});
    histReg.add("hCollsVtx8", "Colls Vtx<8;Count;Entries", kTH1F, {axisCount});
    histReg.add("hCollsVtx10", "Colls Vtx<10;Count;Entries", kTH1F, {axisCount});
  }

  template <typename T, std::size_t N>
  void sortVectorOfArray(std::vector<std::array<T, N>>& myVector, const int& myIDX)
  {
    std::stable_sort(myVector.begin(), myVector.end(), [myIDX](const std::array<T, N>& a, const std::array<T, N>& b) {
      return a[myIDX] < b[myIDX]; // sort at the required index //stable sort keeps the order as is, if index are equal.
    });
  }

  /// @brief Find first index where array[idx][keyIdx] >= targetValue
  /// @return Index of first element >= target, or array.size() if none found
  template <typename T, std::size_t N>
  int binarySearchLower(const std::vector<std::array<T, N>>& arr, int keyIdx, T targetValue)
  {
    int lo = 0;
    int hi = static_cast<int>(arr.size());
    while (lo < hi) {
      int mid = lo + ((hi - lo) >> 1);
      if (arr[mid][keyIdx] < targetValue) {
        lo = mid + 1;
      } else {
        hi = mid;
      }
    }
    return lo;
  }

  /// @brief Find first index where array[idx][keyIdx] > targetValue
  /// @return Index of first element > target, or array.size() if none found
  template <typename T, std::size_t N>
  int binarySearchUpper(const std::vector<std::array<T, N>>& arr, int keyIdx, T targetValue)
  {
    int lo = 0;
    int hi = static_cast<int>(arr.size());
    while (lo < hi) {
      int mid = lo + ((hi - lo) >> 1);
      if (arr[mid][keyIdx] <= targetValue) {
        lo = mid + 1;
      } else {
        hi = mid;
      }
    }
    return lo;
  }

  void getRunInfo(const int& run, int& nBCsPerTF, int64_t& bcSOR)
  {
    auto runDuration = ccdb->getRunDuration(run, true);
    int64_t tsSOR = runDuration.first;
    auto ctpx = ccdb->getForTimeStamp<std::vector<int64_t>>("CTP/Calib/OrbitReset", tsSOR);
    int64_t tsOrbitReset = (*ctpx)[0];
    uint32_t nOrbitsPerTF = run < 534133 ? 128 : 32; // o2-linter: disable=magic-number (algorithm constant)
    int64_t orbitSOR = (tsSOR * 1000 - tsOrbitReset) / o2::constants::lhc::LHCOrbitMUS;
    orbitSOR = orbitSOR / nOrbitsPerTF * nOrbitsPerTF;
    bcSOR = orbitSOR * nBCsPerOrbit + customOrbitOffset * nBCsPerOrbit; // customOrbitOffset is a Configurable
    nBCsPerTF = nOrbitsPerTF * nBCsPerOrbit;
  }

  template <typename T>
  void getTimingInfo(const T& bc, int& lastRun, int32_t& nBCsPerTF, int64_t& bcSOR, uint64_t& time, int64_t& tfIdThis, int64_t& bcInTF)
  {
    int run = bc.runNumber();
    if (run != lastRun) { // update run info
      lastRun = run;
      getRunInfo(run, nBCsPerTF, bcSOR); // update nBCsPerTF && bcSOR
    }
    // update the information
    time = bc.timestamp();
    tfIdThis = (bc.globalBC() - bcSOR) / nBCsPerTF;
    bcInTF = (bc.globalBC() - bcSOR) % nBCsPerTF;
  }

  int runNumber = -1;
  o2::parameters::GRPMagField* grpmag = nullptr;
  o2::base::MatLayerCylSet* lut = nullptr;
  int8_t bFieldKG;

  struct : ConfigurableGroup {
    Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
    Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
    Configurable<std::string> mVtxPath{"mVtxPath", "GLO/Calib/MeanVertex", "Path of the mean vertex file"};
  } cfgB;

  template <typename T>
  void initCCDB(const T& bc)
  {
    if (runNumber == bc.runNumber()) {
      return;
    }

    // load matLUT for this timestamp
    if (!lut) {
      LOG(info) << "Loading material look-up table for timestamp: " << bc.timestamp();
      lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->getForTimeStamp<o2::base::MatLayerCylSet>(cfgB.lutPath, bc.timestamp()));
    } else {
      LOG(info) << "Material look-up table already in place. Not reloading.";
    }

    grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(cfgB.grpmagPath, bc.timestamp());
    bFieldKG = grpmag->getNominalL3Field();
    runNumber = bc.runNumber();
  }

  static constexpr double kB2C = 0.299792458;   // o2-linter: disable=name/constexpr-constant (algorithm constant)
  static constexpr double kAlmost0Field = 1e-4; // o2-linter: disable=name/constexpr-constant (algorithm constant)
  double getCurvature(double signed1Pt, double b) { return signed1Pt * b * kB2C; }
  double getDCA(double trackX, double trackY, double alpha, double snp, double signed1Pt,
                double xRef, double yRef, double b)
  {
    // Handle zero or near-zero magnetic field by simple linear distance
    if (std::abs(b) < kAlmost0Field) {
      return std::hypot(trackX - xRef, trackY - yRef);
    }

    double curvature = getCurvature(signed1Pt, b);

    // Rotate reference point into track frame (rotation by -alpha)
    double sinAlpha = std::sin(alpha);
    double cosAlpha = std::cos(alpha);

    double xRot = xRef * cosAlpha + yRef * sinAlpha;
    double yRot = -xRef * sinAlpha + yRef * cosAlpha;

    // Compute difference vector in local track frame
    double dx = trackX - xRot;
    double dy = trackY - yRot;

    double sqrtTerm = std::sqrt((1.0 - snp) * (1.0 + snp));

    double sn = curvature * dx - snp;
    double cs = curvature * dy + sqrtTerm;

    double numerator = 2.0 * (dx * snp - dy * sqrtTerm) - curvature * (dx * dx + dy * dy);
    double denominator = 1.0 + std::sqrt(sn * sn + cs * cs);

    return -numerator / denominator;
  }

  template <typename T>
  std::vector<int> getCosmicTriggerList(int64_t globalBC, const T& list)
  {
    std::vector<int> triggeredCosmicIndexList;
    for (uint i = 0; i < list.size(); i++) {
      if (std::abs(globalBC - list[i]) < cfgCM.deltaBCForColl) {
        triggeredCosmicIndexList.push_back(i); // it had a cosmic entry;
      }
    }
    return triggeredCosmicIndexList;
  }

  enum VariableType {
    kTrkCountIdx,
    kTrkGI,       // track.globalIndex()
    kTrkFirstGBC, // first global bunch crossing of Track
    kTrkLastGBC,  // second global bunch crossing of Track
    kTrkTfIdThis, // Timeframe if of the Track
    kTrkBCinTF    // Bunch crossing of track in time frame
  };

  enum ProcessFunctionMode {
    kProcessBasic,
    kProcessWithTOFSignal,
    kProcessWithTOFEvTime,
    kProcessWithPidEvTimeFlags,
    kProcessWithTracksQA,
    kProcessWithOccupancies,
    kProcessFull
  };

  enum CosmicDaughterTrkType {
    kUpperTrk = 0,
    kLowerTrk
  };

  enum FloatVariableType {
    kTrkCalcDcaInfo,
    kTrkCollPosZ
  };

  int dfCount = 0;
  int debugCounter = 0;
  int32_t nBCsPerTF = -999;
  int64_t bcSOR = -999;
  uint64_t time = -1;
  int64_t tfIdThis = -1;
  int64_t bcInTF = -1;
  int lastRun = -999;

  std::mt19937 rng{123456789}; // Mersenne Twister generator with fixed random seed
  std::uniform_real_distribution<float> dist{0.f, 1.f};

  std::chrono::high_resolution_clock::time_point start0 = std::chrono::high_resolution_clock::now();

  template <int mode, typename B, typename C, typename T, typename A, typename Q, typename O>
  void executeProcess(const B& BCs, const C& colls, const T& tracks, const A& ambgTracks, const Q& tracksQA, const O& trackMeanOccs, const auto& Origins)
  {
    (void)tracksQA;
    (void)trackMeanOccs;

    dfCount++;
    auto start1 = std::chrono::high_resolution_clock::now();
    if (cfgDebug.printDebugMessages) {
      LOG(info) << "DEBUG :: df_" << dfCount << " :: DF_" << Origins.iteratorAt(0).dataframeID()
                << " :: BCs.size() = " << BCs.size()
                << " :: colls.size() = " << colls.size()
                << " :: tracks.size() = " << tracks.size()
                << " :: ambgTracks.size() = " << ambgTracks.size();
    }

    if (colls.size() == 0)
      return;

    // Step 1: loop over tracks and get good tracks.
    float trackDcaCalc;
    int64_t lastCollId = -999;

    std::vector<int64_t> trackRunNumber;
    std::vector<int64_t> trackGlobalBC;
    std::vector<int64_t> trackTriggerMask;
    std::vector<int64_t> trackTimestamp;
    std::vector<int64_t> trackTFidThis;
    std::vector<int64_t> trackBcInTF;
    std::vector<int64_t> fullTimeFrameIdList;

    std::vector<std::array<int64_t, 6>> upperTrkGIListWithGlobalBCInfo;
    std::vector<std::array<int64_t, 6>> lowerTrkGIListWithGlobalBCInfo;

    std::vector<std::array<float, 2>> upperTrkFloatInfo;
    std::vector<std::array<float, 2>> lowerTrkFloatInfo;

    std::vector<bool> upperTrkBoolInfo;
    std::vector<bool> lowerTrkBoolInfo;

    initCCDB(BCs.begin());
    auto coll = colls.begin();
    auto bc = BCs.begin();
    int64_t firstGlobalBC = 0;
    int64_t lastGlobalBC = 0;
    auto trackLoopStart = std::chrono::high_resolution_clock::now();
    int lTrkCounter = -1;
    int uTrkCounter = -1;
    bool trkHasCollision = false;
    float trkCollVtxZ = 0;

    for (const auto& track : tracks) {
      if (track.pt() < cfgCM.csmMuThrPt) {
        continue;
      }
      if (track.tpcNClsFound() < cfgCM.csmMuTpcNClsFound) {
        continue;
      } // Dynamic Column : filter can't be used here
      trackDcaCalc = getDCA(track.x(), track.y(), track.alpha(), track.snp(), track.signed1Pt(), 0, 0, bFieldKG);
      if (std::abs(trackDcaCalc) < cfgCM.csmMuThrDCAMin) {
        continue;
      }
      if (std::abs(trackDcaCalc) > cfgCM.csmMuThrDCAMax) {
        continue;
      }

      // Get Timing info of the track
      if (track.collisionId() < 0) {
        trkHasCollision = false;
        trkCollVtxZ = 0;

        lastCollId = -999;
        trackRunNumber.clear();
        trackGlobalBC.clear();
        trackTriggerMask.clear();
        trackTimestamp.clear();
        trackTFidThis.clear();
        trackBcInTF.clear();
        auto bcs = track.template ambgTrack_as<A>().template bc_as<B>();
        if (bcs.size() == 0) {
          trackRunNumber.push_back(0);
          trackGlobalBC.push_back(0);
          trackTriggerMask.push_back(0);
          trackTimestamp.push_back(0);
          trackTFidThis.push_back(-1);
          trackBcInTF.push_back(-1);
          firstGlobalBC = -10000;
          lastGlobalBC = -10000;
        } else {
          for (const auto& bc : bcs) {
            getTimingInfo(bc, lastRun, nBCsPerTF, bcSOR, time, tfIdThis, bcInTF);
            trackRunNumber.push_back(bc.runNumber());
            trackGlobalBC.push_back(bc.globalBC());
            trackTriggerMask.push_back(bc.triggerMask());
            trackTimestamp.push_back(time);
            trackTFidThis.push_back(tfIdThis);
            trackBcInTF.push_back(bcInTF);
          }
          if (bcs.size() == 1) {
            firstGlobalBC = trackGlobalBC[0];
            lastGlobalBC = -10000;
          } else {
            firstGlobalBC = trackGlobalBC[0];
            lastGlobalBC = trackGlobalBC[bcs.size() - 1];
          }
        }
      } else {
        coll = track.template collision_as<C>();
        trkHasCollision = true;
        trkCollVtxZ = coll.posZ();

        if (lastCollId != coll.globalIndex()) {
          trackRunNumber.clear();
          trackGlobalBC.clear();
          trackTriggerMask.clear();
          trackTimestamp.clear();
          trackTFidThis.clear();
          trackBcInTF.clear();
          lastCollId = coll.globalIndex();
          bc = coll.template bc_as<B>();
          getTimingInfo(bc, lastRun, nBCsPerTF, bcSOR, time, tfIdThis, bcInTF);
          trackRunNumber.push_back(bc.runNumber());
          trackGlobalBC.push_back(bc.globalBC());
          trackTriggerMask.push_back(bc.triggerMask());
          trackTimestamp.push_back(time);
          trackTFidThis.push_back(tfIdThis);
          trackBcInTF.push_back(bcInTF);
          firstGlobalBC = trackGlobalBC[0];
          lastGlobalBC = -10000;
        }
      }

      fullTimeFrameIdList.push_back(trackTFidThis[0]);
      if (track.alpha() < 0) {
        lTrkCounter++;
        lowerTrkGIListWithGlobalBCInfo.push_back({lTrkCounter, track.globalIndex(), firstGlobalBC, lastGlobalBC, trackTFidThis[0], trackBcInTF[0]});
        lowerTrkFloatInfo.push_back({trackDcaCalc, trkCollVtxZ});
        lowerTrkBoolInfo.push_back(trkHasCollision);
      }
      if (track.alpha() > 0) {
        uTrkCounter++;
        upperTrkGIListWithGlobalBCInfo.push_back({uTrkCounter, track.globalIndex(), firstGlobalBC, lastGlobalBC, trackTFidThis[0], trackBcInTF[0]});
        upperTrkFloatInfo.push_back({trackDcaCalc, trkCollVtxZ});
        upperTrkBoolInfo.push_back(trkHasCollision);
      }
    }

    std::sort(fullTimeFrameIdList.begin(), fullTimeFrameIdList.end());
    auto last = std::unique(fullTimeFrameIdList.begin(), fullTimeFrameIdList.end());
    fullTimeFrameIdList.erase(last, fullTimeFrameIdList.end());

    histReg.fill(HIST("hTiming_nDFsProcessed"), 0.5);
    histReg.fill(HIST("hTiming_nTFsProcessed"), 0.5, fullTimeFrameIdList.size());
    histReg.fill(HIST("hTiming_nTFsPerDataFrame"), fullTimeFrameIdList.size());

    if (cfgDebug.printDebugMessages) {
      LOG(info) << "DEBUG :: upperTrkGIListWithGlobalBCInfo = " << upperTrkGIListWithGlobalBCInfo.size();
      LOG(info) << "DEBUG :: lowerTrkGIListWithGlobalBCInfo = " << lowerTrkGIListWithGlobalBCInfo.size();
      printTime(trackLoopStart, Form("DEBUG :: df_%d :: DF Reading :: Track loop Time :: ", dfCount));
    }

    sortVectorOfArray(lowerTrkGIListWithGlobalBCInfo, kTrkFirstGBC);
    sortVectorOfArray(upperTrkGIListWithGlobalBCInfo, kTrkFirstGBC);

    auto cosmicMuonLoopStart = std::chrono::high_resolution_clock::now();
    std::vector<int64_t> cosmicMuUTrkBC0;
    std::vector<int64_t> cosmicMuUTrkBC1;
    std::vector<int64_t> cosmicMuLTrkBC0;
    std::vector<int64_t> cosmicMuLTrkBC1;
    std::vector<bool> cosmicMuIsSelectedFlag;

    int cosmicPairCounter = 0;
    float fSumPtPair = 0, fSumQPtPair = 0, uTrkDcaCalc = 0, lTrkDcaCalc = 0;

    constexpr float sP0 = 1.f;         // o2-linter: disable=name/constexpr-constant (algorithm constant)
    constexpr float sP3 = 0.006f;      // o2-linter: disable=name/constexpr-constant (algorithm constant)
    constexpr float sP4 = 0.03f;       // o2-linter: disable=name/constexpr-constant (algorithm constant)
    constexpr float sAlpha = 0.01f;    // o2-linter: disable=name/constexpr-constant (algorithm constant)
    constexpr float sixSquared = 36.f; // o2-linter: disable=name/constexpr-constant (algorithm constant) 6*6 //to avoid sqrt compare with square

    float dP0, dP3, dP4, dAlpha, ptMin;
    float resP0, resP3, resP4, resAlpha, sumSq, dP0corr;
    bool ptCut15, isOKSum6, isOK6, tOK0, isSelected, isOK3Corr, isOK3Clean;

    int uTrkOriginalIdx = 0;
    int lTrkOriginalIdx = 0;

    bool uTrkHasCollison = false;
    bool lTrkHasCollison = false;

    float uTrkCollPosZ = 0;
    float lTrkCollPosZ = 0;

    float tpcTime0[2]; // tpc only time0 (mTime0 in TPC track)
    float tpcdEdxNorm[2];
    int16_t tpcdcaR[2];            // tpc only DCAr
    int16_t tpcdcaZ[2];            // tpc only DCAz
    uint8_t tpcClusterByteMask[2]; // tracklet bitmask - track defining 8 tracklets (152=8*19 rows) bit set if nCluster>thr (default 5)
    uint8_t tpcdEdxMax0R[2];       // TPC dEdxQMax -ROC0/dEdx
    uint8_t tpcdEdxMax1R[2];       // TPC dEdxQMax -ROC1/dEdx
    uint8_t tpcdEdxMax2R[2];       // TPC dEdxQMax -ROC2/dEdx
    uint8_t tpcdEdxMax3R[2];       // TPC dEdxQMax -ROC3/dEdx
    uint8_t tpcdEdxTot0R[2];       // TPC dEdxQtot -ROC0/dEdx
    uint8_t tpcdEdxTot1R[2];       // TPC dEdxQtot -ROC1/dEdx
    uint8_t tpcdEdxTot2R[2];       // TPC dEdxQtot -ROC2/dEdx
    uint8_t tpcdEdxTot3R[2];       // TPC dEdxQtot -ROC3/dEdx

    int8_t deltaRefContParamY[2];
    int8_t deltaRefContParamZ[2];
    int8_t deltaRefContParamSnp[2];
    int8_t deltaRefContParamTgl[2];
    int8_t deltaRefContParamQ2Pt[2];
    int8_t deltaRefGloParamY[2];
    int8_t deltaRefGloParamZ[2];
    int8_t deltaRefGloParamSnp[2];
    int8_t deltaRefGloParamTgl[2];
    int8_t deltaRefGloParamQ2Pt[2];
    int8_t deltaTOFdX[2];
    int8_t deltaTOFdZ[2];

    auto fillTrackQA = [&](auto& trk, int kTrk) {
      if (trk.trackQAId() >= 0) {
        const auto& trkQA = trk.template trackQA_as<Q>();
        if (trkQA.trackId() != trk.globalIndex()) {
          LOG(fatal) << "Mismatch!";
        }
        tpcTime0[kTrk] = trkQA.tpcTime0();
        tpcdEdxNorm[kTrk] = trkQA.tpcdEdxNorm();
        tpcdcaR[kTrk] = trkQA.tpcdcaR();
        tpcdcaZ[kTrk] = trkQA.tpcdcaZ();
        tpcClusterByteMask[kTrk] = trkQA.tpcClusterByteMask();
        tpcdEdxMax0R[kTrk] = trkQA.tpcdEdxMax0R();
        tpcdEdxMax1R[kTrk] = trkQA.tpcdEdxMax1R();
        tpcdEdxMax2R[kTrk] = trkQA.tpcdEdxMax2R();
        tpcdEdxMax3R[kTrk] = trkQA.tpcdEdxMax3R();
        tpcdEdxTot0R[kTrk] = trkQA.tpcdEdxTot0R();
        tpcdEdxTot1R[kTrk] = trkQA.tpcdEdxTot1R();
        tpcdEdxTot2R[kTrk] = trkQA.tpcdEdxTot2R();
        tpcdEdxTot3R[kTrk] = trkQA.tpcdEdxTot3R();
        deltaRefContParamY[kTrk] = trkQA.deltaRefContParamY();
        deltaRefContParamZ[kTrk] = trkQA.deltaRefITSParamZ(); // ← getter is deltaRefITSParamZ !
        // deltaRefContParamZ   [kTrk] = trkQA.deltaRefContParamZ();
        deltaRefContParamSnp[kTrk] = trkQA.deltaRefContParamSnp();
        deltaRefContParamTgl[kTrk] = trkQA.deltaRefContParamTgl();
        deltaRefContParamQ2Pt[kTrk] = trkQA.deltaRefContParamQ2Pt();
        deltaRefGloParamY[kTrk] = trkQA.deltaRefGloParamY();
        deltaRefGloParamZ[kTrk] = trkQA.deltaRefGloParamZ();
        deltaRefGloParamSnp[kTrk] = trkQA.deltaRefGloParamSnp();
        deltaRefGloParamTgl[kTrk] = trkQA.deltaRefGloParamTgl();
        deltaRefGloParamQ2Pt[kTrk] = trkQA.deltaRefGloParamQ2Pt();
        deltaTOFdX[kTrk] = trkQA.deltaTOFdX();
        deltaTOFdZ[kTrk] = trkQA.deltaTOFdZ();
      } else {
        tpcTime0[kTrk] = -999.f;
        tpcdEdxNorm[kTrk] = -999.f;
        tpcdcaR[kTrk] = (int16_t)-999;
        tpcdcaZ[kTrk] = (int16_t)-999;
        tpcClusterByteMask[kTrk] = 0x00;
        tpcdEdxMax0R[kTrk] = 0x00;
        tpcdEdxMax1R[kTrk] = 0x00;
        tpcdEdxMax2R[kTrk] = 0x00;
        tpcdEdxMax3R[kTrk] = 0x00;
        tpcdEdxTot0R[kTrk] = 0x00;
        tpcdEdxTot1R[kTrk] = 0x00;
        tpcdEdxTot2R[kTrk] = 0x00;
        tpcdEdxTot3R[kTrk] = 0x00;

        deltaRefContParamY[kTrk] = (int8_t)-128;
        deltaRefContParamZ[kTrk] = (int8_t)-128;
        deltaRefContParamSnp[kTrk] = (int8_t)-128;
        deltaRefContParamTgl[kTrk] = (int8_t)-128;
        deltaRefContParamQ2Pt[kTrk] = (int8_t)-128;
        deltaRefGloParamY[kTrk] = (int8_t)-128;
        deltaRefGloParamZ[kTrk] = (int8_t)-128;
        deltaRefGloParamSnp[kTrk] = (int8_t)-128;
        deltaRefGloParamTgl[kTrk] = (int8_t)-128;
        deltaRefGloParamQ2Pt[kTrk] = (int8_t)-128;
        deltaTOFdX[kTrk] = (int8_t)-128;
        deltaTOFdZ[kTrk] = (int8_t)-128;
      }
    };

    float tmoPrimUnfm80[2];
    float tmoRobustT0V0PrimUnfm80[2];
    float twmoPrimUnfm80[2];
    float twmoRobustT0V0PrimUnfm80[2];

    auto fillOccupancies = [&](auto& trk, int kTrk) {
      if (trk.tmoId() >= 0) {
        const auto& trkTMO = trk.template tmo_as<O>();
        if (trkTMO.trackId() != trk.globalIndex()) {
          LOG(fatal) << "Mismatch!";
        }
        tmoPrimUnfm80[kTrk] = trkTMO.tmoPrimUnfm80();
        tmoRobustT0V0PrimUnfm80[kTrk] = trkTMO.tmoRobustT0V0PrimUnfm80();
        twmoPrimUnfm80[kTrk] = trkTMO.twmoPrimUnfm80();
        twmoRobustT0V0PrimUnfm80[kTrk] = trkTMO.twmoRobustT0V0PrimUnfm80();
      } else {
        tmoPrimUnfm80[kTrk] = -999.f;
        tmoRobustT0V0PrimUnfm80[kTrk] = -999.f;
        twmoPrimUnfm80[kTrk] = -999.f;
        twmoRobustT0V0PrimUnfm80[kTrk] = -999.f;
      }
    };

    for (uint idxUTrk = 0; idxUTrk < upperTrkGIListWithGlobalBCInfo.size(); idxUTrk++) {
      const auto& uTrk = upperTrkGIListWithGlobalBCInfo[idxUTrk];
      const auto& upperTrk = tracks.iteratorAt(uTrk[kTrkGI]);
      uTrkOriginalIdx = upperTrkGIListWithGlobalBCInfo[idxUTrk][kTrkCountIdx];
      const auto& uBC = uTrk[kTrkFirstGBC];
      // ====================================================================
      // BINARY SEARCH: Find time window [uBC - deltaBC, uBC + deltaBC]
      // ====================================================================
      const int idxLow = binarySearchLower(lowerTrkGIListWithGlobalBCInfo, kTrkFirstGBC, uBC - cfgCM.deltaBC);
      const int idxHigh = binarySearchUpper(lowerTrkGIListWithGlobalBCInfo, kTrkFirstGBC, uBC + cfgCM.deltaBC);

      // Check pointer is valid first
      if constexpr (mode == kProcessWithTracksQA || mode == kProcessFull) {
        fillTrackQA(upperTrk, kUpperTrk);
      }

      if constexpr (mode == kProcessWithOccupancies || mode == kProcessFull) {
        fillOccupancies(upperTrk, kUpperTrk);
      }

      // Skip if no candidates in time window
      if (idxLow >= idxHigh)
        continue;

      for (int idxLTrk = idxLow; idxLTrk < idxHigh; ++idxLTrk) {
        const auto& lTrk = lowerTrkGIListWithGlobalBCInfo[idxLTrk];
        const auto& lowerTrk = tracks.iteratorAt(lTrk[kTrkGI]);

        // if(uTrk[kTrkTfIdThis] != lTrk[kTrkTfIdThis]) {continue;}

        fSumPtPair = upperTrk.pt() + lowerTrk.pt();
        fSumQPtPair = upperTrk.signed1Pt() + lowerTrk.signed1Pt();

        lTrkOriginalIdx = lowerTrkGIListWithGlobalBCInfo[idxLTrk][kTrkCountIdx];

        uTrkDcaCalc = upperTrkFloatInfo[uTrkOriginalIdx][kTrkCalcDcaInfo];
        lTrkDcaCalc = lowerTrkFloatInfo[lTrkOriginalIdx][kTrkCalcDcaInfo];

        uTrkCollPosZ = upperTrkFloatInfo[uTrkOriginalIdx][kTrkCollPosZ];
        lTrkCollPosZ = lowerTrkFloatInfo[lTrkOriginalIdx][kTrkCollPosZ];

        uTrkHasCollison = upperTrkBoolInfo[uTrkOriginalIdx];
        lTrkHasCollison = lowerTrkBoolInfo[lTrkOriginalIdx];

        if (std::abs(fSumPtPair) < cfgCM.csmMuSumPtPair) {
          continue;
        }
        if (std::abs(fSumQPtPair) > cfgCM.csmMuSumQPtPair) {
          continue;
        }

        if (!(std::abs(std::abs(upperTrk.signed1Pt()) - std::abs(lowerTrk.signed1Pt())) < cfgCM.csmMuDiffQPtPair)) {
          continue;
        }
        if (!(std::abs(std::abs(upperTrk.tgl()) - std::abs(lowerTrk.tgl())) < cfgCM.csmMuDiffTglPair)) {
          continue;
        }
        if (!(std::abs(std::abs(uTrkDcaCalc) - std::abs(lTrkDcaCalc)) < cfgCM.csmMuDiffDcaCalcPair)) {
          continue;
        }

        cosmicMuUTrkBC0.push_back(uTrk[kTrkFirstGBC]);
        cosmicMuUTrkBC1.push_back(uTrk[kTrkLastGBC]);
        cosmicMuLTrkBC0.push_back(lTrk[kTrkFirstGBC]);
        cosmicMuLTrkBC1.push_back(lTrk[kTrkLastGBC]);

        dP0 = lowerTrk.dcaXY() + upperTrk.dcaXY();
        dP3 = lowerTrk.tgl() + upperTrk.tgl();
        dP4 = lowerTrk.signed1Pt() + upperTrk.signed1Pt();
        dAlpha = upperTrk.alpha() - lowerTrk.alpha() - o2::constants::math::PI;

        ptMin = (lowerTrk.pt() < upperTrk.pt()) ? lowerTrk.pt() : upperTrk.pt();
        ptCut15 = ptMin > 1.5f; // o2-linter: disable=magic-number (algorithm constant) this is residual distortion dependent

        // Compute square of normalised residuals for fast calculation
        resP0 = dP0 / sP0;
        resP3 = dP3 / sP3;
        resP4 = dP4 / sP4;
        resAlpha = dAlpha / sAlpha;
        sumSq = resP0 * resP0 + resP3 * resP3 + resP4 * resP4 + resAlpha * resAlpha;
        isOKSum6 = sumSq < sixSquared;
        isOK6 = std::abs(resP0) < 6 && std::abs(resP3) < 6 && std::abs(resP4) < 6 && std::abs(dAlpha / sAlpha) < 6; // o2-linter: disable=magic-number (algorithm constant)
        tOK0 = std::abs(uTrk[kTrkBCinTF] - lTrk[kTrkBCinTF]) < 4000;                                                // o2-linter: disable=magic-number (algorithm constant)
        isSelected = isOK6 && tOK0 && ptCut15;

        dP0corr = 1.5 * ((dP4 / 667.0) * 85.0 * 85.0) * (bFieldKG / -5.);
        isOK3Corr = std::abs((dP0 - dP0corr) / sP0) < 3;                    // o2-linter: disable=magic-number (algorithm constant)
        isOK3Clean = isOK6 && tOK0 && isOK3Corr && std::abs(dP3 / sP3) < 3; // o2-linter: disable=magic-number (algorithm constant)

        if (cfgCM.csmMuDoBkgDownSampling && !isSelected && dist(rng) > cfgCM.csmMuBkgDownSampling) {
          continue; // Downsampling the background
        }

        cosmicMuIsSelectedFlag.push_back(isSelected);

        if constexpr (mode == kProcessWithTracksQA || mode == kProcessFull) {
          fillTrackQA(lowerTrk, kLowerTrk);
        }

        if constexpr (mode == kProcessWithOccupancies || mode == kProcessFull) {
          fillOccupancies(lowerTrk, kLowerTrk);
        }

        cosmicPairCounter++;

        if (cfgCM.produceTable00)
          genCosmicPairs00(
            // Cosmic pair information variables
            bFieldKG, ptCut15, isOKSum6, isOK6, tOK0, isSelected, dP0corr, isOK3Corr, isOK3Clean);

        if (cfgCM.produceTable01)
          genCosmicPairs01(
            upperTrk.globalIndex(), uTrk[kTrkFirstGBC], uTrk[kTrkLastGBC], uTrk[kTrkTfIdThis], uTrk[kTrkBCinTF]

            ,
            lowerTrk.globalIndex(), lTrk[kTrkFirstGBC], lTrk[kTrkLastGBC], lTrk[kTrkTfIdThis], lTrk[kTrkBCinTF]);

        if (cfgCM.produceTable02)
          genCosmicPairs02(
            upperTrk.x(), upperTrk.alpha(), upperTrk.y(), upperTrk.z(), upperTrk.snp(), upperTrk.tgl(), upperTrk.signed1Pt(), upperTrk.pt(), upperTrk.p(), upperTrk.eta(), upperTrk.phi(),
            lowerTrk.x(), lowerTrk.alpha(), lowerTrk.y(), lowerTrk.z(), lowerTrk.snp(), lowerTrk.tgl(), lowerTrk.signed1Pt(), lowerTrk.pt(), lowerTrk.p(), lowerTrk.eta(), lowerTrk.phi());

        if (cfgCM.produceTable03)
          genCosmicPairs03(
            uTrkDcaCalc, upperTrk.dcaXY(), upperTrk.dcaZ(), upperTrk.sigmaDcaXY2(), upperTrk.sigmaDcaZ2(), lTrkDcaCalc, lowerTrk.dcaXY(), lowerTrk.dcaZ(), lowerTrk.sigmaDcaXY2(), lowerTrk.sigmaDcaZ2());

        if (cfgCM.produceTable04)
          genCosmicPairs04(
            uTrkHasCollison, uTrkCollPosZ, lTrkHasCollison, lTrkCollPosZ);

        if (cfgCM.produceTable05)
          genCosmicPairs05(
            upperTrk.flags(), upperTrk.length(), upperTrk.pidForTracking(), upperTrk.isPVContributor(), upperTrk.trackTime(), upperTrk.trackTimeRes(), upperTrk.detectorMap(),
            lowerTrk.flags(), lowerTrk.length(), lowerTrk.pidForTracking(), lowerTrk.isPVContributor(), lowerTrk.trackTime(), lowerTrk.trackTimeRes(), lowerTrk.detectorMap());

        if (cfgCM.produceTable06)
          genCosmicPairs06(
            upperTrk.tpcInnerParam(),
            upperTrk.tpcNClsFindable(),
            upperTrk.tpcNClsFindableMinusFound(),
            upperTrk.tpcNClsFindableMinusPID(),
            upperTrk.tpcNClsFindableMinusCrossedRows(),
            upperTrk.tpcNClsShared(),
            upperTrk.tpcDeltaTFwd(),
            upperTrk.tpcDeltaTBwd(),
            upperTrk.tpcChi2NCl(),
            upperTrk.tpcSignal(),
            upperTrk.hasTPC(),
            upperTrk.tpcNClsFound(),
            upperTrk.tpcNClsPID(),
            upperTrk.tpcNClsCrossedRows(),
            upperTrk.tpcCrossedRowsOverFindableCls(),
            upperTrk.tpcFoundOverFindableCls(),
            upperTrk.tpcFractionSharedCls(),

            lowerTrk.tpcInnerParam(),
            lowerTrk.tpcNClsFindable(),
            lowerTrk.tpcNClsFindableMinusFound(),
            lowerTrk.tpcNClsFindableMinusPID(),
            lowerTrk.tpcNClsFindableMinusCrossedRows(),
            lowerTrk.tpcNClsShared(),
            lowerTrk.tpcDeltaTFwd(),
            lowerTrk.tpcDeltaTBwd(),
            lowerTrk.tpcChi2NCl(),
            lowerTrk.tpcSignal(),
            lowerTrk.hasTPC(),
            lowerTrk.tpcNClsFound(),
            lowerTrk.tpcNClsPID(),
            lowerTrk.tpcNClsCrossedRows(),
            lowerTrk.tpcCrossedRowsOverFindableCls(),
            lowerTrk.tpcFoundOverFindableCls(),
            lowerTrk.tpcFractionSharedCls());

        if (cfgCM.produceTable07)
          genCosmicPairs07(
            upperTrk.itsClusterSizes(),
            upperTrk.itsChi2NCl(),
            upperTrk.hasITS(),
            upperTrk.itsClusterMap(),
            upperTrk.itsNCls(),
            upperTrk.itsNClsInnerBarrel(),
            // upperTrk.itsClsSizeInLayer(),
            upperTrk.isITSAfterburner(),

            lowerTrk.itsClusterSizes(),
            lowerTrk.itsChi2NCl(),
            lowerTrk.hasITS(),
            lowerTrk.itsClusterMap(),
            lowerTrk.itsNCls(),
            lowerTrk.itsNClsInnerBarrel(),
            // lowerTrk.itsClsSizeInLayer(),
            lowerTrk.isITSAfterburner());

        if (cfgCM.produceTable08)
          genCosmicPairs08(
            upperTrk.trdPattern(),
            upperTrk.trdChi2(),
            upperTrk.trdSignal(),
            upperTrk.hasTRD(),

            lowerTrk.trdPattern(),
            lowerTrk.trdChi2(),
            lowerTrk.trdSignal(),
            lowerTrk.hasTRD());

        if (cfgCM.produceTable09)
          genCosmicPairs09(
            upperTrk.tofChi2(),
            upperTrk.tofExpMom(),
            upperTrk.hasTOF(),
            upperTrk.tofExpTimeEl(),
            upperTrk.tofExpTimeMu(),
            upperTrk.tofExpTimePi(),
            upperTrk.tofExpTimeKa(),
            upperTrk.tofExpTimePr(),
            upperTrk.tofExpTimeDe(),
            upperTrk.tofExpTimeTr(),
            upperTrk.tofExpTimeHe(),
            upperTrk.tofExpTimeAl(),

            lowerTrk.tofChi2(),
            lowerTrk.tofExpMom(),
            lowerTrk.hasTOF(),
            lowerTrk.tofExpTimeEl(),
            lowerTrk.tofExpTimeMu(),
            lowerTrk.tofExpTimePi(),
            lowerTrk.tofExpTimeKa(),
            lowerTrk.tofExpTimePr(),
            lowerTrk.tofExpTimeDe(),
            lowerTrk.tofExpTimeTr(),
            lowerTrk.tofExpTimeHe(),
            lowerTrk.tofExpTimeAl());

        if (cfgCM.produceTable10)
          genCosmicPairs10(
            upperTrk.trackEtaEmcal(),
            upperTrk.trackPhiEmcal(),

            lowerTrk.trackEtaEmcal(),
            lowerTrk.trackPhiEmcal());

        if constexpr (mode == kProcessWithTOFSignal || mode == kProcessFull) {
          if (cfgCM.produceTable11)
            genCosmicPairs11(
              upperTrk.goodTOFMatch(),
              upperTrk.tofSignal(),
              // upperTrk.eventCollisionTime(),

              lowerTrk.goodTOFMatch(),
              lowerTrk.tofSignal()
              // lowerTrk.eventCollisionTime()
            );
        }
        if constexpr (mode == kProcessWithPidEvTimeFlags) {
          if (cfgCM.produceTable12)
            genCosmicPairs12(
              upperTrk.tofFlags(),
              upperTrk.isEvTimeDefined(),
              upperTrk.isEvTimeTOF(),
              upperTrk.isEvTimeT0AC(),
              upperTrk.isEvTimeTOFT0AC(),

              lowerTrk.tofFlags(),
              lowerTrk.isEvTimeDefined(),
              lowerTrk.isEvTimeTOF(),
              lowerTrk.isEvTimeT0AC(),
              lowerTrk.isEvTimeTOFT0AC());
        }

        if constexpr (mode == kProcessWithTOFEvTime || mode == kProcessFull) {
          if (cfgCM.produceTable13)
            genCosmicPairs13(
              upperTrk.tofEvTime(),
              upperTrk.tofEvTimeErr(),
              lowerTrk.tofEvTime(),
              lowerTrk.tofEvTimeErr());
        }

        if constexpr (mode == kProcessWithTracksQA || mode == kProcessFull) {
          if (cfgCM.produceTable14)
            genCosmicPairs14(
              tpcTime0[kUpperTrk],
              tpcdEdxNorm[kUpperTrk],
              tpcdcaR[kUpperTrk],
              tpcdcaZ[kUpperTrk],
              tpcClusterByteMask[kUpperTrk],
              tpcdEdxMax0R[kUpperTrk],
              tpcdEdxMax1R[kUpperTrk],
              tpcdEdxMax2R[kUpperTrk],
              tpcdEdxMax3R[kUpperTrk],
              tpcdEdxTot0R[kUpperTrk],
              tpcdEdxTot1R[kUpperTrk],
              tpcdEdxTot2R[kUpperTrk],
              tpcdEdxTot3R[kUpperTrk],

              tpcTime0[kLowerTrk],
              tpcdEdxNorm[kLowerTrk],
              tpcdcaR[kLowerTrk],
              tpcdcaZ[kLowerTrk],
              tpcClusterByteMask[kLowerTrk],
              tpcdEdxMax0R[kLowerTrk],
              tpcdEdxMax1R[kLowerTrk],
              tpcdEdxMax2R[kLowerTrk],
              tpcdEdxMax3R[kLowerTrk],
              tpcdEdxTot0R[kLowerTrk],
              tpcdEdxTot1R[kLowerTrk],
              tpcdEdxTot2R[kLowerTrk],
              tpcdEdxTot3R[kLowerTrk]);

          if (cfgCM.produceTable15)
            genCosmicPairs15(
              deltaRefContParamY[kUpperTrk],
              deltaRefContParamZ[kUpperTrk],
              deltaRefContParamSnp[kUpperTrk],
              deltaRefContParamTgl[kUpperTrk],
              deltaRefContParamQ2Pt[kUpperTrk],
              deltaRefGloParamY[kUpperTrk],
              deltaRefGloParamZ[kUpperTrk],
              deltaRefGloParamSnp[kUpperTrk],
              deltaRefGloParamTgl[kUpperTrk],
              deltaRefGloParamQ2Pt[kUpperTrk],
              deltaTOFdX[kUpperTrk],
              deltaTOFdZ[kUpperTrk],

              deltaRefContParamY[kLowerTrk],
              deltaRefContParamZ[kLowerTrk],
              deltaRefContParamSnp[kLowerTrk],
              deltaRefContParamTgl[kLowerTrk],
              deltaRefContParamQ2Pt[kLowerTrk],
              deltaRefGloParamY[kLowerTrk],
              deltaRefGloParamZ[kLowerTrk],
              deltaRefGloParamSnp[kLowerTrk],
              deltaRefGloParamTgl[kLowerTrk],
              deltaRefGloParamQ2Pt[kLowerTrk],
              deltaTOFdX[kLowerTrk],
              deltaTOFdZ[kLowerTrk]);
        }

        if constexpr (mode == kProcessWithOccupancies || mode == kProcessFull) {
          if (cfgCM.produceTable16)
            genCosmicPairs16(
              tmoPrimUnfm80[kUpperTrk],
              tmoRobustT0V0PrimUnfm80[kUpperTrk],
              twmoPrimUnfm80[kUpperTrk],
              twmoRobustT0V0PrimUnfm80[kUpperTrk],

              tmoPrimUnfm80[kLowerTrk],
              tmoRobustT0V0PrimUnfm80[kLowerTrk],
              twmoPrimUnfm80[kLowerTrk],
              twmoRobustT0V0PrimUnfm80[kLowerTrk]);
        }
      }
    }

    int isSelCosmicPairCounter = std::count(cosmicMuIsSelectedFlag.begin(), cosmicMuIsSelectedFlag.end(), true);

    histReg.fill(HIST("hCosmicPairs_nCosmicsPerDF"), cosmicPairCounter);
    histReg.fill(HIST("hCosmicPairs_nIsSelCosmicsPerDF"), isSelCosmicPairCounter);

    if (cfgDebug.printDebugMessages) {
      printTime(cosmicMuonLoopStart, Form("DEBUG :: df_%d :: DF Reading :: cosmicMuonLoop Time :: ", dfCount));
      LOG(info) << "DEBUG :: df_" << dfCount << " :: cosmic counted = " << cosmicPairCounter;
    }

    // Collision looping and flag table creator
    auto collLoopStart = std::chrono::high_resolution_clock::now();
    std::unordered_set<int> finalCosmicIndexSet; // to keep only unique entries.
    std::vector<int> tempCosmicIndexList;
    std::vector<int> finalCosmicIndexList;
    bool flagHasCosmicMuon = false;
    bool flagHasIsSelectedCosmic = false;

    int counterCollVtx10 = 0;
    int counterCollVtx8 = 0;
    int counterFlagHasCosmic = 0;
    int counterFlagHasIsSelectedCosmic = 0;
    int counterVtx10FlagHasCosmic = 0;
    int counterVtx10FlagHasIsSelectedCosmic = 0;
    int counterVtx8FlagHasCosmic = 0;
    int counterVtx8FlagHasIsSelectedCosmic = 0;

    auto addIndices = [](const std::vector<int>& indices, std::unordered_set<int>& indexSet, bool& flagHasCosmic) {
      if (!indices.empty()) {
        flagHasCosmic = true;
        for (const auto& idx : indices)
          indexSet.insert(idx);
      }
    };

    for (const auto& coll : colls) {
      const auto& bc = coll.template bc_as<B>();
      getTimingInfo(bc, lastRun, nBCsPerTF, bcSOR, time, tfIdThis, bcInTF);

      if (nBCsPerTF > nBCinTF) {
        LOG(fatal) << "DEBUG :: FATAL ERROR :: nBCsPerTF > nBCinTF i.e " << nBCsPerTF << " > " << nBCinTF << " will cause crash in further process";
        return;
      }

      flagHasCosmicMuon = false;
      finalCosmicIndexSet.clear();

      addIndices(getCosmicTriggerList(bc.globalBC(), cosmicMuUTrkBC0), finalCosmicIndexSet, flagHasCosmicMuon);
      addIndices(getCosmicTriggerList(bc.globalBC(), cosmicMuUTrkBC1), finalCosmicIndexSet, flagHasCosmicMuon);
      addIndices(getCosmicTriggerList(bc.globalBC(), cosmicMuLTrkBC0), finalCosmicIndexSet, flagHasCosmicMuon);
      addIndices(getCosmicTriggerList(bc.globalBC(), cosmicMuLTrkBC1), finalCosmicIndexSet, flagHasCosmicMuon);

      finalCosmicIndexList.clear();
      finalCosmicIndexList.assign(finalCosmicIndexSet.begin(), finalCosmicIndexSet.end());

      flagHasIsSelectedCosmic = false;
      for (const auto& idx : finalCosmicIndexSet) {
        if (cosmicMuIsSelectedFlag[idx]) {
          flagHasIsSelectedCosmic = true;
          break;
        }
      }

      const float absZ = std::abs(coll.posZ());
      const bool z8 = absZ < 8.0f;   // o2-linter: disable=magic-number (algorithm constant)
      const bool z10 = absZ < 10.0f; // o2-linter: disable=magic-number (algorithm constant)

      counterCollVtx8 += z8;
      counterCollVtx10 += z10;

      if (flagHasCosmicMuon) {
        counterFlagHasCosmic++;
        counterVtx8FlagHasCosmic += z8;
        counterVtx10FlagHasCosmic += z10;
      }

      if (flagHasIsSelectedCosmic) {
        counterFlagHasIsSelectedCosmic++;
        counterVtx8FlagHasIsSelectedCosmic += z8;
        counterVtx10FlagHasIsSelectedCosmic += z10;
      }

      genCollisionCosmicFlags(flagHasCosmicMuon, flagHasIsSelectedCosmic);
      genCollToCosmic(finalCosmicIndexList);
    }

    auto safeRatio = [](double num, double den) { return den > 0 ? num / den : 0.; };

    const double nAll = colls.size(), nV8 = counterCollVtx8, nV10 = counterCollVtx10;
    const double nCos = counterFlagHasCosmic, nSel = counterFlagHasIsSelectedCosmic;
    const double nCosV8 = counterVtx8FlagHasCosmic, nCosV10 = counterVtx10FlagHasCosmic;
    const double nSelV8 = counterVtx8FlagHasIsSelectedCosmic, nSelV10 = counterVtx10FlagHasIsSelectedCosmic;

    histReg.fill(HIST("hRatio"), safeRatio(nCos, nAll));
    histReg.fill(HIST("hRatioVtx8"), safeRatio(nCosV8, nV8));
    histReg.fill(HIST("hRatioVtx10"), safeRatio(nCosV10, nV10));

    histReg.fill(HIST("hRatioIsSel"), safeRatio(nSel, nAll));
    histReg.fill(HIST("hRatioIsSelVtx8"), safeRatio(nSelV8, nV8));
    histReg.fill(HIST("hRatioIsSelVtx10"), safeRatio(nSelV10, nV10));

    histReg.fill(HIST("hCosmicFlag"), nCos);
    histReg.fill(HIST("hCosmicFlagVtx8"), nCosV8);
    histReg.fill(HIST("hCosmicFlagVtx10"), nCosV10);
    histReg.fill(HIST("hCosmicFlagIsSel"), nSel);
    histReg.fill(HIST("hCosmicFlagIsSelVtx8"), nSelV8);
    histReg.fill(HIST("hCosmicFlagIsSelVtx10"), nSelV10);
    histReg.fill(HIST("hColls"), nAll);
    histReg.fill(HIST("hCollsVtx8"), nV8);
    histReg.fill(HIST("hCollsVtx10"), nV10);

    if (cfgDebug.printDebugMessages) {
      printTime(collLoopStart, Form("DEBUG :: df_%d :: DF Reading :: coll loop Time :: ", dfCount));
      printTime(start1, Form("DEBUG :: df_%d :: DF End    :: DF Read Time :: ", dfCount));
      printTime(start0, Form("DEBUG :: df_%d :: DF End    :: Elapsed Time :: ", dfCount));
      LOG(info) << "DEBUG ::";
    }

  } // Process function ends

  using MyTracks = soa::Join<aod::Tracks, aod::TrackToAmbgTrk, aod::TracksExtra, aod::TracksDCA, aod::TracksDCACov>;
  using MyTracks1 = soa::Join<aod::Tracks, aod::TrackToAmbgTrk, aod::TracksExtra, aod::TracksDCA, aod::TracksDCACov, aod::pidTOFFlags, aod::TOFSignal>;
  using MyTracks2 = soa::Join<aod::Tracks, aod::TrackToAmbgTrk, aod::TracksExtra, aod::TracksDCA, aod::TracksDCACov, aod::pidTOFFlags, aod::TOFSignal, aod::TOFEvTime>;
  using MyTracks3 = soa::Join<aod::Tracks, aod::TrackToAmbgTrk, aod::TracksExtra, aod::TracksDCA, aod::TracksDCACov, aod::pidTOFFlags, aod::TOFSignal, aod::TOFEvTime, aod::pidEvTimeFlags>;
  using MyTracks4 = soa::Join<aod::Tracks, aod::TrackToAmbgTrk, aod::TracksExtra, aod::TracksDCA, aod::TracksDCACov, aod::pidTOFFlags, aod::TOFSignal, aod::TrackToTracksQA>;
  using MyTracks5 = soa::Join<aod::Tracks, aod::TrackToAmbgTrk, aod::TracksExtra, aod::TracksDCA, aod::TracksDCACov, aod::pidTOFFlags, aod::TOFSignal, aod::TrackToTracksQA, aod::TrackToTmo>;
  using MyTracks6 = soa::Join<aod::Tracks, aod::TrackToAmbgTrk, aod::TracksExtra, aod::TracksDCA, aod::TracksDCACov, aod::pidTOFFlags, aod::TOFSignal, aod::TOFEvTime, aod::TrackToTracksQA, aod::TrackToTmo>;

  using MyTrackMeanOccs = soa::Join<aod::TmoTrackIds, aod::TmoPrim, aod::TmoRT0V0Prim, aod::TwmoPrim, aod::TwmoRT0V0Prim>;

  // Mode 0: bare tracks, no TOF, no QA, no occupancy
  void processBasic(aod::BCsWithTimestamps const& BCs,
                    aod::Collisions const& collisions, MyTracks const& tracks, o2::aod::AmbiguousTracks const& ambgTracks, o2::aod::Origins const& Origins)
  {
    executeProcess<kProcessBasic>(BCs, collisions, tracks, ambgTracks, nullptr, nullptr, Origins);
  }
  PROCESS_SWITCH(CosmicMuonToCollisionAssociator, processBasic, "processBasic", false);

  // Mode 1: with TOF signal + good match flag
  void processWithTOFSignal(aod::BCsWithTimestamps const& BCs,
                            aod::Collisions const& collisions, MyTracks1 const& tracks, o2::aod::AmbiguousTracks const& ambgTracks, o2::aod::Origins const& Origins)
  {
    executeProcess<kProcessWithTOFSignal>(BCs, collisions, tracks, ambgTracks, nullptr, nullptr, Origins);
  }
  PROCESS_SWITCH(CosmicMuonToCollisionAssociator, processWithTOFSignal, "processWithTOFSignal", false);

  // Mode 2: with TOF event time (tofEvTime + tofEvTimeErr)
  void processWithTOFEvTime(aod::BCsWithTimestamps const& BCs,
                            aod::Collisions const& collisions, MyTracks2 const& tracks, o2::aod::AmbiguousTracks const& ambgTracks, o2::aod::Origins const& Origins)
  {
    executeProcess<kProcessWithTOFEvTime>(BCs, collisions, tracks, ambgTracks, nullptr, nullptr, Origins);
  }
  PROCESS_SWITCH(CosmicMuonToCollisionAssociator, processWithTOFEvTime, "processWithTOFEvTime", false);

  // Mode 3: with pidEvTimeFlags (NOTE: currently failing — table size mismatch) need to be fixed
  void processWithPidEvTimeFlags(aod::BCsWithTimestamps const& BCs,
                                 aod::Collisions const& collisions, MyTracks3 const& tracks, o2::aod::AmbiguousTracks const& ambgTracks, o2::aod::Origins const& Origins)
  {
    executeProcess<kProcessWithPidEvTimeFlags>(BCs, collisions, tracks, ambgTracks, nullptr, nullptr, Origins);
  }
  PROCESS_SWITCH(CosmicMuonToCollisionAssociator, processWithPidEvTimeFlags, "processWithPidEvTimeFlags", false);

  // Mode 4: with TrackQA (TPC time0, dEdx, DCA, deltaRef params)
  void processWithTrackQA(aod::BCsWithTimestamps const& BCs,
                          aod::Collisions const& collisions, MyTracks4 const& tracks, o2::aod::AmbiguousTracks const& ambgTracks, o2::aod::TracksQAVersion const& tracksQA, o2::aod::Origins const& Origins)
  {
    executeProcess<kProcessWithTracksQA>(BCs, collisions, tracks, ambgTracks, tracksQA, nullptr, Origins);
  }
  PROCESS_SWITCH(CosmicMuonToCollisionAssociator, processWithTrackQA, "processWithTrackQA", false);

  void processWithOccupancies(aod::BCsWithTimestamps const& BCs,
                              aod::Collisions const& collisions, MyTracks5 const& tracks, o2::aod::AmbiguousTracks const& ambgTracks, MyTrackMeanOccs const& trackMeanOccs, o2::aod::Origins const& Origins)
  {
    executeProcess<kProcessWithOccupancies>(BCs, collisions, tracks, ambgTracks, nullptr, trackMeanOccs, Origins);
  }
  PROCESS_SWITCH(CosmicMuonToCollisionAssociator, processWithOccupancies, "processWithOccupancies", false);

  void processFullMode(aod::BCsWithTimestamps const& BCs,
                       aod::Collisions const& collisions, MyTracks6 const& tracks, o2::aod::AmbiguousTracks const& ambgTracks, o2::aod::TracksQAVersion const& tracksQA, MyTrackMeanOccs const& trackMeanOccs, o2::aod::Origins const& Origins)
  {
    executeProcess<kProcessFull>(BCs, collisions, tracks, ambgTracks, tracksQA, trackMeanOccs, Origins);
  }
  PROCESS_SWITCH(CosmicMuonToCollisionAssociator, processFullMode, "processFullMode", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CosmicMuonToCollisionAssociator>(cfgc)};
}
