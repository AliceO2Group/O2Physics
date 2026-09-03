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
/// \file   TauThreeProngEventTables.h
/// \author Adam Matyja <adam.tomasz.matyja@cern.ch>, IFJ PAN, Cracow
/// \since  2025-09-06
/// \brief  A table to store information about events preselected to be candidates for UPC gammagamma->tautau in 1+3 ot 3+3 topology
///

#ifndef PWGUD_DATAMODEL_TAUTHREEPRONGEVENTTABLES_H_
#define PWGUD_DATAMODEL_TAUTHREEPRONGEVENTTABLES_H_

#include <Framework/AnalysisDataModel.h>

#include <cstdint>
// derived tables for tautau->4 (=1+3) tracks
// derived tables for tautau->6 (=3+3) tracks
namespace o2::aod
{
namespace tautree
{
// event info
DECLARE_SOA_COLUMN(RunNumber, runNumber, int32_t);
DECLARE_SOA_COLUMN(Bc, bc, int);
DECLARE_SOA_COLUMN(TotalTracks, totalTracks, int);
DECLARE_SOA_COLUMN(NumContrib, numContrib, int8_t);
DECLARE_SOA_COLUMN(RctOk, rctOk, int);
// DECLARE_SOA_COLUMN(GlobalNonPVtracks, globalNonPVtracks, int);
// DECLARE_SOA_COLUMN(PosX, posX, float);
// DECLARE_SOA_COLUMN(PosY, posY, float);
DECLARE_SOA_COLUMN(PosZ, posZ, float);
DECLARE_SOA_COLUMN(FlagUPC, flagUPC, int8_t);
DECLARE_SOA_COLUMN(OccupancyInTime, occupancyInTime, int);
DECLARE_SOA_COLUMN(HadronicRate, hadronicRate, double);

// DECLARE_SOA_COLUMN(Trs, trs, int8_t);
// DECLARE_SOA_COLUMN(Trofs, trofs, int8_t);
// DECLARE_SOA_COLUMN(Hmpr, hmpr, int8_t);
// DECLARE_SOA_COLUMN(Tfb, tfb, int8_t);
// DECLARE_SOA_COLUMN(ItsRofb, itsRofb, int8_t);
// DECLARE_SOA_COLUMN(Sbp, sbp, int8_t);
// DECLARE_SOA_COLUMN(ZvtxFT0vsPv, zvtxFT0vsPv, int8_t);
// DECLARE_SOA_COLUMN(VtxITSTPC, vtxITSTPC, int8_t);
DECLARE_SOA_COLUMN(BcSelBits, bcSelBits, uint8_t);

DECLARE_SOA_COLUMN(ZdcAenergy, zdcAenergy, float);
DECLARE_SOA_COLUMN(ZdcCenergy, zdcCenergy, float);
DECLARE_SOA_COLUMN(ZdcAtime, zdcAtime, float);
DECLARE_SOA_COLUMN(ZdcCtime, zdcCtime, float);
// DECLARE_SOA_COLUMN(Qtot, qtot, int8_t);
// FIT info
DECLARE_SOA_COLUMN(TotalFT0AmplitudeA, totalFT0AmplitudeA, float);
DECLARE_SOA_COLUMN(TotalFT0AmplitudeC, totalFT0AmplitudeC, float);
DECLARE_SOA_COLUMN(TotalFV0AmplitudeA, totalFV0AmplitudeA, float);
// DECLARE_SOA_COLUMN(TimeFT0A, timeFT0A, float);
// DECLARE_SOA_COLUMN(TimeFT0C, timeFT0C, float);
// DECLARE_SOA_COLUMN(TimeFV0A, timeFV0A, float);
// tracks - 4 tracks
DECLARE_SOA_COLUMN(TrkPx, trkPx, float[4]);
DECLARE_SOA_COLUMN(TrkPy, trkPy, float[4]);
DECLARE_SOA_COLUMN(TrkPz, trkPz, float[4]);
DECLARE_SOA_COLUMN(TrkSign, trkSign, int[4]);
DECLARE_SOA_COLUMN(TrkDCAxy, trkDCAxy, float[4]);
DECLARE_SOA_COLUMN(TrkDCAz, trkDCAz, float[4]);
DECLARE_SOA_COLUMN(TrkTPCcr, trkTPCcr, int[4]);
DECLARE_SOA_COLUMN(TrkTPCfind, trkTPCfind, int[4]);
DECLARE_SOA_COLUMN(TrkTPCchi2, trkTPCchi2, float[4]);
DECLARE_SOA_COLUMN(TrkITSchi2, trkITSchi2, float[4]);
DECLARE_SOA_COLUMN(TrkITScl, trkITScl, int[4]);
// PID - 4 tracks
DECLARE_SOA_COLUMN(TrkTPCsignal, trkTPCsignal, float[4]);
DECLARE_SOA_COLUMN(TrkTPCnSigmaEl, trkTPCnSigmaEl, float[4]);
DECLARE_SOA_COLUMN(TrkTPCnSigmaPi, trkTPCnSigmaPi, float[4]);
DECLARE_SOA_COLUMN(TrkTPCnSigmaKa, trkTPCnSigmaKa, float[4]);
DECLARE_SOA_COLUMN(TrkTPCnSigmaPr, trkTPCnSigmaPr, float[4]);
DECLARE_SOA_COLUMN(TrkTPCnSigmaMu, trkTPCnSigmaMu, float[4]);
DECLARE_SOA_COLUMN(TrkTOFbeta, trkTOFbeta, float[4]);
DECLARE_SOA_COLUMN(TrkTOFnSigmaEl, trkTOFnSigmaEl, float[4]);
DECLARE_SOA_COLUMN(TrkTOFnSigmaPi, trkTOFnSigmaPi, float[4]);
DECLARE_SOA_COLUMN(TrkTOFnSigmaKa, trkTOFnSigmaKa, float[4]);
DECLARE_SOA_COLUMN(TrkTOFnSigmaPr, trkTOFnSigmaPr, float[4]);
DECLARE_SOA_COLUMN(TrkTOFnSigmaMu, trkTOFnSigmaMu, float[4]);
DECLARE_SOA_COLUMN(TrkTOFchi2, trkTOFchi2, float[4]);

// tracks - 6 tracks
DECLARE_SOA_COLUMN(Trk6Px, trk6Px, float[6]);
DECLARE_SOA_COLUMN(Trk6Py, trk6Py, float[6]);
DECLARE_SOA_COLUMN(Trk6Pz, trk6Pz, float[6]);
DECLARE_SOA_COLUMN(Trk6Sign, trk6Sign, int[6]);
DECLARE_SOA_COLUMN(Trk6DCAxy, trk6DCAxy, float[6]);
DECLARE_SOA_COLUMN(Trk6DCAz, trk6DCAz, float[6]);
DECLARE_SOA_COLUMN(Trk6TPCcr, trk6TPCcr, int[6]);
DECLARE_SOA_COLUMN(Trk6TPCfind, trk6TPCfind, int[6]);
DECLARE_SOA_COLUMN(Trk6TPCchi2, trk6TPCchi2, float[6]);
DECLARE_SOA_COLUMN(Trk6ITSchi2, trk6ITSchi2, float[6]);
DECLARE_SOA_COLUMN(Trk6ITScl, trk6ITScl, int[6]);
// PID - 6 tracks
DECLARE_SOA_COLUMN(Trk6TPCsignal, trk6TPCsignal, float[6]);
DECLARE_SOA_COLUMN(Trk6TPCnSigmaEl, trk6TPCnSigmaEl, float[6]);
DECLARE_SOA_COLUMN(Trk6TPCnSigmaPi, trk6TPCnSigmaPi, float[6]);
DECLARE_SOA_COLUMN(Trk6TPCnSigmaKa, trk6TPCnSigmaKa, float[6]);
DECLARE_SOA_COLUMN(Trk6TPCnSigmaPr, trk6TPCnSigmaPr, float[6]);
DECLARE_SOA_COLUMN(Trk6TPCnSigmaMu, trk6TPCnSigmaMu, float[6]);
DECLARE_SOA_COLUMN(Trk6TOFbeta, trk6TOFbeta, float[6]);
DECLARE_SOA_COLUMN(Trk6TOFnSigmaEl, trk6TOFnSigmaEl, float[6]);
DECLARE_SOA_COLUMN(Trk6TOFnSigmaPi, trk6TOFnSigmaPi, float[6]);
DECLARE_SOA_COLUMN(Trk6TOFnSigmaKa, trk6TOFnSigmaKa, float[6]);
DECLARE_SOA_COLUMN(Trk6TOFnSigmaPr, trk6TOFnSigmaPr, float[6]);
DECLARE_SOA_COLUMN(Trk6TOFnSigmaMu, trk6TOFnSigmaMu, float[6]);
DECLARE_SOA_COLUMN(Trk6TOFchi2, trk6TOFchi2, float[6]);

// truth event
DECLARE_SOA_COLUMN(TrueChannel, trueChannel, int);
// DECLARE_SOA_COLUMN(TrueHasRecoColl, trueHasRecoColl, bool);
// DECLARE_SOA_COLUMN(TruePosX, truePosX, float);
// DECLARE_SOA_COLUMN(TruePosY, truePosY, float);
DECLARE_SOA_COLUMN(TruePosZ, truePosZ, float);
// // truth tau particles // index 0: tau+ // index 1: tau -
// DECLARE_SOA_COLUMN(TrueTauPx, trueTauPx, float[2]);
// DECLARE_SOA_COLUMN(TrueTauPy, trueTauPy, float[2]);
// DECLARE_SOA_COLUMN(TrueTauPz, trueTauPz, float[2]);
// truth tau daughter particles - 4 particles
DECLARE_SOA_COLUMN(TrueDaugPx, trueDaugPx, float[4]);
DECLARE_SOA_COLUMN(TrueDaugPy, trueDaugPy, float[4]);
DECLARE_SOA_COLUMN(TrueDaugPz, trueDaugPz, float[4]);
DECLARE_SOA_COLUMN(TrueDaugPdgCode, trueDaugPdgCode, int[4]);
// truth tau daughter particles - 6 particles
DECLARE_SOA_COLUMN(True6DaugPx, true6DaugPx, float[6]);
DECLARE_SOA_COLUMN(True6DaugPy, true6DaugPy, float[6]);
DECLARE_SOA_COLUMN(True6DaugPz, true6DaugPz, float[6]);
DECLARE_SOA_COLUMN(True6DaugPdgCode, true6DaugPdgCode, int[6]);
DECLARE_SOA_COLUMN(Problem, problem, int8_t);
DECLARE_SOA_COLUMN(IsRec, isRec, int8_t);
} // namespace tautree

DECLARE_SOA_TABLE(DataTauFourTracks, "AOD", "TAUFOURTRACK",
                  tautree::RunNumber, tautree::Bc, tautree::TotalTracks, tautree::NumContrib,
                  tautree::RctOk,
                  // tautree::GlobalNonPVtracks,
                  // tautree::PosX, tautree::PosY,
                  tautree::PosZ,
                  tautree::FlagUPC, tautree::OccupancyInTime, tautree::HadronicRate,
                  //
                  tautree::BcSelBits,
                  // tautree::Trs, tautree::Trofs, tautree::Hmpr,
                  // tautree::Tfb, tautree::ItsRofb, tautree::Sbp, tautree::ZvtxFT0vsPv, tautree::VtxITSTPC,
                  //
                  tautree::ZdcAenergy, tautree::ZdcCenergy,
                  tautree::ZdcAtime, tautree::ZdcCtime,
                  // tautree::Qtot,
                  tautree::TotalFT0AmplitudeA, tautree::TotalFT0AmplitudeC, tautree::TotalFV0AmplitudeA,
                  // tautree::TimeFT0A, tautree::TimeFT0C, tautree::TimeFV0A,
                  tautree::TrkPx, tautree::TrkPy, tautree::TrkPz,
                  tautree::TrkSign,
                  tautree::TrkDCAxy, tautree::TrkDCAz,
                  tautree::TrkTPCcr,
                  tautree::TrkTPCfind, tautree::TrkTPCchi2, tautree::TrkITSchi2, tautree::TrkITScl,
                  tautree::TrkTPCsignal, tautree::TrkTPCnSigmaEl, tautree::TrkTPCnSigmaPi, tautree::TrkTPCnSigmaKa, tautree::TrkTPCnSigmaPr, tautree::TrkTPCnSigmaMu,
                  tautree::TrkTOFbeta, tautree::TrkTOFnSigmaEl, tautree::TrkTOFnSigmaPi, tautree::TrkTOFnSigmaKa, tautree::TrkTOFnSigmaPr, tautree::TrkTOFnSigmaMu,
                  tautree::TrkTOFchi2);

DECLARE_SOA_TABLE(TrueTauFourTracks, "AOD", "TRUETAU",
                  tautree::RunNumber, tautree::Bc, tautree::TotalTracks, tautree::NumContrib,
                  tautree::RctOk,
                  // tautree::GlobalNonPVtracks,
                  // tautree::PosX, tautree::PosY,
                  tautree::PosZ,
                  tautree::FlagUPC, tautree::OccupancyInTime, tautree::HadronicRate,
                  //
                  tautree::BcSelBits,
                  // tautree::Trs, tautree::Trofs, tautree::Hmpr,
                  // tautree::Tfb, tautree::ItsRofb, tautree::Sbp, tautree::ZvtxFT0vsPv, tautree::VtxITSTPC,
                  //
                  // zdc information do not exist in MC
                  // tautree::ZdcAenergy, tautree::ZdcCenergy,
                  // tautree::ZdcAtime, tautree::ZdcCtime,
                  // tautree::Qtot,
                  tautree::TotalFT0AmplitudeA, tautree::TotalFT0AmplitudeC, tautree::TotalFV0AmplitudeA,
                  // tautree::TimeFT0A, tautree::TimeFT0C, tautree::TimeFV0A,
                  tautree::TrkPx, tautree::TrkPy, tautree::TrkPz,
                  tautree::TrkSign,
                  tautree::TrkDCAxy, tautree::TrkDCAz,
                  tautree::TrkTPCcr,
                  tautree::TrkTPCfind, tautree::TrkTPCchi2, tautree::TrkITSchi2, tautree::TrkITScl,
                  tautree::TrkTPCsignal, tautree::TrkTPCnSigmaEl, tautree::TrkTPCnSigmaPi, tautree::TrkTPCnSigmaKa, tautree::TrkTPCnSigmaPr, tautree::TrkTPCnSigmaMu,
                  tautree::TrkTOFbeta, tautree::TrkTOFnSigmaEl, tautree::TrkTOFnSigmaPi, tautree::TrkTOFnSigmaKa, tautree::TrkTOFnSigmaPr, tautree::TrkTOFnSigmaMu,
                  tautree::TrkTOFchi2,
                  tautree::TrueChannel,
                  // tautree::TrueHasRecoColl,
                  tautree::TruePosZ,
                  // tautree::TrueTauPx, tautree::TrueTauPy, tautree::TrueTauPz,
                  tautree::TrueDaugPx, tautree::TrueDaugPy, tautree::TrueDaugPz,
                  tautree::TrueDaugPdgCode,
                  tautree::Problem);

DECLARE_SOA_TABLE(GenTauFourTracks, "AOD", "GENTAU",
                  tautree::TrueChannel,
                  tautree::TruePosZ,
                  // tautree::TrueTauPx, tautree::TrueTauPy, tautree::TrueTauPz,
                  tautree::TrueDaugPx, tautree::TrueDaugPy, tautree::TrueDaugPz,
                  tautree::TrueDaugPdgCode,
                  tautree::Problem,
                  tautree::IsRec);

DECLARE_SOA_TABLE(DataTauSixTracks, "AOD", "TAUSIXTRACK",
                  tautree::RunNumber, tautree::Bc, tautree::TotalTracks, tautree::NumContrib,
                  tautree::RctOk,
                  // tautree::GlobalNonPVtracks,
                  // tautree::PosX, tautree::PosY,
                  tautree::PosZ,
                  tautree::FlagUPC, tautree::OccupancyInTime, tautree::HadronicRate,
                  //
                  tautree::BcSelBits,
                  // tautree::Trs, tautree::Trofs, tautree::Hmpr,
                  // tautree::Tfb, tautree::ItsRofb, tautree::Sbp, tautree::ZvtxFT0vsPv, tautree::VtxITSTPC,
                  //
                  tautree::ZdcAenergy, tautree::ZdcCenergy,
                  tautree::ZdcAtime, tautree::ZdcCtime,
                  // tautree::Qtot,
                  tautree::TotalFT0AmplitudeA, tautree::TotalFT0AmplitudeC, tautree::TotalFV0AmplitudeA,
                  // tautree::TimeFT0A, tautree::TimeFT0C, tautree::TimeFV0A,
                  tautree::Trk6Px, tautree::Trk6Py, tautree::Trk6Pz,
                  tautree::Trk6Sign,
                  tautree::Trk6DCAxy, tautree::Trk6DCAz,
                  tautree::Trk6TPCcr,
                  tautree::Trk6TPCfind, tautree::Trk6TPCchi2, tautree::Trk6ITSchi2, tautree::Trk6ITScl,
                  tautree::Trk6TPCsignal, tautree::Trk6TPCnSigmaEl, tautree::Trk6TPCnSigmaPi, tautree::Trk6TPCnSigmaKa, tautree::Trk6TPCnSigmaPr, tautree::Trk6TPCnSigmaMu,
                  tautree::Trk6TOFbeta, tautree::Trk6TOFnSigmaEl, tautree::Trk6TOFnSigmaPi, tautree::Trk6TOFnSigmaKa, tautree::Trk6TOFnSigmaPr, tautree::Trk6TOFnSigmaMu,
                  tautree::Trk6TOFchi2);

DECLARE_SOA_TABLE(TrueTauSixTracks, "AOD", "TRUETAUSIX",
                  tautree::RunNumber, tautree::Bc, tautree::TotalTracks, tautree::NumContrib,
                  tautree::RctOk,
                  // tautree::GlobalNonPVtracks,
                  // tautree::PosX, tautree::PosY,
                  tautree::PosZ,
                  tautree::FlagUPC, tautree::OccupancyInTime, tautree::HadronicRate,
                  //
                  tautree::BcSelBits,
                  // tautree::Trs, tautree::Trofs, tautree::Hmpr,
                  // tautree::Tfb, tautree::ItsRofb, tautree::Sbp, tautree::ZvtxFT0vsPv, tautree::VtxITSTPC,
                  //
                  // ZDC information do not exist in MC
                  // tautree::ZdcAenergy, tautree::ZdcCenergy,
                  // tautree::ZdcAtime, tautree::ZdcCtime,
                  // tautree::Qtot,
                  tautree::TotalFT0AmplitudeA, tautree::TotalFT0AmplitudeC, tautree::TotalFV0AmplitudeA,
                  // tautree::TimeFT0A, tautree::TimeFT0C, tautree::TimeFV0A,
                  tautree::Trk6Px, tautree::Trk6Py, tautree::Trk6Pz,
                  tautree::Trk6Sign,
                  tautree::Trk6DCAxy, tautree::Trk6DCAz,
                  tautree::Trk6TPCcr,
                  tautree::Trk6TPCfind, tautree::Trk6TPCchi2, tautree::Trk6ITSchi2, tautree::Trk6ITScl,
                  tautree::Trk6TPCsignal, tautree::Trk6TPCnSigmaEl, tautree::Trk6TPCnSigmaPi, tautree::Trk6TPCnSigmaKa, tautree::Trk6TPCnSigmaPr, tautree::Trk6TPCnSigmaMu,
                  tautree::Trk6TOFbeta, tautree::Trk6TOFnSigmaEl, tautree::Trk6TOFnSigmaPi, tautree::Trk6TOFnSigmaKa, tautree::Trk6TOFnSigmaPr, tautree::Trk6TOFnSigmaMu,
                  tautree::Trk6TOFchi2,
                  tautree::TrueChannel,
                  // tautree::TrueHasRecoColl,
                  tautree::TruePosZ,
                  // tautree::TrueTauPx, tautree::TrueTauPy, tautree::TrueTauPz,
                  tautree::True6DaugPx, tautree::True6DaugPy, tautree::True6DaugPz,
                  tautree::True6DaugPdgCode,
                  tautree::Problem);

DECLARE_SOA_TABLE(GenTauSixTracks, "AOD", "GENTAUSIX",
                  tautree::TrueChannel,
                  tautree::TruePosZ,
                  // tautree::TrueTauPx, tautree::TrueTauPy, tautree::TrueTauPz,
                  tautree::True6DaugPx, tautree::True6DaugPy, tautree::True6DaugPz,
                  tautree::True6DaugPdgCode,
                  tautree::Problem,
                  tautree::IsRec);

} // namespace o2::aod

#endif // PWGUD_DATAMODEL_TAUTHREEPRONGEVENTTABLES_H_
