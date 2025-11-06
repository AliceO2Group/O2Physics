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

#include "Framework/AnalysisDataModel.h"
// derived tables for tautau->4 (=1+3) tracks
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
DECLARE_SOA_COLUMN(Trs, trs, int8_t);
DECLARE_SOA_COLUMN(Trofs, trofs, int8_t);
DECLARE_SOA_COLUMN(Hmpr, hmpr, int8_t);
DECLARE_SOA_COLUMN(Tfb, tfb, int8_t);
DECLARE_SOA_COLUMN(ItsRofb, itsRofb, int8_t);
DECLARE_SOA_COLUMN(Sbp, sbp, int8_t);
DECLARE_SOA_COLUMN(ZvtxFT0vsPv, zvtxFT0vsPv, int8_t);
DECLARE_SOA_COLUMN(VtxITSTPC, vtxITSTPC, int8_t);
DECLARE_SOA_COLUMN(ZdcAenergy, zdcAenergy, float);
DECLARE_SOA_COLUMN(ZdcCenergy, zdcCenergy, float);
// DECLARE_SOA_COLUMN(Qtot, qtot, int8_t);
// FIT info
DECLARE_SOA_COLUMN(TotalFT0AmplitudeA, totalFT0AmplitudeA, float);
DECLARE_SOA_COLUMN(TotalFT0AmplitudeC, totalFT0AmplitudeC, float);
DECLARE_SOA_COLUMN(TotalFV0AmplitudeA, totalFV0AmplitudeA, float);
// DECLARE_SOA_COLUMN(TimeFT0A, timeFT0A, float);
// DECLARE_SOA_COLUMN(TimeFT0C, timeFT0C, float);
// DECLARE_SOA_COLUMN(TimeFV0A, timeFV0A, float);
// tracks
DECLARE_SOA_COLUMN(TrkPx, trkPx, float[6]);
DECLARE_SOA_COLUMN(TrkPy, trkPy, float[6]);
DECLARE_SOA_COLUMN(TrkPz, trkPz, float[6]);
DECLARE_SOA_COLUMN(TrkSign, trkSign, int[6]);
DECLARE_SOA_COLUMN(TrkDCAxy, trkDCAxy, float[6]);
DECLARE_SOA_COLUMN(TrkDCAz, trkDCAz, float[6]);
DECLARE_SOA_COLUMN(TrkTPCcr, trkTPCcr, int[6]);
DECLARE_SOA_COLUMN(TrkTPCfind, trkTPCfind, int[6]);
DECLARE_SOA_COLUMN(TrkTPCchi2, trkTPCchi2, float[6]);
DECLARE_SOA_COLUMN(TrkITSchi2, trkITSchi2, float[6]);
DECLARE_SOA_COLUMN(TrkITScl, trkITScl, int[6]);
// PID
DECLARE_SOA_COLUMN(TrkTPCsignal, trkTPCsignal, float[6]);
DECLARE_SOA_COLUMN(TrkTPCnSigmaEl, trkTPCnSigmaEl, float[6]);
DECLARE_SOA_COLUMN(TrkTPCnSigmaPi, trkTPCnSigmaPi, float[6]);
DECLARE_SOA_COLUMN(TrkTPCnSigmaKa, trkTPCnSigmaKa, float[6]);
DECLARE_SOA_COLUMN(TrkTPCnSigmaPr, trkTPCnSigmaPr, float[6]);
DECLARE_SOA_COLUMN(TrkTPCnSigmaMu, trkTPCnSigmaMu, float[6]);
DECLARE_SOA_COLUMN(TrkTOFbeta, trkTOFbeta, float[6]);
DECLARE_SOA_COLUMN(TrkTOFnSigmaEl, trkTOFnSigmaEl, float[6]);
DECLARE_SOA_COLUMN(TrkTOFnSigmaPi, trkTOFnSigmaPi, float[6]);
DECLARE_SOA_COLUMN(TrkTOFnSigmaKa, trkTOFnSigmaKa, float[6]);
DECLARE_SOA_COLUMN(TrkTOFnSigmaPr, trkTOFnSigmaPr, float[6]);
DECLARE_SOA_COLUMN(TrkTOFnSigmaMu, trkTOFnSigmaMu, float[6]);
DECLARE_SOA_COLUMN(TrkTOFchi2, trkTOFchi2, float[6]);
// truth event
DECLARE_SOA_COLUMN(TrueChannel, trueChannel, int);
DECLARE_SOA_COLUMN(TrueHasRecoColl, trueHasRecoColl, bool);
// DECLARE_SOA_COLUMN(TruePosX, truePosX, float);
// DECLARE_SOA_COLUMN(TruePosY, truePosY, float);
DECLARE_SOA_COLUMN(TruePosZ, truePosZ, float);
// truth tau particles // index 0: tau+ // index 1: tau -
DECLARE_SOA_COLUMN(TrueTauPx, trueTauPx, float[2]);
DECLARE_SOA_COLUMN(TrueTauPy, trueTauPy, float[2]);
DECLARE_SOA_COLUMN(TrueTauPz, trueTauPz, float[2]);
// truth tau daughter particles
DECLARE_SOA_COLUMN(TrueDaugPx, trueDaugPx, float[6]);
DECLARE_SOA_COLUMN(TrueDaugPy, trueDaugPy, float[6]);
DECLARE_SOA_COLUMN(TrueDaugPz, trueDaugPz, float[6]);
DECLARE_SOA_COLUMN(TrueDaugPdgCode, trueDaugPdgCode, int[6]);
DECLARE_SOA_COLUMN(Problem, problem, int8_t);
} // namespace tautree

DECLARE_SOA_TABLE(DataTauFourTracks, "AOD", "TAUFOURTRACK",
                  tautree::RunNumber, tautree::Bc, tautree::TotalTracks, tautree::NumContrib,
                  tautree::RctOk,
                  // tautree::GlobalNonPVtracks,
                  // tautree::PosX, tautree::PosY,
                  tautree::PosZ,
                  tautree::FlagUPC, tautree::OccupancyInTime, tautree::HadronicRate,
                  //
                  tautree::Trs, tautree::Trofs, tautree::Hmpr,
                  tautree::Tfb, tautree::ItsRofb, tautree::Sbp, tautree::ZvtxFT0vsPv, tautree::VtxITSTPC,
                  tautree::ZdcAenergy, tautree::ZdcCenergy,
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
                  tautree::Trs, tautree::Trofs, tautree::Hmpr,
                  tautree::Tfb, tautree::ItsRofb, tautree::Sbp, tautree::ZvtxFT0vsPv, tautree::VtxITSTPC,
                  tautree::ZdcAenergy, tautree::ZdcCenergy,
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
                  tautree::TrueHasRecoColl,
                  tautree::TruePosZ,
                  tautree::TrueTauPx, tautree::TrueTauPy, tautree::TrueTauPz,
                  tautree::TrueDaugPx, tautree::TrueDaugPy, tautree::TrueDaugPz,
                  tautree::TrueDaugPdgCode,
                  tautree::Problem);

} // namespace o2::aod

#endif // PWGUD_DATAMODEL_TAUTHREEPRONGEVENTTABLES_H_
