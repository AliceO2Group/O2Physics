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
/// \file   TauEventTables.h
/// \author Roman Lavička
/// \since  2025-04-23
/// \brief  A table to store information about events preselected to be candidates for UPC gammagamma->tautau
///

#ifndef ALISW_TAUEVENTTABLES_H
#define ALISW_TAUEVENTTABLES_H

#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace tau_tree
{
// event info
DECLARE_SOA_COLUMN(RunNumber, runNumber, int32_t);
DECLARE_SOA_COLUMN(Bc, bc, int);
DECLARE_SOA_COLUMN(TotalTracks, totalTracks, int);
DECLARE_SOA_COLUMN(NumContrib, numContrib, int);
DECLARE_SOA_COLUMN(GlobalNonPVtracks, globalNonPVtracks, int);
DECLARE_SOA_COLUMN(PosX, posX, float);
DECLARE_SOA_COLUMN(PosY, posY, float);
DECLARE_SOA_COLUMN(PosZ, posZ, float);
DECLARE_SOA_COLUMN(RecoMode, recoMode, int);
DECLARE_SOA_COLUMN(OccupancyInTime, occupancyInTime, int);
DECLARE_SOA_COLUMN(HadronicRate, hadronicRate, double);
DECLARE_SOA_COLUMN(Trs, trs, int);
DECLARE_SOA_COLUMN(Trofs, trofs, int);
DECLARE_SOA_COLUMN(Hmpr, hmpr, int);
DECLARE_SOA_COLUMN(Tfb, tfb, int);
DECLARE_SOA_COLUMN(ItsRofb, itsRofb, int);
DECLARE_SOA_COLUMN(Sbp, sbp, int);
DECLARE_SOA_COLUMN(ZvtxFT0vsPv, zvtxFT0vsPv, int);
DECLARE_SOA_COLUMN(VtxITSTPC, vtxITSTPC, int);
// FIT info
DECLARE_SOA_COLUMN(TotalFT0AmplitudeA, totalFT0AmplitudeA, float);
DECLARE_SOA_COLUMN(TotalFT0AmplitudeC, totalFT0AmplitudeC, float);
DECLARE_SOA_COLUMN(TotalFV0AmplitudeA, totalFV0AmplitudeA, float);
DECLARE_SOA_COLUMN(EnergyCommonZNA, energyCommonZNA, float);
DECLARE_SOA_COLUMN(EnergyCommonZNC, energyCommonZNC, float);
DECLARE_SOA_COLUMN(TimeFT0A, timeFT0A, float);
DECLARE_SOA_COLUMN(TimeFT0C, timeFT0C, float);
DECLARE_SOA_COLUMN(TimeFV0A, timeFV0A, float);
DECLARE_SOA_COLUMN(TimeZNA, timeZNA, float);
DECLARE_SOA_COLUMN(TimeZNC, timeZNC, float);
// tracks
DECLARE_SOA_COLUMN(TrkPx, trkPx, float[2]);
DECLARE_SOA_COLUMN(TrkPy, trkPy, float[2]);
DECLARE_SOA_COLUMN(TrkPz, trkPz, float[2]);
DECLARE_SOA_COLUMN(TrkSign, trkSign, int[2]);
DECLARE_SOA_COLUMN(TrkDCAxy, trkDCAxy, float[2]);
DECLARE_SOA_COLUMN(TrkDCAz, trkDCAz, float[2]);
DECLARE_SOA_COLUMN(TrkTimeRes, trkTimeRes, float[2]);
DECLARE_SOA_COLUMN(Trk1ITSclusterSizes, trk1ITSclusterSizes, uint32_t);
DECLARE_SOA_COLUMN(Trk2ITSclusterSizes, trk2ITSclusterSizes, uint32_t);
DECLARE_SOA_COLUMN(TrkTPCsignal, trkTPCsignal, float[2]);
DECLARE_SOA_COLUMN(TrkTPCnSigmaEl, trkTPCnSigmaEl, float[2]);
DECLARE_SOA_COLUMN(TrkTPCnSigmaMu, trkTPCnSigmaMu, float[2]);
DECLARE_SOA_COLUMN(TrkTPCnSigmaPi, trkTPCnSigmaPi, float[2]);
DECLARE_SOA_COLUMN(TrkTPCnSigmaKa, trkTPCnSigmaKa, float[2]);
DECLARE_SOA_COLUMN(TrkTPCnSigmaPr, trkTPCnSigmaPr, float[2]);
DECLARE_SOA_COLUMN(TrkTPCinnerParam, trkTPCinnerParam, float[2]);
DECLARE_SOA_COLUMN(TrkTOFsignal, trkTOFsignal, float[2]);
DECLARE_SOA_COLUMN(TrkTOFnSigmaEl, trkTOFnSigmaEl, float[2]);
DECLARE_SOA_COLUMN(TrkTOFnSigmaMu, trkTOFnSigmaMu, float[2]);
DECLARE_SOA_COLUMN(TrkTOFnSigmaPi, trkTOFnSigmaPi, float[2]);
DECLARE_SOA_COLUMN(TrkTOFnSigmaKa, trkTOFnSigmaKa, float[2]);
DECLARE_SOA_COLUMN(TrkTOFnSigmaPr, trkTOFnSigmaPr, float[2]);
DECLARE_SOA_COLUMN(TrkTOFexpMom, trkTOFexpMom, float[2]);
// truth event
DECLARE_SOA_COLUMN(TrueChannel, trueChannel, int);
DECLARE_SOA_COLUMN(TrueHasRecoColl, trueHasRecoColl, bool);
DECLARE_SOA_COLUMN(TruePosX, truePosX, float);
DECLARE_SOA_COLUMN(TruePosY, truePosY, float);
DECLARE_SOA_COLUMN(TruePosZ, truePosZ, float);
// truth particles
DECLARE_SOA_COLUMN(TrueTauPx, trueTauPx, float[2]);
DECLARE_SOA_COLUMN(TrueTauPy, trueTauPy, float[2]);
DECLARE_SOA_COLUMN(TrueTauPz, trueTauPz, float[2]);
DECLARE_SOA_COLUMN(TrueDaugPx, trueDaugPx, float[2]);
DECLARE_SOA_COLUMN(TrueDaugPy, trueDaugPy, float[2]);
DECLARE_SOA_COLUMN(TrueDaugPz, trueDaugPz, float[2]);
DECLARE_SOA_COLUMN(TrueDaugPdgCode, trueDaugPdgCode, int[2]);
// additional info
DECLARE_SOA_COLUMN(ProblematicEvent, problematicEvent, bool);

} // namespace tau_tree
DECLARE_SOA_TABLE(TauTwoTracks, "AOD", "TAUTWOTRACK",
                  tau_tree::RunNumber,
                  tau_tree::Bc,
                  tau_tree::TotalTracks,
                  tau_tree::NumContrib,
                  tau_tree::GlobalNonPVtracks,
                  tau_tree::PosX,
                  tau_tree::PosY,
                  tau_tree::PosZ,
                  tau_tree::RecoMode,
                  tau_tree::OccupancyInTime,
                  tau_tree::HadronicRate,
                  tau_tree::Trs,
                  tau_tree::Trofs,
                  tau_tree::Hmpr,
                  tau_tree::Tfb,
                  tau_tree::ItsRofb,
                  tau_tree::Sbp,
                  tau_tree::ZvtxFT0vsPv,
                  tau_tree::VtxITSTPC,
                  tau_tree::TotalFT0AmplitudeA,
                  tau_tree::TotalFT0AmplitudeC,
                  tau_tree::TotalFV0AmplitudeA,
                  tau_tree::EnergyCommonZNA,
                  tau_tree::EnergyCommonZNC,
                  tau_tree::TimeFT0A,
                  tau_tree::TimeFT0C,
                  tau_tree::TimeFV0A,
                  tau_tree::TimeZNA,
                  tau_tree::TimeZNC,
                  tau_tree::TrkPx,
                  tau_tree::TrkPy,
                  tau_tree::TrkPz,
                  tau_tree::TrkSign,
                  tau_tree::TrkDCAxy,
                  tau_tree::TrkDCAz,
                  tau_tree::TrkTimeRes,
                  tau_tree::Trk1ITSclusterSizes,
                  tau_tree::Trk2ITSclusterSizes,
                  tau_tree::TrkTPCsignal,
                  tau_tree::TrkTPCnSigmaEl,
                  tau_tree::TrkTPCnSigmaMu,
                  tau_tree::TrkTPCnSigmaPi,
                  tau_tree::TrkTPCnSigmaKa,
                  tau_tree::TrkTPCnSigmaPr,
                  tau_tree::TrkTPCinnerParam,
                  tau_tree::TrkTOFsignal,
                  tau_tree::TrkTOFnSigmaEl,
                  tau_tree::TrkTOFnSigmaMu,
                  tau_tree::TrkTOFnSigmaPi,
                  tau_tree::TrkTOFnSigmaKa,
                  tau_tree::TrkTOFnSigmaPr,
                  tau_tree::TrkTOFexpMom);

DECLARE_SOA_TABLE(TrueTauTwoTracks, "AOD", "TRUETAUTWOTRACK",
                  tau_tree::RunNumber,
                  tau_tree::Bc,
                  tau_tree::TotalTracks,
                  tau_tree::NumContrib,
                  tau_tree::GlobalNonPVtracks,
                  tau_tree::PosX,
                  tau_tree::PosY,
                  tau_tree::PosZ,
                  tau_tree::RecoMode,
                  tau_tree::OccupancyInTime,
                  tau_tree::HadronicRate,
                  tau_tree::Trs,
                  tau_tree::Trofs,
                  tau_tree::Hmpr,
                  tau_tree::Tfb,
                  tau_tree::ItsRofb,
                  tau_tree::Sbp,
                  tau_tree::ZvtxFT0vsPv,
                  tau_tree::VtxITSTPC,
                  tau_tree::TotalFT0AmplitudeA,
                  tau_tree::TotalFT0AmplitudeC,
                  tau_tree::TotalFV0AmplitudeA,
                  tau_tree::EnergyCommonZNA,
                  tau_tree::EnergyCommonZNC,
                  tau_tree::TimeFT0A,
                  tau_tree::TimeFT0C,
                  tau_tree::TimeFV0A,
                  tau_tree::TimeZNA,
                  tau_tree::TimeZNC,
                  tau_tree::TrkPx,
                  tau_tree::TrkPy,
                  tau_tree::TrkPz,
                  tau_tree::TrkSign,
                  tau_tree::TrkDCAxy,
                  tau_tree::TrkDCAz,
                  tau_tree::TrkTimeRes,
                  tau_tree::Trk1ITSclusterSizes,
                  tau_tree::Trk2ITSclusterSizes,
                  tau_tree::TrkTPCsignal,
                  tau_tree::TrkTPCnSigmaEl,
                  tau_tree::TrkTPCnSigmaMu,
                  tau_tree::TrkTPCnSigmaPi,
                  tau_tree::TrkTPCnSigmaKa,
                  tau_tree::TrkTPCnSigmaPr,
                  tau_tree::TrkTPCinnerParam,
                  tau_tree::TrkTOFsignal,
                  tau_tree::TrkTOFnSigmaEl,
                  tau_tree::TrkTOFnSigmaMu,
                  tau_tree::TrkTOFnSigmaPi,
                  tau_tree::TrkTOFnSigmaKa,
                  tau_tree::TrkTOFnSigmaPr,
                  tau_tree::TrkTOFexpMom,
                  tau_tree::TrueChannel,
                  tau_tree::TrueHasRecoColl,
                  tau_tree::TruePosX,
                  tau_tree::TruePosY,
                  tau_tree::TruePosZ,
                  tau_tree::TrueTauPx,
                  tau_tree::TrueTauPy,
                  tau_tree::TrueTauPz,
                  tau_tree::TrueDaugPx,
                  tau_tree::TrueDaugPy,
                  tau_tree::TrueDaugPz,
                  tau_tree::TrueDaugPdgCode,
                  tau_tree::ProblematicEvent);

} // namespace o2::aod

#endif // ALISW_TAUEVENTTABLES_H
