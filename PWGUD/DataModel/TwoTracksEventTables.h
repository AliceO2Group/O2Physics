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
/// \file   TwoTracksEventTables.h
/// \author Roman Lavička
/// \since  2025-06-20
/// \brief  A table to store information about events preselected to have exactly two tracks.
/// \brief  Good for UPC gammagamma (->elel,mumu,tautau) and gammalead (vector mesons)
/// \brief  If MC, careful with filling the mother
///

#ifndef PWGUD_DATAMODEL_TWOTRACKSEVENTTABLES_H_
#define PWGUD_DATAMODEL_TWOTRACKSEVENTTABLES_H_

#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace two_tracks_tree
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
DECLARE_SOA_COLUMN(TrueMotherPx, trueMotherPx, float[2]);
DECLARE_SOA_COLUMN(TrueMotherPy, trueMotherPy, float[2]);
DECLARE_SOA_COLUMN(TrueMotherPz, trueMotherPz, float[2]);
DECLARE_SOA_COLUMN(TrueDaugPx, trueDaugPx, float[2]);
DECLARE_SOA_COLUMN(TrueDaugPy, trueDaugPy, float[2]);
DECLARE_SOA_COLUMN(TrueDaugPz, trueDaugPz, float[2]);
DECLARE_SOA_COLUMN(TrueDaugPdgCode, trueDaugPdgCode, int[2]);
// additional info
DECLARE_SOA_COLUMN(ProblematicEvent, problematicEvent, bool);

} // namespace two_tracks_tree
DECLARE_SOA_TABLE(TwoTracks, "AOD", "TWOTRACK",
                  two_tracks_tree::RunNumber,
                  two_tracks_tree::Bc,
                  two_tracks_tree::TotalTracks,
                  two_tracks_tree::NumContrib,
                  two_tracks_tree::GlobalNonPVtracks,
                  two_tracks_tree::PosX,
                  two_tracks_tree::PosY,
                  two_tracks_tree::PosZ,
                  two_tracks_tree::RecoMode,
                  two_tracks_tree::OccupancyInTime,
                  two_tracks_tree::HadronicRate,
                  two_tracks_tree::Trs,
                  two_tracks_tree::Trofs,
                  two_tracks_tree::Hmpr,
                  two_tracks_tree::Tfb,
                  two_tracks_tree::ItsRofb,
                  two_tracks_tree::Sbp,
                  two_tracks_tree::ZvtxFT0vsPv,
                  two_tracks_tree::VtxITSTPC,
                  two_tracks_tree::TotalFT0AmplitudeA,
                  two_tracks_tree::TotalFT0AmplitudeC,
                  two_tracks_tree::TotalFV0AmplitudeA,
                  two_tracks_tree::EnergyCommonZNA,
                  two_tracks_tree::EnergyCommonZNC,
                  two_tracks_tree::TimeFT0A,
                  two_tracks_tree::TimeFT0C,
                  two_tracks_tree::TimeFV0A,
                  two_tracks_tree::TimeZNA,
                  two_tracks_tree::TimeZNC,
                  two_tracks_tree::TrkPx,
                  two_tracks_tree::TrkPy,
                  two_tracks_tree::TrkPz,
                  two_tracks_tree::TrkSign,
                  two_tracks_tree::TrkDCAxy,
                  two_tracks_tree::TrkDCAz,
                  two_tracks_tree::TrkTimeRes,
                  two_tracks_tree::Trk1ITSclusterSizes,
                  two_tracks_tree::Trk2ITSclusterSizes,
                  two_tracks_tree::TrkTPCsignal,
                  two_tracks_tree::TrkTPCnSigmaEl,
                  two_tracks_tree::TrkTPCnSigmaMu,
                  two_tracks_tree::TrkTPCnSigmaPi,
                  two_tracks_tree::TrkTPCnSigmaKa,
                  two_tracks_tree::TrkTPCnSigmaPr,
                  two_tracks_tree::TrkTPCinnerParam,
                  two_tracks_tree::TrkTOFsignal,
                  two_tracks_tree::TrkTOFnSigmaEl,
                  two_tracks_tree::TrkTOFnSigmaMu,
                  two_tracks_tree::TrkTOFnSigmaPi,
                  two_tracks_tree::TrkTOFnSigmaKa,
                  two_tracks_tree::TrkTOFnSigmaPr,
                  two_tracks_tree::TrkTOFexpMom);

DECLARE_SOA_TABLE(TrueTwoTracks, "AOD", "TRUETWOTRACK",
                  two_tracks_tree::RunNumber,
                  two_tracks_tree::Bc,
                  two_tracks_tree::TotalTracks,
                  two_tracks_tree::NumContrib,
                  two_tracks_tree::GlobalNonPVtracks,
                  two_tracks_tree::PosX,
                  two_tracks_tree::PosY,
                  two_tracks_tree::PosZ,
                  two_tracks_tree::RecoMode,
                  two_tracks_tree::OccupancyInTime,
                  two_tracks_tree::HadronicRate,
                  two_tracks_tree::Trs,
                  two_tracks_tree::Trofs,
                  two_tracks_tree::Hmpr,
                  two_tracks_tree::Tfb,
                  two_tracks_tree::ItsRofb,
                  two_tracks_tree::Sbp,
                  two_tracks_tree::ZvtxFT0vsPv,
                  two_tracks_tree::VtxITSTPC,
                  two_tracks_tree::TotalFT0AmplitudeA,
                  two_tracks_tree::TotalFT0AmplitudeC,
                  two_tracks_tree::TotalFV0AmplitudeA,
                  two_tracks_tree::EnergyCommonZNA,
                  two_tracks_tree::EnergyCommonZNC,
                  two_tracks_tree::TimeFT0A,
                  two_tracks_tree::TimeFT0C,
                  two_tracks_tree::TimeFV0A,
                  two_tracks_tree::TimeZNA,
                  two_tracks_tree::TimeZNC,
                  two_tracks_tree::TrkPx,
                  two_tracks_tree::TrkPy,
                  two_tracks_tree::TrkPz,
                  two_tracks_tree::TrkSign,
                  two_tracks_tree::TrkDCAxy,
                  two_tracks_tree::TrkDCAz,
                  two_tracks_tree::TrkTimeRes,
                  two_tracks_tree::Trk1ITSclusterSizes,
                  two_tracks_tree::Trk2ITSclusterSizes,
                  two_tracks_tree::TrkTPCsignal,
                  two_tracks_tree::TrkTPCnSigmaEl,
                  two_tracks_tree::TrkTPCnSigmaMu,
                  two_tracks_tree::TrkTPCnSigmaPi,
                  two_tracks_tree::TrkTPCnSigmaKa,
                  two_tracks_tree::TrkTPCnSigmaPr,
                  two_tracks_tree::TrkTPCinnerParam,
                  two_tracks_tree::TrkTOFsignal,
                  two_tracks_tree::TrkTOFnSigmaEl,
                  two_tracks_tree::TrkTOFnSigmaMu,
                  two_tracks_tree::TrkTOFnSigmaPi,
                  two_tracks_tree::TrkTOFnSigmaKa,
                  two_tracks_tree::TrkTOFnSigmaPr,
                  two_tracks_tree::TrkTOFexpMom,
                  two_tracks_tree::TrueChannel,
                  two_tracks_tree::TrueHasRecoColl,
                  two_tracks_tree::TruePosX,
                  two_tracks_tree::TruePosY,
                  two_tracks_tree::TruePosZ,
                  two_tracks_tree::TrueMotherPx,
                  two_tracks_tree::TrueMotherPy,
                  two_tracks_tree::TrueMotherPz,
                  two_tracks_tree::TrueDaugPx,
                  two_tracks_tree::TrueDaugPy,
                  two_tracks_tree::TrueDaugPz,
                  two_tracks_tree::TrueDaugPdgCode,
                  two_tracks_tree::ProblematicEvent);

} // namespace o2::aod

#endif // PWGUD_DATAMODEL_TWOTRACKSEVENTTABLES_H_
