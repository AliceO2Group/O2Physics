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
/// \file   TwoFwdTracksEventTables.h
/// \author Roman Lavička
/// \since  2026-10-02
/// \brief  A table to store information about forward UPC candidates preselected to have exactly two muon tracks.
/// \brief  Forward counterpart of TwoTracksEventTables.h, filled from the output of upcCandProducerGlobalMuon
/// \brief  (MCH/MID/MFT tracks, FV0A and ZDC instead of ITS/TPC/TOF tracks, FT0 and event-selection bits).
///

#ifndef PWGUD_DATAMODEL_TWOFWDTRACKSEVENTTABLES_H_
#define PWGUD_DATAMODEL_TWOFWDTRACKSEVENTTABLES_H_

#include <Framework/AnalysisDataModel.h>

#include <cstdint>
#include <vector>

namespace o2::aod
{
namespace two_fwd_tracks_tree
{
// event info
DECLARE_SOA_COLUMN(RunNumber, runNumber, int32_t);
DECLARE_SOA_COLUMN(Bc, bc, uint64_t);
DECLARE_SOA_COLUMN(NumContrib, numContrib, int);
DECLARE_SOA_COLUMN(PosX, posX, float);
DECLARE_SOA_COLUMN(PosY, posY, float);
DECLARE_SOA_COLUMN(PosZ, posZ, float);
// FV0A and ZDC info
DECLARE_SOA_COLUMN(TotalFV0AmplitudeA, totalFV0AmplitudeA, float);
DECLARE_SOA_COLUMN(AmplitudesV0A, amplitudesV0A, std::vector<float>);
DECLARE_SOA_COLUMN(AmpRelBCsV0A, ampRelBCsV0A, std::vector<int8_t>);
DECLARE_SOA_COLUMN(EnergyCommonZNA, energyCommonZNA, float);
DECLARE_SOA_COLUMN(EnergyCommonZNC, energyCommonZNC, float);
DECLARE_SOA_COLUMN(TimeZNA, timeZNA, float);
DECLARE_SOA_COLUMN(TimeZNC, timeZNC, float);
// tracks
DECLARE_SOA_COLUMN(TrkPx, trkPx, float[2]);
DECLARE_SOA_COLUMN(TrkPy, trkPy, float[2]);
DECLARE_SOA_COLUMN(TrkPz, trkPz, float[2]);
DECLARE_SOA_COLUMN(TrkSign, trkSign, int[2]);
DECLARE_SOA_COLUMN(TrkTime, trkTime, float[2]);
DECLARE_SOA_COLUMN(TrkTimeRes, trkTimeRes, float[2]);
DECLARE_SOA_COLUMN(TrkType, trkType, int[2]);
DECLARE_SOA_COLUMN(TrkNClusters, trkNClusters, int[2]);
DECLARE_SOA_COLUMN(TrkPDca, trkPDca, float[2]);
DECLARE_SOA_COLUMN(TrkRAtAbsorberEnd, trkRAtAbsorberEnd, float[2]);
DECLARE_SOA_COLUMN(TrkChi2, trkChi2, float[2]);
DECLARE_SOA_COLUMN(TrkChi2MatchMCHMID, trkChi2MatchMCHMID, float[2]);
DECLARE_SOA_COLUMN(TrkChi2MatchMCHMFT, trkChi2MatchMCHMFT, float[2]);
DECLARE_SOA_COLUMN(TrkMCHBitMap, trkMCHBitMap, int[2]);
DECLARE_SOA_COLUMN(TrkMIDBitMap, trkMIDBitMap, int[2]);
DECLARE_SOA_COLUMN(Trk1MIDBoards, trk1MIDBoards, uint32_t);
DECLARE_SOA_COLUMN(Trk2MIDBoards, trk2MIDBoards, uint32_t);
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

} // namespace two_fwd_tracks_tree

DECLARE_SOA_TABLE(TwoFwdTracks, "AOD", "TWOFWDTRACK",
                  two_fwd_tracks_tree::RunNumber,
                  two_fwd_tracks_tree::Bc,
                  two_fwd_tracks_tree::NumContrib,
                  two_fwd_tracks_tree::PosX,
                  two_fwd_tracks_tree::PosY,
                  two_fwd_tracks_tree::PosZ,
                  two_fwd_tracks_tree::TotalFV0AmplitudeA,
                  two_fwd_tracks_tree::AmplitudesV0A,
                  two_fwd_tracks_tree::AmpRelBCsV0A,
                  two_fwd_tracks_tree::EnergyCommonZNA,
                  two_fwd_tracks_tree::EnergyCommonZNC,
                  two_fwd_tracks_tree::TimeZNA,
                  two_fwd_tracks_tree::TimeZNC,
                  two_fwd_tracks_tree::TrkPx,
                  two_fwd_tracks_tree::TrkPy,
                  two_fwd_tracks_tree::TrkPz,
                  two_fwd_tracks_tree::TrkSign,
                  two_fwd_tracks_tree::TrkTime,
                  two_fwd_tracks_tree::TrkTimeRes,
                  two_fwd_tracks_tree::TrkType,
                  two_fwd_tracks_tree::TrkNClusters,
                  two_fwd_tracks_tree::TrkPDca,
                  two_fwd_tracks_tree::TrkRAtAbsorberEnd,
                  two_fwd_tracks_tree::TrkChi2,
                  two_fwd_tracks_tree::TrkChi2MatchMCHMID,
                  two_fwd_tracks_tree::TrkChi2MatchMCHMFT,
                  two_fwd_tracks_tree::TrkMCHBitMap,
                  two_fwd_tracks_tree::TrkMIDBitMap,
                  two_fwd_tracks_tree::Trk1MIDBoards,
                  two_fwd_tracks_tree::Trk2MIDBoards);

DECLARE_SOA_TABLE(TrueTwoFwdTracks, "AOD", "TRUETWOFWDTRACK",
                  two_fwd_tracks_tree::RunNumber,
                  two_fwd_tracks_tree::Bc,
                  two_fwd_tracks_tree::NumContrib,
                  two_fwd_tracks_tree::PosX,
                  two_fwd_tracks_tree::PosY,
                  two_fwd_tracks_tree::PosZ,
                  two_fwd_tracks_tree::TotalFV0AmplitudeA,
                  two_fwd_tracks_tree::AmplitudesV0A,
                  two_fwd_tracks_tree::AmpRelBCsV0A,
                  two_fwd_tracks_tree::EnergyCommonZNA,
                  two_fwd_tracks_tree::EnergyCommonZNC,
                  two_fwd_tracks_tree::TimeZNA,
                  two_fwd_tracks_tree::TimeZNC,
                  two_fwd_tracks_tree::TrkPx,
                  two_fwd_tracks_tree::TrkPy,
                  two_fwd_tracks_tree::TrkPz,
                  two_fwd_tracks_tree::TrkSign,
                  two_fwd_tracks_tree::TrkTime,
                  two_fwd_tracks_tree::TrkTimeRes,
                  two_fwd_tracks_tree::TrkType,
                  two_fwd_tracks_tree::TrkNClusters,
                  two_fwd_tracks_tree::TrkPDca,
                  two_fwd_tracks_tree::TrkRAtAbsorberEnd,
                  two_fwd_tracks_tree::TrkChi2,
                  two_fwd_tracks_tree::TrkChi2MatchMCHMID,
                  two_fwd_tracks_tree::TrkChi2MatchMCHMFT,
                  two_fwd_tracks_tree::TrkMCHBitMap,
                  two_fwd_tracks_tree::TrkMIDBitMap,
                  two_fwd_tracks_tree::Trk1MIDBoards,
                  two_fwd_tracks_tree::Trk2MIDBoards,
                  two_fwd_tracks_tree::TrueChannel,
                  two_fwd_tracks_tree::TrueHasRecoColl,
                  two_fwd_tracks_tree::TruePosX,
                  two_fwd_tracks_tree::TruePosY,
                  two_fwd_tracks_tree::TruePosZ,
                  two_fwd_tracks_tree::TrueMotherPx,
                  two_fwd_tracks_tree::TrueMotherPy,
                  two_fwd_tracks_tree::TrueMotherPz,
                  two_fwd_tracks_tree::TrueDaugPx,
                  two_fwd_tracks_tree::TrueDaugPy,
                  two_fwd_tracks_tree::TrueDaugPz,
                  two_fwd_tracks_tree::TrueDaugPdgCode,
                  two_fwd_tracks_tree::ProblematicEvent);

} // namespace o2::aod

#endif // PWGUD_DATAMODEL_TWOFWDTRACKSEVENTTABLES_H_
