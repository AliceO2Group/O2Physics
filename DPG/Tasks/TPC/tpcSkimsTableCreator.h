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

/// \file tpcSkimsTableCreator.h
/// \author Annalena Kalteyer <annalena.sophie.kalteyer@cern.ch>
/// \author Christian Sonnabend <christian.sonnabend@cern.ch>
/// \author Jeremy Wilkinson <jeremy.wilkinson@cern.ch>
/// \author Ana Marin <ana.marin@cern.ch>
/// \author Oleksii Lubynets <oleksii.lubynets@cern.ch>
/// \brief  Creates clean samples of particles for PID fits

#ifndef DPG_TASKS_TPC_TPCSKIMSTABLECREATOR_H_
#define DPG_TASKS_TPC_TPCSKIMSTABLECREATOR_H_

#include "Common/Core/trackUtilities.h"

#include <Framework/AnalysisDataModel.h>
#include <ReconstructionDataFormats/PID.h>

#include <cstdint>

namespace o2::aod
{
namespace tpcskims
{
DECLARE_SOA_COLUMN(InvDeDxExpTPC, invDeDxExpTPC, float);
DECLARE_SOA_COLUMN(Mass, mass, float);
DECLARE_SOA_COLUMN(BetaGamma, betaGamma, float);
DECLARE_SOA_COLUMN(NormMultTPC, normMultTPC, float);
DECLARE_SOA_COLUMN(NormNClustersTPC, normNClustersTPC, float);
DECLARE_SOA_COLUMN(NormNClustersTPCPID, normNClustersTPCPID, float);
DECLARE_SOA_COLUMN(PidIndex, pidIndex, uint8_t);
DECLARE_SOA_COLUMN(NSigTPC, nSigTPC, float);
DECLARE_SOA_COLUMN(NSigTOF, nSigTOF, float);
DECLARE_SOA_COLUMN(NSigITS, nSigITS, float);
DECLARE_SOA_COLUMN(AlphaV0, alphaV0, float);
DECLARE_SOA_COLUMN(QtV0, qtV0, float);
DECLARE_SOA_COLUMN(CosPAV0, cosPAV0, float);
DECLARE_SOA_COLUMN(PtV0, ptV0, float);
DECLARE_SOA_COLUMN(RadiusV0, radiusV0, float);
DECLARE_SOA_COLUMN(GammaPsiPair, gammaPsiPair, float);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
DECLARE_SOA_COLUMN(TrackOcc, trackOcc, float);
DECLARE_SOA_COLUMN(Ft0Occ, ft0Occ, float);
DECLARE_SOA_COLUMN(OccMedianTime, occMedianTime, float);
DECLARE_SOA_COLUMN(HadronicRate, hadronicRate, float);
DECLARE_SOA_COLUMN(BcGlobalIndex, bcGlobalIndex, int);
DECLARE_SOA_COLUMN(BcTimeFrameId, bcTimeFrameId, int);
DECLARE_SOA_COLUMN(BcBcInTimeFrame, bcBcInTimeFrame, int);
} // namespace tpcskims

#define TPCSKIMS_COLUMNS_BASE      \
  o2::aod::track::TPCSignal,       \
    tpcskims::InvDeDxExpTPC,       \
    o2::aod::track::TPCInnerParam, \
    o2::aod::track::Tgl,           \
    o2::aod::track::Signed1Pt,     \
    o2::aod::track::Eta,           \
    o2::aod::track::Phi,           \
    o2::aod::track::Y,             \
    tpcskims::Mass,                \
    tpcskims::BetaGamma,           \
    tpcskims::NormMultTPC,         \
    tpcskims::NormNClustersTPC,    \
    tpcskims::NormNClustersTPCPID, \
    tpcskims::PidIndex,            \
    tpcskims::NSigTPC,             \
    tpcskims::NSigTOF,             \
    tpcskims::NSigITS,             \
    tpcskims::RunNumber,           \
    tpcskims::TrackOcc,            \
    tpcskims::Ft0Occ,              \
    tpcskims::OccMedianTime,       \
    tpcskims::HadronicRate,        \
    o2::aod::trackqa::TPCdEdxNorm

#define TPCSKIMS_COLUMNS_V0 \
  TPCSKIMS_COLUMNS_BASE,    \
    tpcskims::AlphaV0,      \
    tpcskims::QtV0,         \
    tpcskims::CosPAV0,      \
    tpcskims::PtV0,         \
    tpcskims::RadiusV0,     \
    tpcskims::GammaPsiPair

#define TPCSKIMS_COLUMNS_TOF \
  TPCSKIMS_COLUMNS_BASE

#define TPCSKIMS_COLUMNS_TRACK_QA                   \
  tpcskims::BcGlobalIndex,                          \
    tpcskims::BcTimeFrameId,                        \
    tpcskims::BcBcInTimeFrame,                      \
    o2::aod::trackqa::TPCClusterByteMask,           \
    o2::aod::trackqa::TPCdEdxMax0R,                 \
    o2::aod::trackqa::TPCdEdxMax1R,                 \
    o2::aod::trackqa::TPCdEdxMax2R,                 \
    o2::aod::trackqa::TPCdEdxMax3R,                 \
    o2::aod::trackqa::TPCdEdxTot0R,                 \
    o2::aod::trackqa::TPCdEdxTot1R,                 \
    o2::aod::trackqa::TPCdEdxTot2R,                 \
    o2::aod::trackqa::TPCdEdxTot3R,                 \
    o2::aod::trackmeanocc::TmoPrimUnfm80,           \
    o2::aod::trackmeanocc::TmoFV0AUnfm80,           \
    o2::aod::trackmeanocc::TmoFT0AUnfm80,           \
    o2::aod::trackmeanocc::TmoFT0CUnfm80,           \
    o2::aod::trackmeanocc::TmoRobustT0V0PrimUnfm80, \
    o2::aod::trackmeanocc::TwmoPrimUnfm80,          \
    o2::aod::trackmeanocc::TwmoFV0AUnfm80,          \
    o2::aod::trackmeanocc::TwmoFT0AUnfm80,          \
    o2::aod::trackmeanocc::TwmoFT0CUnfm80,          \
    o2::aod::trackmeanocc::TwmoRobustT0V0PrimUnfm80

DECLARE_SOA_TABLE(SkimmedTPCV0Tree, "AOD", "TPCSKIMV0TREE",
                  TPCSKIMS_COLUMNS_V0);

DECLARE_SOA_TABLE(SkimmedTPCV0TreeWithTrkQA, "AOD", "TPCSKIMV0WQA",
                  TPCSKIMS_COLUMNS_V0,
                  TPCSKIMS_COLUMNS_TRACK_QA);

DECLARE_SOA_TABLE(SkimmedTPCTOFTree, "AOD", "TPCTOFSKIMTREE",
                  TPCSKIMS_COLUMNS_TOF);

DECLARE_SOA_TABLE(SkimmedTPCTOFTreeWithTrkQA, "AOD", "TPCTOFSKIMWQA",
                  TPCSKIMS_COLUMNS_TOF,
                  TPCSKIMS_COLUMNS_TRACK_QA);

#undef TPCSKIMS_COLUMNS_TRACK_QA
#undef TPCSKIMS_COLUMNS_TOF
#undef TPCSKIMS_COLUMNS_V0
#undef TPCSKIMS_COLUMNS_BASE
} // namespace o2::aod

namespace o2::dpg_tpcskimstablecreator
{
enum {
  ModeStandard = 0,
  ModeWithdEdxTrkQA,
  ModeWithTrkQA
};

constexpr o2::track::PID::ID PidElectron{o2::track::PID::Electron};
constexpr o2::track::PID::ID PidPion{o2::track::PID::Pion};
constexpr o2::track::PID::ID PidKaon{o2::track::PID::Kaon};
constexpr o2::track::PID::ID PidProton{o2::track::PID::Proton};
constexpr o2::track::PID::ID PidDeuteron{o2::track::PID::Deuteron};
constexpr o2::track::PID::ID PidTriton{o2::track::PID::Triton};

constexpr int UndefValueInt{-999};
constexpr float UndefValueFloat{-999.f};
constexpr double UndefValueDouble{-999.};

// an arbitrary big value to convert multiplicity into a value between 0 and 1
constexpr double MultiplicityNorm{11000.};

constexpr float OneToKilo{1e-3f};
} // namespace o2::dpg_tpcskimstablecreator
#endif // DPG_TASKS_TPC_TPCSKIMSTABLECREATOR_H_
