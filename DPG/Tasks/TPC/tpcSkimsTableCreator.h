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

#ifndef DPG_TASKS_TPC_TPCSKIMSTABLECREATOR_H_
#define DPG_TASKS_TPC_TPCSKIMSTABLECREATOR_H_

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/PIDResponse.h"

namespace o2::aod
{
namespace tpcskims
{
DECLARE_SOA_COLUMN(InvDeDxExpTPC, invdEdxExpTPC, float);
DECLARE_SOA_COLUMN(Mass, mass, float);
DECLARE_SOA_COLUMN(BetaGamma, bg, float);
DECLARE_SOA_COLUMN(NormMultTPC, normMultTPC, float);
DECLARE_SOA_COLUMN(NormNClustersTPC, normNClustersTPC, float);
DECLARE_SOA_COLUMN(PidIndex, pidIndexTPC, uint8_t);
DECLARE_SOA_COLUMN(NSigTPC, nsigTPC, float);
DECLARE_SOA_COLUMN(NSigTOF, nsigTOF, float);
DECLARE_SOA_COLUMN(AlphaV0, alphaV0, float);
DECLARE_SOA_COLUMN(QtV0, qtV0, float);
DECLARE_SOA_COLUMN(CosPAV0, cosPAV0, float);
DECLARE_SOA_COLUMN(PtV0, ptV0, float);
DECLARE_SOA_COLUMN(RadiusV0, radiusV0, float);
DECLARE_SOA_COLUMN(GammaPsiPair, gammaPsiPair, float);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
DECLARE_SOA_COLUMN(TrackOcc, trackOcc, float);
DECLARE_SOA_COLUMN(Ft0Occ, ft0Occ, float);
DECLARE_SOA_COLUMN(HadronicRate, hadronicRate, float);
} // namespace tpcskims
DECLARE_SOA_TABLE(SkimmedTPCV0Tree, "AOD", "TPCSKIMV0TREE",
                  o2::aod::track::TPCSignal,
                  tpcskims::InvDeDxExpTPC,
                  o2::aod::track::TPCInnerParam,
                  o2::aod::track::Tgl,
                  o2::aod::track::Signed1Pt,
                  o2::aod::track::Eta,
                  o2::aod::track::Phi,
                  o2::aod::track::Y,
                  tpcskims::Mass,
                  tpcskims::BetaGamma,
                  tpcskims::NormMultTPC,
                  tpcskims::NormNClustersTPC,
                  tpcskims::PidIndex,
                  tpcskims::NSigTPC,
                  tpcskims::NSigTOF,
                  tpcskims::AlphaV0,
                  tpcskims::QtV0,
                  tpcskims::CosPAV0,
                  tpcskims::PtV0,
                  tpcskims::RadiusV0,
                  tpcskims::GammaPsiPair,
                  tpcskims::RunNumber,
                  tpcskims::TrackOcc,
                  tpcskims::Ft0Occ,
                  tpcskims::HadronicRate);

DECLARE_SOA_TABLE(SkimmedTPCTOFTree, "AOD", "TPCTOFSKIMTREE",
                  o2::aod::track::TPCSignal,
                  tpcskims::InvDeDxExpTPC,
                  o2::aod::track::TPCInnerParam,
                  o2::aod::track::Tgl,
                  o2::aod::track::Signed1Pt,
                  o2::aod::track::Eta,
                  o2::aod::track::Phi,
                  o2::aod::track::Y,
                  tpcskims::Mass,
                  tpcskims::BetaGamma,
                  tpcskims::NormMultTPC,
                  tpcskims::NormNClustersTPC,
                  tpcskims::PidIndex,
                  tpcskims::NSigTPC,
                  tpcskims::NSigTOF,
                  tpcskims::RunNumber,
                  tpcskims::TrackOcc,
                  tpcskims::Ft0Occ,
                  tpcskims::HadronicRate);
} // namespace o2::aod
#endif // DPG_TASKS_TPC_TPCSKIMSTABLECREATOR_H_
