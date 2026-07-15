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

/// \file treeCreatorPidTpcQa.h
/// \author Ana Marin <ana.marin@cern.ch>
/// \author Oleksii Lubynets <oleksii.lubynets@cern.ch>
/// \brief  Creates trees with PID QA variables along with variables used for NN training

#ifndef DPG_TASKS_AOTTRACK_PID_TPC_TREECREATORPIDTPCQA_H_
#define DPG_TASKS_AOTTRACK_PID_TPC_TREECREATORPIDTPCQA_H_

#include "DPG/Tasks/TPC/tpcSkimsTableCreator.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>

namespace o2::aod
{
namespace dpg_tpcpidqa
{
DECLARE_SOA_COLUMN(NSigmaTpc, nSigmaTpc, float);
DECLARE_SOA_COLUMN(DedxExpected, dedxExpected, float);
DECLARE_SOA_COLUMN(DedxDiff, dedxDiff, float);
DECLARE_SOA_COLUMN(ExpSigma, expSigma, float);
DECLARE_SOA_COLUMN(NSigmaTof, nSigmaTof, float);
} // namespace dpg_tpcpidqa

DECLARE_SOA_TABLE(QaPidTpc, "AOD", "QAPIDTPC",
                  tpcskims::PidIndex,
                  tpcskims::Ft0Occ,
                  tpcskims::HadronicRate,
                  tpcskims::NormMultTPC,
                  tpcskims::NormNClustersTPC,
                  tpcskims::NormNClustersTPCPID,
                  o2::aod::track::Phi,
                  o2::aod::track::Tgl,
                  o2::aod::track::TPCInnerParam,
                  o2::aod::track::Y,
                  o2::aod::track::P,
                  o2::aod::track::Signed1Pt,
                  dpg_tpcpidqa::NSigmaTpc,
                  dpg_tpcpidqa::DedxExpected,
                  dpg_tpcpidqa::DedxDiff,
                  dpg_tpcpidqa::ExpSigma,
                  dpg_tpcpidqa::NSigmaTof)
} // namespace o2::aod

#endif // DPG_TASKS_AOTTRACK_PID_TPC_TREECREATORPIDTPCQA_H_
