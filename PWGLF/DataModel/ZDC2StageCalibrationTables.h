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

/// \author Sourav Kundu <sourav.kundu@cern.ch>

#ifndef PWGLF_DATAMODEL_ZDC2STAGECALIBRATIONTABLES_H_
#define PWGLF_DATAMODEL_ZDC2STAGECALIBRATIONTABLES_H_

#include <Framework/AnalysisDataModel.h>

#include <cstdint>

namespace o2::aod
{
namespace zdc2stagecalib
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_COLUMN(TriggerSP, triggerSP, bool);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int32_t);
DECLARE_SOA_COLUMN(GlobalBC, globalBC, uint64_t);
DECLARE_SOA_COLUMN(Centrality, centrality, float);
DECLARE_SOA_COLUMN(Vx, vx, float);
DECLARE_SOA_COLUMN(Vy, vy, float);
DECLARE_SOA_COLUMN(Vz, vz, float);
DECLARE_SOA_COLUMN(QxZDCA, qxZDCA, float);
DECLARE_SOA_COLUMN(QxZDCC, qxZDCC, float);
DECLARE_SOA_COLUMN(QyZDCA, qyZDCA, float);
DECLARE_SOA_COLUMN(QyZDCC, qyZDCC, float);
} // namespace zdc2stagecalib

DECLARE_SOA_TABLE(ZDC2StageCalibs, "AOD", "ZDC2STGCAL",
                  zdc2stagecalib::CollisionId,
                  zdc2stagecalib::TriggerSP,
                  zdc2stagecalib::RunNumber,
                  zdc2stagecalib::GlobalBC,
                  zdc2stagecalib::Centrality,
                  zdc2stagecalib::Vx,
                  zdc2stagecalib::Vy,
                  zdc2stagecalib::Vz,
                  zdc2stagecalib::QxZDCA,
                  zdc2stagecalib::QxZDCC,
                  zdc2stagecalib::QyZDCA,
                  zdc2stagecalib::QyZDCC);
} // namespace o2::aod

#endif // PWGLF_DATAMODEL_ZDC2STAGECALIBRATIONTABLES_H_
