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

#ifndef DPG_TASKS_AOTTRACK_PID_HMPID_TABLEHMPID_H_
#define DPG_TASKS_AOTTRACK_PID_HMPID_TABLEHMPID_H_

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{

inline constexpr int kDimPhotonsCharge = 10;

namespace variables_table
{
DECLARE_SOA_COLUMN(ChAngle, chAngle, float);
DECLARE_SOA_COLUMN(MomentumHmpid, momentumHmpid, float);
DECLARE_SOA_COLUMN(MomentumTrack, momentumTrack, float);
DECLARE_SOA_COLUMN(XTrack, xTrack, float);
DECLARE_SOA_COLUMN(YTrack, yTrack, float);
DECLARE_SOA_COLUMN(XMip, xMip, float);
DECLARE_SOA_COLUMN(YMip, yMip, float);
DECLARE_SOA_COLUMN(NPhotons, nPhotons, float);
DECLARE_SOA_COLUMN(ChargeMip, chargeMip, float);
DECLARE_SOA_COLUMN(ClusterSize, clusterSize, float);
DECLARE_SOA_COLUMN(Chamber, chamber, float);
DECLARE_SOA_COLUMN(PhotonsCharge, photonsCharge, float[kDimPhotonsCharge]);

DECLARE_SOA_COLUMN(EtaTrack, etaTrack, float);
DECLARE_SOA_COLUMN(PhiTrack, phiTrack, float);
DECLARE_SOA_COLUMN(Px, px, float);
DECLARE_SOA_COLUMN(Py, py, float);
DECLARE_SOA_COLUMN(Pz, pz, float);

DECLARE_SOA_COLUMN(ItsNCluster, itsNCluster, float);
DECLARE_SOA_COLUMN(TpcNCluster, tpcNCluster, float);
DECLARE_SOA_COLUMN(TpcNClsCrossedRows, tpcNClsCrossedRows, float);
DECLARE_SOA_COLUMN(TpcChi2, tpcChi2, float);
DECLARE_SOA_COLUMN(ItsChi2, itsChi2, float);

DECLARE_SOA_COLUMN(DcaXY, dcaXY, float);
DECLARE_SOA_COLUMN(DcaZ, dcaZ, float);

DECLARE_SOA_COLUMN(TpcNSigmaPi, tpcNSigmaPi, float);
DECLARE_SOA_COLUMN(TofNSigmaPi, tofNSigmaPi, float);
DECLARE_SOA_COLUMN(TpcNSigmaKa, tpcNSigmaKa, float);
DECLARE_SOA_COLUMN(TofNSigmaKa, tofNSigmaKa, float);
DECLARE_SOA_COLUMN(TpcNSigmaPr, tpcNSigmaPr, float);
DECLARE_SOA_COLUMN(TofNSigmaPr, tofNSigmaPr, float);
DECLARE_SOA_COLUMN(TpcNSigmaDe, tpcNSigmaDe, float);
DECLARE_SOA_COLUMN(TofNSigmaDe, tofNSigmaDe, float);
DECLARE_SOA_COLUMN(Centrality, centrality, float);

} // namespace variables_table

DECLARE_SOA_TABLE(HmpidAnalysis, "AOD", "HMPIDANALYSIS",
                  variables_table::ChAngle,
                  variables_table::MomentumHmpid,
                  variables_table::MomentumTrack,
                  variables_table::XTrack,
                  variables_table::YTrack,
                  variables_table::XMip,
                  variables_table::YMip,
                  variables_table::NPhotons,
                  variables_table::ChargeMip,
                  variables_table::ClusterSize,
                  variables_table::Chamber,
                  variables_table::PhotonsCharge,
                  variables_table::EtaTrack,
                  variables_table::PhiTrack,
                  variables_table::Px,
                  variables_table::Py,
                  variables_table::Pz,
                  variables_table::ItsNCluster,
                  variables_table::TpcNCluster,
                  variables_table::TpcNClsCrossedRows,
                  variables_table::TpcChi2,
                  variables_table::ItsChi2,
                  variables_table::DcaXY,
                  variables_table::DcaZ,
                  variables_table::TpcNSigmaPi,
                  variables_table::TofNSigmaPi,
                  variables_table::TpcNSigmaKa,
                  variables_table::TofNSigmaKa,
                  variables_table::TpcNSigmaPr,
                  variables_table::TofNSigmaPr,
                  variables_table::TpcNSigmaDe,
                  variables_table::TofNSigmaDe,
                  variables_table::Centrality);

} // namespace o2::aod

#endif // DPG_TASKS_AOTTRACK_PID_HMPID_TABLEHMPID_H_
