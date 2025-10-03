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

#ifndef O2_ANALYSIS_HMPIDTABLEPB_H
#define O2_ANALYSIS_HMPIDTABLEPB_H

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

namespace o2::aod
{

namespace variables_table // declaration of columns to create
{
DECLARE_SOA_COLUMN(ChAngle, chAngle, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(MomentumHMPID, momentumHMPID, float);
DECLARE_SOA_COLUMN(MomentumTrack, momentumTrack, float);
DECLARE_SOA_COLUMN(Xtrack, xtrack, float);
DECLARE_SOA_COLUMN(Ytrack, ytrack, float);
DECLARE_SOA_COLUMN(Xmip, xmip, float);
DECLARE_SOA_COLUMN(Ymip, ymip, float);
DECLARE_SOA_COLUMN(Nphotons, nphotons, float);
DECLARE_SOA_COLUMN(ChargeMIP, chargeMIP, float);
DECLARE_SOA_COLUMN(ClusterSize, clustersize, float);
DECLARE_SOA_COLUMN(Chamber, chamber, float);
DECLARE_SOA_COLUMN(Photons_charge, photons_charge, float[10]);

DECLARE_SOA_COLUMN(EtaTrack, etatrack, float);
DECLARE_SOA_COLUMN(PhiTrack, phitrack, float);
DECLARE_SOA_COLUMN(Px, px, float);
DECLARE_SOA_COLUMN(Py, py, float);
DECLARE_SOA_COLUMN(Pz, pz, float);

DECLARE_SOA_COLUMN(ITSNcluster, itsNcluster, float);
DECLARE_SOA_COLUMN(TPCNcluster, tpcNcluster, float);
DECLARE_SOA_COLUMN(TPCNClsCrossedRows, tpcNClsCrossedRows, float);
DECLARE_SOA_COLUMN(TPCchi2, tpcChi2, float);
DECLARE_SOA_COLUMN(ITSchi2, itsChi2, float);

DECLARE_SOA_COLUMN(DCAxy, dcaxy, float);
DECLARE_SOA_COLUMN(DCAz, dcaz, float);

DECLARE_SOA_COLUMN(TPCNSigmaPi, tpcNsigmaPi, float);
DECLARE_SOA_COLUMN(TOFNSigmaPi, tofNsigmaPi, float);
DECLARE_SOA_COLUMN(TPCNSigmaKa, tpcNsigmaKa, float);
DECLARE_SOA_COLUMN(TOFNSigmaKa, tofNsigmaKa, float);
DECLARE_SOA_COLUMN(TPCNSigmaPr, tpcNsigmaPr, float);
DECLARE_SOA_COLUMN(TOFNSigmaPr, tofNsigmaPr, float);
DECLARE_SOA_COLUMN(TPCNSigmaDe, tpcNsigmaDe, float);
DECLARE_SOA_COLUMN(TOFNSigmaDe, tofNsigmaDe, float);

DECLARE_SOA_COLUMN(Centrality, centrality, float);

} // namespace variables_table

DECLARE_SOA_TABLE(HMPID_analysisPb, "AOD", "HMPIDANALYSISPB",
                  variables_table::ChAngle, variables_table::Phi, variables_table::Eta, variables_table::MomentumHMPID,
                  variables_table::MomentumTrack, variables_table::Xtrack, variables_table::Ytrack, variables_table::Xmip,
                  variables_table::Ymip, variables_table::Nphotons, variables_table::ChargeMIP, variables_table::ClusterSize,
                  variables_table::Chamber, variables_table::Photons_charge, variables_table::EtaTrack, variables_table::PhiTrack,
                  variables_table::Px, variables_table::Py, variables_table::Pz,
                  variables_table::ITSNcluster, variables_table::TPCNcluster, variables_table::TPCNClsCrossedRows,
                  variables_table::TPCchi2, variables_table::ITSchi2, variables_table::DCAxy, variables_table::DCAz,
                  variables_table::TPCNSigmaPi, variables_table::TOFNSigmaPi, variables_table::TPCNSigmaKa, variables_table::TOFNSigmaKa,
                  variables_table::TPCNSigmaPr, variables_table::TOFNSigmaPr, variables_table::TPCNSigmaDe, variables_table::TOFNSigmaDe,
                  variables_table::Centrality);
} // namespace o2::aod

#endif // O2_ANALYSIS_HMPIDTABLEPB_H
