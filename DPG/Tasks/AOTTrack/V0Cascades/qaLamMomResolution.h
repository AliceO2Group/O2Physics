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

/// \file qaLamMomResolution.h

#ifndef DPG_TASKS_AOTTRACK_V0CASCADES_QALAMMOMRESOLUTION_H_
#define DPG_TASKS_AOTTRACK_V0CASCADES_QALAMMOMRESOLUTION_H_

#include <Framework/AnalysisDataModel.h>

namespace o2::aod
{
namespace qaLamMomResolution
{
DECLARE_SOA_COLUMN(CollID, collid, int);
// Lambda
DECLARE_SOA_COLUMN(MLambda, mLambda, float);
DECLARE_SOA_COLUMN(RLambda, rLambda, float);
DECLARE_SOA_COLUMN(PtLambda, ptLambda, float);

// charge
DECLARE_SOA_COLUMN(ChargeProton, chargeProton, int);
DECLARE_SOA_COLUMN(ChargePion, chargePion, int);
// eta
DECLARE_SOA_COLUMN(EtaProton, etaProton, float);
DECLARE_SOA_COLUMN(EtaPion, etaPion, float);
// phi
DECLARE_SOA_COLUMN(PhiProton, phiProton, float);
DECLARE_SOA_COLUMN(PhiPion, phiPion, float);
// nTPCclusters
DECLARE_SOA_COLUMN(NTPCClusProton, nTPCclusProton, int);
DECLARE_SOA_COLUMN(NTPCClusPion, nTPCclusPion, int);
// momentum
DECLARE_SOA_COLUMN(PxProton, pxProton, float);
DECLARE_SOA_COLUMN(PyProton, pyProton, float);
DECLARE_SOA_COLUMN(PzProton, pzProton, float);
DECLARE_SOA_COLUMN(PxPion, pxPion, float);
DECLARE_SOA_COLUMN(PyPion, pyPion, float);
DECLARE_SOA_COLUMN(PzPion, pzPion, float);
DECLARE_SOA_COLUMN(PxProtonIU, pxProtonIU, float);
DECLARE_SOA_COLUMN(PyProtonIU, pyProtonIU, float);
DECLARE_SOA_COLUMN(PzProtonIU, pzProtonIU, float);
DECLARE_SOA_COLUMN(PxPionIU, pxPionIU, float);
DECLARE_SOA_COLUMN(PyPionIU, pyPionIU, float);
DECLARE_SOA_COLUMN(PzPionIU, pzPionIU, float);
// momentum uncertainties
DECLARE_SOA_COLUMN(PxProtonErr, pxProtonErr, float);
DECLARE_SOA_COLUMN(PyProtonErr, pyProtonErr, float);
DECLARE_SOA_COLUMN(PzProtonErr, pzProtonErr, float);
DECLARE_SOA_COLUMN(PxPionErr, pxPionErr, float);
DECLARE_SOA_COLUMN(PyPionErr, pyPionErr, float);
DECLARE_SOA_COLUMN(PzPionErr, pzPionErr, float);
DECLARE_SOA_COLUMN(PxProtonIUErr, pxProtonIUErr, float);
DECLARE_SOA_COLUMN(PyProtonIUErr, pyProtonIUErr, float);
DECLARE_SOA_COLUMN(PzProtonIUErr, pzProtonIUErr, float);
DECLARE_SOA_COLUMN(PxPionIUErr, pxPionIUErr, float);
DECLARE_SOA_COLUMN(PyPionIUErr, pyPionIUErr, float);
DECLARE_SOA_COLUMN(PzPionIUErr, pzPionIUErr, float);
DECLARE_SOA_COLUMN(Sigma1PtProtonIU, sigma1ptProtonIU, float);
DECLARE_SOA_COLUMN(Sigma1PtPionIU, sigma1ptPionIU, float);
// IU position
DECLARE_SOA_COLUMN(XProtonIU, xProtonIU, float);
DECLARE_SOA_COLUMN(YProtonIU, yProtonIU, float);
DECLARE_SOA_COLUMN(ZProtonIU, zProtonIU, float);
DECLARE_SOA_COLUMN(XPionIU, xPionIU, float);
DECLARE_SOA_COLUMN(YPionIU, yPionIU, float);
DECLARE_SOA_COLUMN(ZPionIU, zPionIU, float);
// IU position uncertainties
DECLARE_SOA_COLUMN(XProtonIUErr, xProtonIUErr, float);
DECLARE_SOA_COLUMN(YProtonIUErr, yProtonIUErr, float);
DECLARE_SOA_COLUMN(ZProtonIUErr, zProtonIUErr, float);
DECLARE_SOA_COLUMN(XPionIUErr, xPionIUErr, float);
DECLARE_SOA_COLUMN(YPionIUErr, yPionIUErr, float);
DECLARE_SOA_COLUMN(ZPionIUErr, zPionIUErr, float);
// DCA
DECLARE_SOA_COLUMN(DCAxyProton, dcaxyProton, float);
DECLARE_SOA_COLUMN(DCAxyPion, dcaxyPion, float);
DECLARE_SOA_COLUMN(DCAzProton, dcazProton, float);
DECLARE_SOA_COLUMN(DCAzPion, dcazPion, float);
// DCA uncertainties
DECLARE_SOA_COLUMN(DCAxyProtonErr, dcaxyProtonErr, float);
DECLARE_SOA_COLUMN(DCAxyPionErr, dcaxyPionErr, float);
DECLARE_SOA_COLUMN(DCAzProtonErr, dcazProtonErr, float);
DECLARE_SOA_COLUMN(DCAzPionErr, dcazPionErr, float);
// MC
DECLARE_SOA_COLUMN(PxProtonMC, pxProtonMC, float);
DECLARE_SOA_COLUMN(PyProtonMC, pyProtonMC, float);
DECLARE_SOA_COLUMN(PzProtonMC, pzProtonMC, float);
DECLARE_SOA_COLUMN(PxPionMC, pxPionMC, float);
DECLARE_SOA_COLUMN(PyPionMC, pyPionMC, float);
DECLARE_SOA_COLUMN(PzPionMC, pzPionMC, float);
} // namespace qaLamMomResolution

DECLARE_SOA_TABLE(LamDaughters, "AOD", "LAMDAUGHTERS",
                  qaLamMomResolution::CollID,
                  qaLamMomResolution::MLambda,
                  qaLamMomResolution::RLambda,
                  qaLamMomResolution::PtLambda,
                  qaLamMomResolution::ChargeProton,
                  qaLamMomResolution::ChargePion,
                  qaLamMomResolution::EtaProton,
                  qaLamMomResolution::EtaPion,
                  qaLamMomResolution::PhiProton,
                  qaLamMomResolution::PhiPion,
                  qaLamMomResolution::NTPCClusProton,
                  qaLamMomResolution::NTPCClusPion,
                  qaLamMomResolution::PxProton,
                  qaLamMomResolution::PyProton,
                  qaLamMomResolution::PzProton,
                  qaLamMomResolution::PxProtonErr,
                  qaLamMomResolution::PyProtonErr,
                  qaLamMomResolution::PzProtonErr,
                  qaLamMomResolution::PxPion,
                  qaLamMomResolution::PyPion,
                  qaLamMomResolution::PzPion,
                  qaLamMomResolution::PxPionErr,
                  qaLamMomResolution::PyPionErr,
                  qaLamMomResolution::PzPionErr,
                  qaLamMomResolution::PxProtonIU,
                  qaLamMomResolution::PyProtonIU,
                  qaLamMomResolution::PzProtonIU,
                  qaLamMomResolution::PxProtonIUErr,
                  qaLamMomResolution::PyProtonIUErr,
                  qaLamMomResolution::PzProtonIUErr,
                  qaLamMomResolution::PxPionIU,
                  qaLamMomResolution::PyPionIU,
                  qaLamMomResolution::PzPionIU,
                  qaLamMomResolution::PxPionIUErr,
                  qaLamMomResolution::PyPionIUErr,
                  qaLamMomResolution::PzPionIUErr,
                  qaLamMomResolution::PxProtonMC,
                  qaLamMomResolution::PyProtonMC,
                  qaLamMomResolution::PzProtonMC,
                  qaLamMomResolution::PxPionMC,
                  qaLamMomResolution::PyPionMC,
                  qaLamMomResolution::PzPionMC,
                  qaLamMomResolution::Sigma1PtProtonIU,
                  qaLamMomResolution::Sigma1PtPionIU,
                  qaLamMomResolution::XProtonIU,
                  qaLamMomResolution::YProtonIU,
                  qaLamMomResolution::ZProtonIU,
                  qaLamMomResolution::XProtonIUErr,
                  qaLamMomResolution::YProtonIUErr,
                  qaLamMomResolution::ZProtonIUErr,
                  qaLamMomResolution::XPionIU,
                  qaLamMomResolution::YPionIU,
                  qaLamMomResolution::ZPionIU,
                  qaLamMomResolution::XPionIUErr,
                  qaLamMomResolution::YPionIUErr,
                  qaLamMomResolution::ZPionIUErr,
                  qaLamMomResolution::DCAxyProton,
                  qaLamMomResolution::DCAzProton,
                  qaLamMomResolution::DCAxyProtonErr,
                  qaLamMomResolution::DCAzProtonErr,
                  qaLamMomResolution::DCAxyPion,
                  qaLamMomResolution::DCAzPion,
                  qaLamMomResolution::DCAxyPionErr,
                  qaLamMomResolution::DCAzPionErr);
} // namespace o2::aod

#endif // DPG_TASKS_AOTTRACK_V0CASCADES_QALAMMOMRESOLUTION_H_
