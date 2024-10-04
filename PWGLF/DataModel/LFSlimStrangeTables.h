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

#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"

#ifndef PWGLF_DATAMODEL_LFSLIMSTRANGETABLES_H_
#define PWGLF_DATAMODEL_LFSLIMSTRANGETABLES_H_

namespace o2::aod
{

namespace SlimLambdaTables
{
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Mass, mass, float);
DECLARE_SOA_COLUMN(Ct, ct, float);
DECLARE_SOA_COLUMN(Len, len, float);
DECLARE_SOA_COLUMN(Radius, radius, float);
DECLARE_SOA_COLUMN(DcaV0PV, dcaV0Pv, float);
DECLARE_SOA_COLUMN(DcaPosPV, dcaPosPv, float);
DECLARE_SOA_COLUMN(DcaNegPV, dcaNegPv, float);
DECLARE_SOA_COLUMN(DcaV0Tracks, dcaV0tracks, float);
DECLARE_SOA_COLUMN(CosPA, cosPa, double);
DECLARE_SOA_COLUMN(AlphaAP, alphaAP, double);
DECLARE_SOA_COLUMN(QtAP, qtAP, double);
DECLARE_SOA_COLUMN(TpcNsigmaPos, tpcNsigmaPos, float);
DECLARE_SOA_COLUMN(TpcNsigmaNeg, tpcNsigmaNeg, float);
DECLARE_SOA_COLUMN(PxPos, pxPos, float);
DECLARE_SOA_COLUMN(PyPos, pyPos, float);
DECLARE_SOA_COLUMN(PzPos, pzPos, float);
DECLARE_SOA_COLUMN(PxNeg, pxNeg, float);
DECLARE_SOA_COLUMN(PyNeg, pyNeg, float);
DECLARE_SOA_COLUMN(PzNeg, pzNeg, float);
DECLARE_SOA_COLUMN(PxPosMC, pxPosMC, float);
DECLARE_SOA_COLUMN(PyPosMC, pyPosMC, float);
DECLARE_SOA_COLUMN(PzPosMC, pzPosMC, float);
DECLARE_SOA_COLUMN(PxNegMC, pxNegMC, float);
DECLARE_SOA_COLUMN(PyNegMC, pyNegMC, float);
DECLARE_SOA_COLUMN(PzNegMC, pzNegMC, float);
DECLARE_SOA_COLUMN(GenPt, gentPt, float);
DECLARE_SOA_COLUMN(GenEta, genEta, float);
DECLARE_SOA_COLUMN(GenCt, genCt, float);
DECLARE_SOA_COLUMN(GenLen, genLen, float);
DECLARE_SOA_COLUMN(PDGCodeDauPos, pdgCodeDauPos, int);
DECLARE_SOA_COLUMN(PDGCodeMotherDauPos, pdgCodeMotherDauPos, int);
DECLARE_SOA_COLUMN(PDGCodeDauNeg, pdgCodeDauNeg, int);
DECLARE_SOA_COLUMN(PDGCodeMotherDauNeg, pdgCodeMotherDauNeg, int);
DECLARE_SOA_COLUMN(PDGCode, pdgCode, int);
DECLARE_SOA_COLUMN(PDGCodeMother, pdgCodeMother, int);
DECLARE_SOA_COLUMN(IsReco, isReco, bool);
DECLARE_SOA_COLUMN(IsFD, isFD, uint8_t);
DECLARE_SOA_COLUMN(PDGMatchMotherSecondMother, pdgMatchMotherSecondMother, int);
} // namespace SlimLambdaTables

DECLARE_SOA_TABLE(LambdaTableML, "AOD", "LAMBDATABLEML",
                  SlimLambdaTables::Pt,
                  SlimLambdaTables::Eta,
                  SlimLambdaTables::Mass,
                  SlimLambdaTables::Ct,
                  SlimLambdaTables::Radius,
                  SlimLambdaTables::DcaV0PV,
                  SlimLambdaTables::DcaPosPV,
                  SlimLambdaTables::DcaNegPV,
                  SlimLambdaTables::DcaV0Tracks,
                  SlimLambdaTables::CosPA,
                  SlimLambdaTables::AlphaAP,
                  SlimLambdaTables::QtAP,
                  SlimLambdaTables::TpcNsigmaPos,
                  SlimLambdaTables::TpcNsigmaNeg,
                  SlimLambdaTables::IsFD);

DECLARE_SOA_TABLE(McLambdaTableML, "AOD", "MCLAMBDATABLEML",
                  SlimLambdaTables::Pt,
                  SlimLambdaTables::Eta,
                  SlimLambdaTables::Mass,
                  SlimLambdaTables::Ct,
                  SlimLambdaTables::Radius,
                  SlimLambdaTables::DcaV0PV,
                  SlimLambdaTables::DcaPosPV,
                  SlimLambdaTables::DcaNegPV,
                  SlimLambdaTables::DcaV0Tracks,
                  SlimLambdaTables::CosPA,
                  SlimLambdaTables::AlphaAP,
                  SlimLambdaTables::QtAP,
                  SlimLambdaTables::TpcNsigmaPos,
                  SlimLambdaTables::TpcNsigmaNeg,
                  SlimLambdaTables::GenPt,
                  SlimLambdaTables::GenEta,
                  SlimLambdaTables::GenCt,
                  SlimLambdaTables::PDGCodeDauPos,
                  SlimLambdaTables::PDGCodeMotherDauPos,
                  SlimLambdaTables::PDGCodeDauNeg,
                  SlimLambdaTables::PDGCodeMotherDauNeg,
                  SlimLambdaTables::PDGCode,
                  SlimLambdaTables::PDGCodeMother,
                  SlimLambdaTables::IsReco,
                  SlimLambdaTables::PDGMatchMotherSecondMother);

DECLARE_SOA_TABLE(V0TableAP, "AOD", "V0TABLEAP",
                  SlimLambdaTables::Eta,
                  SlimLambdaTables::Len,
                  SlimLambdaTables::PxPos,
                  SlimLambdaTables::PyPos,
                  SlimLambdaTables::PzPos,
                  SlimLambdaTables::PxNeg,
                  SlimLambdaTables::PyNeg,
                  SlimLambdaTables::PzNeg,
                  SlimLambdaTables::Radius,
                  SlimLambdaTables::DcaV0PV,
                  SlimLambdaTables::DcaPosPV,
                  SlimLambdaTables::DcaNegPV,
                  SlimLambdaTables::DcaV0Tracks,
                  SlimLambdaTables::CosPA);

DECLARE_SOA_TABLE(McV0TableAP, "AOD", "MCV0TABLEAP",
                  SlimLambdaTables::Eta,
                  SlimLambdaTables::Len,
                  SlimLambdaTables::PxPos,
                  SlimLambdaTables::PyPos,
                  SlimLambdaTables::PzPos,
                  SlimLambdaTables::PxNeg,
                  SlimLambdaTables::PyNeg,
                  SlimLambdaTables::PzNeg,
                  SlimLambdaTables::PxPosMC,
                  SlimLambdaTables::PyPosMC,
                  SlimLambdaTables::PzPosMC,
                  SlimLambdaTables::PxNegMC,
                  SlimLambdaTables::PyNegMC,
                  SlimLambdaTables::PzNegMC,
                  SlimLambdaTables::Radius,
                  SlimLambdaTables::DcaV0PV,
                  SlimLambdaTables::DcaPosPV,
                  SlimLambdaTables::DcaNegPV,
                  SlimLambdaTables::DcaV0Tracks,
                  SlimLambdaTables::CosPA,
                  SlimLambdaTables::GenEta,
                  SlimLambdaTables::GenLen,
                  SlimLambdaTables::PDGCode,
                  SlimLambdaTables::IsReco);

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFSLIMSTRANGETABLES_H_
