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
DECLARE_SOA_COLUMN(CentFT0C, centFT0C, float);
DECLARE_SOA_COLUMN(IsMatter, isMatter, bool);
DECLARE_SOA_COLUMN(Mass, mass, float);
DECLARE_SOA_COLUMN(Ct, ct, float);
DECLARE_SOA_COLUMN(Radius, radius, float);
DECLARE_SOA_COLUMN(DcaV0PV, dcaV0Pv, float);
DECLARE_SOA_COLUMN(DcaPiPV, dcaPiPv, float);
DECLARE_SOA_COLUMN(DcaPrPV, dcaPrPv, float);
DECLARE_SOA_COLUMN(DcaV0Tracks, dcaV0tracks, float);
DECLARE_SOA_COLUMN(CosPA, cosPa, double);
DECLARE_SOA_COLUMN(TpcNsigmaPi, tpcNsigmaPi, float);
DECLARE_SOA_COLUMN(TpcNsigmaPr, tpcNsigmaPr, float);
DECLARE_SOA_COLUMN(GenPt, gentPt, float);
DECLARE_SOA_COLUMN(GenEta, genEta, float);
DECLARE_SOA_COLUMN(GenCt, genCt, float);
DECLARE_SOA_COLUMN(PDGCode, pdgCode, int);
DECLARE_SOA_COLUMN(IsReco, isReco, bool);
} // namespace SlimLambdaTables

DECLARE_SOA_TABLE(LambdaTableML, "AOD", "LAMBDATABLEML",
                  SlimLambdaTables::Pt,
                  SlimLambdaTables::Eta,
                  SlimLambdaTables::CentFT0C,
                  SlimLambdaTables::IsMatter,
                  SlimLambdaTables::Mass,
                  SlimLambdaTables::Ct,
                  SlimLambdaTables::Radius,
                  SlimLambdaTables::DcaV0PV,
                  SlimLambdaTables::DcaPiPV,
                  SlimLambdaTables::DcaPrPV,
                  SlimLambdaTables::DcaV0Tracks,
                  SlimLambdaTables::CosPA,
                  SlimLambdaTables::TpcNsigmaPi,
                  SlimLambdaTables::TpcNsigmaPr);

DECLARE_SOA_TABLE(McLambdaTableML, "AOD", "MCLAMBDATABLEML",
                  SlimLambdaTables::Pt,
                  SlimLambdaTables::Eta,
                  SlimLambdaTables::CentFT0C,
                  SlimLambdaTables::IsMatter,
                  SlimLambdaTables::Mass,
                  SlimLambdaTables::Ct,
                  SlimLambdaTables::Radius,
                  SlimLambdaTables::DcaV0PV,
                  SlimLambdaTables::DcaPiPV,
                  SlimLambdaTables::DcaPrPV,
                  SlimLambdaTables::DcaV0Tracks,
                  SlimLambdaTables::CosPA,
                  SlimLambdaTables::TpcNsigmaPi,
                  SlimLambdaTables::TpcNsigmaPr,
                  SlimLambdaTables::GenPt,
                  SlimLambdaTables::GenEta,
                  SlimLambdaTables::GenCt,
                  SlimLambdaTables::PDGCode,
                  SlimLambdaTables::IsReco);

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFSLIMSTRANGETABLES_H_
