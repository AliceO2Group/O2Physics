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

#ifndef PWGLF_DATAMODEL_LFDOUBLEOMEGATABLES_H_
#define PWGLF_DATAMODEL_LFDOUBLEOMEGATABLES_H_

#include <Framework/AnalysisDataModel.h>

#include <cstdint>

namespace o2::aod
{

namespace DoubleOmegaTables
{
DECLARE_SOA_COLUMN(PtCasc1, ptCasc1, float); // signed pt of the cascade
DECLARE_SOA_COLUMN(EtaCasc1, etaCasc1, float);
DECLARE_SOA_COLUMN(PhiCasc1, phiCasc1, float);
DECLARE_SOA_COLUMN(CascDecLength1, cascDecLength1, float);
DECLARE_SOA_COLUMN(OmegaMassCasc1, omegaMassCasc1, float);
DECLARE_SOA_COLUMN(XiMassCasc1, xiMassCasc1, float);
DECLARE_SOA_COLUMN(CosPACasc1, cosPACasc1, float);
DECLARE_SOA_COLUMN(DcaBachPVCasc1, dcaBachPVCasc1, float);
DECLARE_SOA_COLUMN(DcaV0BachCasc1, dcaV0BachCasc1, float);
DECLARE_SOA_COLUMN(NSigmaKBach1, nSigmaKBach1, float);

DECLARE_SOA_COLUMN(PtLambda, ptLambda, float);
DECLARE_SOA_COLUMN(EtaLambda, etaLambda, float);
DECLARE_SOA_COLUMN(PhiLambda, phiLambda, float);
DECLARE_SOA_COLUMN(LambdaDecLength, lambdaDecLength, float);
DECLARE_SOA_COLUMN(LambdaMass, lambdaMass, float);
DECLARE_SOA_COLUMN(CosPALambda, cosPALambda, float);
DECLARE_SOA_COLUMN(DcaPosPVLambda, dcaPosPVLambda, float);
DECLARE_SOA_COLUMN(DcaNegPVLambda, dcaNegPVLambda, float);
DECLARE_SOA_COLUMN(DcaLambdaDaughters, dcaLambdaDaughters, float);
DECLARE_SOA_COLUMN(NSigmaPrLambda, nSigmaPrLambda, float);
DECLARE_SOA_COLUMN(NSigmaPiLambda, nSigmaPiLambda, float);

DECLARE_SOA_COLUMN(PtKaon, ptKaon, float);
DECLARE_SOA_COLUMN(EtaKaon, etaKaon, float);
DECLARE_SOA_COLUMN(PhiKaon, phiKaon, float);
DECLARE_SOA_COLUMN(ChargeKaon, chargeKaon, int8_t);
DECLARE_SOA_COLUMN(DcaXYKaon, dcaXYKaon, float);
DECLARE_SOA_COLUMN(DcaZKaon, dcaZKaon, float);
DECLARE_SOA_COLUMN(NSigmaKKaon, nSigmaKKaon, float);

DECLARE_SOA_COLUMN(DoubleOmegaMass, doubleOmegaMass, float);

DECLARE_SOA_COLUMN(GenPt, genPt, float);
DECLARE_SOA_COLUMN(GenEta, genEta, float);
DECLARE_SOA_COLUMN(GenPhi, genPhi, float);
DECLARE_SOA_COLUMN(GenDecayLength, genDecayLength, float);
DECLARE_SOA_COLUMN(PdgDoubleOmega, pdgDoubleOmega, int);
DECLARE_SOA_COLUMN(IsReco, isReco, bool);
} // namespace DoubleOmegaTables

#define DOUBLE_OMEGA_RECO_COLUMNS          \
  DoubleOmegaTables::PtCasc1,              \
    DoubleOmegaTables::EtaCasc1,           \
    DoubleOmegaTables::PhiCasc1,           \
    DoubleOmegaTables::CascDecLength1,     \
    DoubleOmegaTables::OmegaMassCasc1,     \
    DoubleOmegaTables::XiMassCasc1,        \
    DoubleOmegaTables::CosPACasc1,         \
    DoubleOmegaTables::DcaBachPVCasc1,     \
    DoubleOmegaTables::DcaV0BachCasc1,     \
    DoubleOmegaTables::NSigmaKBach1,       \
    DoubleOmegaTables::PtLambda,           \
    DoubleOmegaTables::EtaLambda,          \
    DoubleOmegaTables::PhiLambda,          \
    DoubleOmegaTables::LambdaDecLength,    \
    DoubleOmegaTables::LambdaMass,         \
    DoubleOmegaTables::CosPALambda,        \
    DoubleOmegaTables::DcaPosPVLambda,     \
    DoubleOmegaTables::DcaNegPVLambda,     \
    DoubleOmegaTables::DcaLambdaDaughters, \
    DoubleOmegaTables::NSigmaPrLambda,     \
    DoubleOmegaTables::NSigmaPiLambda,     \
    DoubleOmegaTables::PtKaon,             \
    DoubleOmegaTables::EtaKaon,            \
    DoubleOmegaTables::PhiKaon,            \
    DoubleOmegaTables::ChargeKaon,         \
    DoubleOmegaTables::DcaXYKaon,          \
    DoubleOmegaTables::DcaZKaon,           \
    DoubleOmegaTables::NSigmaKKaon,        \
    DoubleOmegaTables::DoubleOmegaMass

DECLARE_SOA_TABLE(DoubleOmegaTable, "AOD", "DOUBLEOMEGATABLE",
                  DOUBLE_OMEGA_RECO_COLUMNS);

DECLARE_SOA_TABLE(DoubleOmegaTableMC, "AOD", "DBLOMEGAMCTABLE",
                  DOUBLE_OMEGA_RECO_COLUMNS,
                  DoubleOmegaTables::GenPt,
                  DoubleOmegaTables::GenEta,
                  DoubleOmegaTables::GenPhi,
                  DoubleOmegaTables::GenDecayLength,
                  DoubleOmegaTables::PdgDoubleOmega,
                  DoubleOmegaTables::IsReco);

#undef DOUBLE_OMEGA_RECO_COLUMNS

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFDOUBLEOMEGATABLES_H_
