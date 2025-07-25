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

#ifndef PWGLF_DATAMODEL_LFDOUBLECASCTABLES_H_
#define PWGLF_DATAMODEL_LFDOUBLECASCTABLES_H_

namespace o2::aod
{

namespace DoubleCascTables
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

DECLARE_SOA_COLUMN(PtCasc2, ptCasc2, float);
DECLARE_SOA_COLUMN(EtaCasc2, etaCasc2, float);
DECLARE_SOA_COLUMN(PhiCasc2, phiCasc2, float);
DECLARE_SOA_COLUMN(CascDecLength2, cascDecLength2, float);
DECLARE_SOA_COLUMN(OmegaMassCasc2, omegaMassCasc2, float);
DECLARE_SOA_COLUMN(XiMassCasc2, xiMassCasc2, float);
DECLARE_SOA_COLUMN(CosPACasc2, cosPACasc2, float);
DECLARE_SOA_COLUMN(DcaBachPVCasc2, dcaBachPVCasc2, float);
DECLARE_SOA_COLUMN(DcaV0BachCasc2, dcaV0BachCasc2, float);
DECLARE_SOA_COLUMN(NSigmaKBach2, nSigmaKBach2, float);

DECLARE_SOA_COLUMN(DoubleOmegaMass, doubleOmegaMass, float);
} // namespace DoubleCascTables

DECLARE_SOA_TABLE(DoubleCascTable, "AOD", "DOUBLECASCTABLE",
                  DoubleCascTables::PtCasc1,
                  DoubleCascTables::EtaCasc1,
                  DoubleCascTables::PhiCasc1,
                  DoubleCascTables::CascDecLength1,
                  DoubleCascTables::OmegaMassCasc1,
                  DoubleCascTables::XiMassCasc1,
                  DoubleCascTables::CosPACasc1,
                  DoubleCascTables::DcaBachPVCasc1,
                  DoubleCascTables::DcaV0BachCasc1,
                  DoubleCascTables::NSigmaKBach1,
                  DoubleCascTables::PtCasc2,
                  DoubleCascTables::EtaCasc2,
                  DoubleCascTables::PhiCasc2,
                  DoubleCascTables::CascDecLength2,
                  DoubleCascTables::OmegaMassCasc2,
                  DoubleCascTables::XiMassCasc2,
                  DoubleCascTables::CosPACasc2,
                  DoubleCascTables::DcaBachPVCasc2,
                  DoubleCascTables::DcaV0BachCasc2,
                  DoubleCascTables::NSigmaKBach2,
                  DoubleCascTables::DoubleOmegaMass);

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFDOUBLECASCTABLES_H_
