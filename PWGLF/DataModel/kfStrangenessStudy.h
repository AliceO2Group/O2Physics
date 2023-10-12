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

/// \file kfStrangenessStudy.h

#ifndef PWGLF_DATAMODEL_KFSTRANGENESSSTUDY_H_
#define PWGLF_DATAMODEL_KFSTRANGENESSSTUDY_H_

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/Core/trackUtilities.h"

namespace o2::aod
{
namespace kfStrangenessStudy
{
DECLARE_SOA_COLUMN(CollID, collid, int);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(PtKF, ptkf, float);
DECLARE_SOA_COLUMN(MassXi, massxi, float);
DECLARE_SOA_COLUMN(MassXiKF, massxikf, float);
DECLARE_SOA_COLUMN(MassOmega, massomega, float);
DECLARE_SOA_COLUMN(MassOmegaKF, massomegakf, float);
DECLARE_SOA_COLUMN(MassLambda, masslambda, float);
DECLARE_SOA_COLUMN(MassLambdaKF, masslambdakf, float);
DECLARE_SOA_COLUMN(DCAXYCascToPV, dcaxycasctopv, float);
DECLARE_SOA_COLUMN(DCAXYCascToPVKF, dcaxycasctopvkf, float);
DECLARE_SOA_COLUMN(DCAZCascToPV, dcazcasctopv, float);
DECLARE_SOA_COLUMN(DCAZCascToPVKF, dcazcasctopvkf, float);
DECLARE_SOA_COLUMN(DCACascDau, dcacascdau, float);
DECLARE_SOA_COLUMN(DCACascDauKF, dcacascdaukf, float);
DECLARE_SOA_COLUMN(DCAV0Dau, dcav0dau, float);
DECLARE_SOA_COLUMN(DCAV0DauKF, dcav0daukf, float);
DECLARE_SOA_COLUMN(DCAPosToPV, dcapostopv, float);
DECLARE_SOA_COLUMN(DCAPosToPVKF, dcapostopvkf, float);
DECLARE_SOA_COLUMN(DCANegToPV, dcanegtopv, float);
DECLARE_SOA_COLUMN(DCANegToPVKF, dcanegtopvkf, float);
DECLARE_SOA_COLUMN(DCABachToPV, dcabachtopv, float);
DECLARE_SOA_COLUMN(DCABachToPVKF, dcabachtopvkf, float);
DECLARE_SOA_COLUMN(PACascToPV, pacasctopv, float);
DECLARE_SOA_COLUMN(PACascToPVKF, pacasctopvkf, float);
DECLARE_SOA_COLUMN(CascRad, cascrad, float);
DECLARE_SOA_COLUMN(CascRadKF, cascradkf, float);
DECLARE_SOA_COLUMN(VtxX, vtxx, float);
DECLARE_SOA_COLUMN(VtxY, vtxy, float);
DECLARE_SOA_COLUMN(VtxZ, vtxz, float);
DECLARE_SOA_COLUMN(VtxXErr, vtxxerr, float);
DECLARE_SOA_COLUMN(VtxYErr, vtxyerr, float);
DECLARE_SOA_COLUMN(VtxZErr, vtxzerr, float);
DECLARE_SOA_COLUMN(VtxXKF, vtxxkf, float);
DECLARE_SOA_COLUMN(VtxYKF, vtxykf, float);
DECLARE_SOA_COLUMN(VtxZKF, vtxzkf, float);
DECLARE_SOA_COLUMN(VtxXErrKF, vtxxerrkf, float);
DECLARE_SOA_COLUMN(VtxYErrKF, vtxyerrkf, float);
DECLARE_SOA_COLUMN(VtxZErrKF, vtxzerrkf, float);
DECLARE_SOA_COLUMN(PAV0ToPV, pav0topv, float);
DECLARE_SOA_COLUMN(PAV0ToPVKF, pav0topvkf, float);
DECLARE_SOA_COLUMN(V0Rad, v0rad, float);
DECLARE_SOA_COLUMN(V0RadKF, v0radkf, float);
DECLARE_SOA_COLUMN(Charge, charge, float);
DECLARE_SOA_COLUMN(IsDCAFitter, isdcafitter, int);
DECLARE_SOA_COLUMN(IsKF, iskf, int);

DECLARE_SOA_COLUMN(PtRec, ptrec, float);
DECLARE_SOA_COLUMN(PtRecKF, ptreckf, float);
DECLARE_SOA_COLUMN(PtGen, ptgen, float);
DECLARE_SOA_COLUMN(VtxXrec, vtxxrec, float);
DECLARE_SOA_COLUMN(VtxYrec, vtxyrec, float);
DECLARE_SOA_COLUMN(VtxZrec, vtxzrec, float);
DECLARE_SOA_COLUMN(VtxXrecErr, vtxxrecerr, float);
DECLARE_SOA_COLUMN(VtxYrecErr, vtxyrecerr, float);
DECLARE_SOA_COLUMN(VtxZrecErr, vtxzrecerr, float);
DECLARE_SOA_COLUMN(VtxXrecKF, vtxxreckf, float);
DECLARE_SOA_COLUMN(VtxYrecKF, vtxyreckf, float);
DECLARE_SOA_COLUMN(VtxZrecKF, vtxzreckf, float);
DECLARE_SOA_COLUMN(VtxXrecErrKF, vtxxrecerrkf, float);
DECLARE_SOA_COLUMN(VtxYrecErrKF, vtxyrecerrkf, float);
DECLARE_SOA_COLUMN(VtxZrecErrKF, vtxzrecerrkf, float);
DECLARE_SOA_COLUMN(VtxXMC, vtxxmc, float);
DECLARE_SOA_COLUMN(VtxYMC, vtxymc, float);
DECLARE_SOA_COLUMN(VtxZMC, vtxzmc, float);
DECLARE_SOA_COLUMN(IsTrueCasc, istruecasc, int);



} // namespace kfStrangenessStudy

DECLARE_SOA_TABLE(CascCand, "AOD", "CASCCAND",
                  kfStrangenessStudy::CollID,
                  kfStrangenessStudy::Pt,
                  kfStrangenessStudy::PtKF,
                  kfStrangenessStudy::MassXi,
                  kfStrangenessStudy::MassXiKF,
                  kfStrangenessStudy::MassOmega,
                  kfStrangenessStudy::MassOmegaKF,
                  kfStrangenessStudy::CascRad,
                  kfStrangenessStudy::CascRadKF,
                  kfStrangenessStudy::VtxX,
                  kfStrangenessStudy::VtxY,
                  kfStrangenessStudy::VtxZ,
                  kfStrangenessStudy::VtxXErr,
                  kfStrangenessStudy::VtxYErr,
                  kfStrangenessStudy::VtxZErr,
                  kfStrangenessStudy::VtxXKF,
                  kfStrangenessStudy::VtxYKF,
                  kfStrangenessStudy::VtxZKF,
                  kfStrangenessStudy::VtxXErrKF,
                  kfStrangenessStudy::VtxYErrKF,
                  kfStrangenessStudy::VtxZErrKF,
                  kfStrangenessStudy::DCAXYCascToPV,
                  kfStrangenessStudy::DCAXYCascToPVKF,
                  kfStrangenessStudy::DCAZCascToPV,
                  kfStrangenessStudy::DCAZCascToPVKF,
                  kfStrangenessStudy::DCACascDau,
                  kfStrangenessStudy::DCACascDauKF,
                  kfStrangenessStudy::PACascToPV,
                  kfStrangenessStudy::PACascToPVKF,
                  kfStrangenessStudy::Charge,
                  kfStrangenessStudy::MassLambda,
                  kfStrangenessStudy::MassLambdaKF,
                  kfStrangenessStudy::V0Rad,
                  kfStrangenessStudy::V0RadKF,
                  kfStrangenessStudy::DCAV0Dau,
                  kfStrangenessStudy::DCAV0DauKF,
                  kfStrangenessStudy::DCAPosToPV,
                  kfStrangenessStudy::DCAPosToPVKF,
                  kfStrangenessStudy::DCANegToPV,
                  kfStrangenessStudy::DCANegToPVKF,
                  kfStrangenessStudy::DCABachToPV,
                  kfStrangenessStudy::DCABachToPVKF,
                  kfStrangenessStudy::PAV0ToPV,
                  kfStrangenessStudy::PAV0ToPVKF,
                  kfStrangenessStudy::IsDCAFitter,
                  kfStrangenessStudy::IsKF);

DECLARE_SOA_TABLE(CascCandMC, "AOD", "CASCCANDMC",
                  kfStrangenessStudy::CollID,
                  kfStrangenessStudy::PtRec,
                  kfStrangenessStudy::PtRecKF,
                  kfStrangenessStudy::PtGen,
                  kfStrangenessStudy::MassXi,
                  kfStrangenessStudy::MassXiKF,
                  kfStrangenessStudy::MassOmega,
                  kfStrangenessStudy::MassOmegaKF,
                  kfStrangenessStudy::CascRad,
                  kfStrangenessStudy::CascRadKF,
                  kfStrangenessStudy::VtxXrec,
                  kfStrangenessStudy::VtxYrec,
                  kfStrangenessStudy::VtxZrec,
                  kfStrangenessStudy::VtxXrecErr,
                  kfStrangenessStudy::VtxYrecErr,
                  kfStrangenessStudy::VtxZrecErr,
                  kfStrangenessStudy::VtxXrecKF,
                  kfStrangenessStudy::VtxYrecKF,
                  kfStrangenessStudy::VtxZrecKF,
                  kfStrangenessStudy::VtxXrecErrKF,
                  kfStrangenessStudy::VtxYrecErrKF,
                  kfStrangenessStudy::VtxZrecErrKF,
                  kfStrangenessStudy::VtxXMC,
                  kfStrangenessStudy::VtxYMC,
                  kfStrangenessStudy::VtxZMC,
                  kfStrangenessStudy::DCAXYCascToPV,
                  kfStrangenessStudy::DCAXYCascToPVKF,
                  kfStrangenessStudy::DCAZCascToPV,
                  kfStrangenessStudy::DCAZCascToPVKF,
                  kfStrangenessStudy::DCACascDau,
                  kfStrangenessStudy::DCACascDauKF,
                  kfStrangenessStudy::PACascToPV,
                  kfStrangenessStudy::PACascToPVKF,
                  kfStrangenessStudy::Charge,
                  kfStrangenessStudy::MassLambda,
                  kfStrangenessStudy::MassLambdaKF,
                  kfStrangenessStudy::V0Rad,
                  kfStrangenessStudy::V0RadKF,
                  kfStrangenessStudy::DCAV0Dau,
                  kfStrangenessStudy::DCAV0DauKF,
                  kfStrangenessStudy::DCAPosToPV,
                  kfStrangenessStudy::DCAPosToPVKF,
                  kfStrangenessStudy::DCANegToPV,
                  kfStrangenessStudy::DCANegToPVKF,
                  kfStrangenessStudy::DCABachToPV,
                  kfStrangenessStudy::DCABachToPVKF,
                  kfStrangenessStudy::PAV0ToPV,
                  kfStrangenessStudy::PAV0ToPVKF,
                  kfStrangenessStudy::IsDCAFitter,
                  kfStrangenessStudy::IsKF,
                  kfStrangenessStudy::IsTrueCasc);
} // namespace o2::aod

#endif // PWGLF_DATAMODEL_KFSTRANGENESSSTUDY_H_
