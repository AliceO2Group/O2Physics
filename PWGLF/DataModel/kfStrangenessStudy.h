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
// DECLARE_SOA_COLUMN(DCAV0ToPV, dcav0topv, float);
// DECLARE_SOA_COLUMN(DCAV0ToPVKF, dcav0topvkf, float);
DECLARE_SOA_COLUMN(DCAV0Dau, dcav0dau, float);
DECLARE_SOA_COLUMN(DCAV0DauKF, dcav0daukf, float);
DECLARE_SOA_COLUMN(DCAPosToPV, dcapostopv, float);
DECLARE_SOA_COLUMN(DCAPosToPVKF, dcapostopvkf, float);
DECLARE_SOA_COLUMN(DCANegToPV, dcanegtopv, float);
DECLARE_SOA_COLUMN(DCANegToPVKF, dcanegtopvkf, float);
DECLARE_SOA_COLUMN(DCABachToPV, dcabachtopv, float);
DECLARE_SOA_COLUMN(DCABachToPVKF, dcabachtopvkf, float);
DECLARE_SOA_COLUMN(CosPACascToPV, cospacasctopv, float);
DECLARE_SOA_COLUMN(CosPACascToPVKF, cospacasctopvkf, float);
DECLARE_SOA_COLUMN(CascRad, cascrad, float);
DECLARE_SOA_COLUMN(CascRadKF, cascradkf, float);
DECLARE_SOA_COLUMN(V0Rad, v0rad, float);
DECLARE_SOA_COLUMN(V0RadKF, v0radkf, float);
DECLARE_SOA_COLUMN(Charge, charge, float);
DECLARE_SOA_COLUMN(IsDCAFitter, isdcafitter, int);
DECLARE_SOA_COLUMN(IsKF, iskf, int);
} // namespace kfStrangenessStudy

DECLARE_SOA_TABLE(CascCand, "AOD", "CASCCAND",
                  kfStrangenessStudy::CollID,
                  kfStrangenessStudy::Pt,
                  kfStrangenessStudy::PtKF,
                  kfStrangenessStudy::MassXi,
                  kfStrangenessStudy::MassXiKF,
                  kfStrangenessStudy::MassOmega,
                  kfStrangenessStudy::MassOmegaKF,
                  kfStrangenessStudy::MassLambda,
                  kfStrangenessStudy::MassLambdaKF,
                  kfStrangenessStudy::DCAXYCascToPV,
                  kfStrangenessStudy::DCAXYCascToPVKF,
                  kfStrangenessStudy::DCAZCascToPV,
                  kfStrangenessStudy::DCAZCascToPVKF,
                  kfStrangenessStudy::DCACascDau,
                  kfStrangenessStudy::DCACascDauKF,
                  // kfStrangenessStudy::DCAV0ToPV,
                  // kfStrangenessStudy::DCAV0ToPVKF,
                  kfStrangenessStudy::DCAV0Dau,
                  kfStrangenessStudy::DCAV0DauKF,
                  kfStrangenessStudy::DCAPosToPV,
                  kfStrangenessStudy::DCAPosToPVKF,
                  kfStrangenessStudy::DCANegToPV,
                  kfStrangenessStudy::DCANegToPVKF,
                  kfStrangenessStudy::DCABachToPV,
                  kfStrangenessStudy::DCABachToPVKF,
                  kfStrangenessStudy::CosPACascToPV,
                  kfStrangenessStudy::CosPACascToPVKF,
                  kfStrangenessStudy::CascRad,
                  kfStrangenessStudy::CascRadKF,
                  kfStrangenessStudy::V0Rad,
                  kfStrangenessStudy::V0RadKF,
                  kfStrangenessStudy::Charge, 
                  kfStrangenessStudy::IsDCAFitter,
                  kfStrangenessStudy::IsKF);
} // namespace o2::aod

#endif // PWGLF_DATAMODEL_KFSTRANGENESSSTUDY_H_
