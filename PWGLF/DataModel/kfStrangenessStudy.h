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

// ********* Cascade **********
// momentum
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(PtKF, ptkf, float);
DECLARE_SOA_COLUMN(PtRec, ptrec, float);
DECLARE_SOA_COLUMN(PtRecKF, ptreckf, float);
DECLARE_SOA_COLUMN(PtGen, ptgen, float);
// mass
DECLARE_SOA_COLUMN(MassXi, massxi, float);
DECLARE_SOA_COLUMN(MassXiKF, massxikf, float);
// DECLARE_SOA_COLUMN(MassOmega, massomega, float);
// DECLARE_SOA_COLUMN(MassOmegaKF, massomegakf, float);
// DCA
DECLARE_SOA_COLUMN(DCAXYCascToPV, dcaxycasctopv, float);
DECLARE_SOA_COLUMN(DCAXYCascToPVKF, dcaxycasctopvkf, float);
DECLARE_SOA_COLUMN(DCAZCascToPV, dcazcasctopv, float);
DECLARE_SOA_COLUMN(DCAZCascToPVKF, dcazcasctopvkf, float);
DECLARE_SOA_COLUMN(DCACascDau, dcacascdau, float);
DECLARE_SOA_COLUMN(DCACascDauKF, dcacascdaukf, float);
// pointing angle
DECLARE_SOA_COLUMN(PACascToPV, pacasctopv, float);
DECLARE_SOA_COLUMN(PACascToPVKF, pacasctopvkf, float);
// radius
DECLARE_SOA_COLUMN(CascRad, cascrad, float);
DECLARE_SOA_COLUMN(CascRadKF, cascradkf, float);
// decay vertex
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
// production vertex
DECLARE_SOA_COLUMN(ProdVtxXMC, prodvtxxmc, float);
DECLARE_SOA_COLUMN(ProdVtxYMC, prodvtxymc, float);
DECLARE_SOA_COLUMN(ProdVtxZMC, prodvtxzmc, float);
// charge
DECLARE_SOA_COLUMN(Charge, charge, float);

// ********* V0 **********
// momentum
DECLARE_SOA_COLUMN(V0Pt, v0pt, float);
DECLARE_SOA_COLUMN(V0PtKF, v0ptkf, float);
DECLARE_SOA_COLUMN(V0PtRec, v0ptrec, float);
DECLARE_SOA_COLUMN(V0PtRecKF, v0ptreckf, float);
DECLARE_SOA_COLUMN(V0PtGen, v0ptgen, float);
// mass
DECLARE_SOA_COLUMN(MassLambda, masslambda, float);
DECLARE_SOA_COLUMN(MassLambdaKF, masslambdakf, float);
// DCA
DECLARE_SOA_COLUMN(DCAV0Dau, dcav0dau, float);
DECLARE_SOA_COLUMN(DCAV0DauKF, dcav0daukf, float);
// pointing angle
DECLARE_SOA_COLUMN(PAV0ToPV, pav0topv, float);
DECLARE_SOA_COLUMN(PAV0ToPVKF, pav0topvkf, float);
// radius
DECLARE_SOA_COLUMN(V0Rad, v0rad, float);
DECLARE_SOA_COLUMN(V0RadKF, v0radkf, float);
// decay vertex
DECLARE_SOA_COLUMN(V0VtxX, v0vtxx, float);
DECLARE_SOA_COLUMN(V0VtxY, v0vtxy, float);
DECLARE_SOA_COLUMN(V0VtxZ, v0vtxz, float);
DECLARE_SOA_COLUMN(V0VtxXErr, v0vtxxerr, float);
DECLARE_SOA_COLUMN(V0VtxYErr, v0vtxyerr, float);
DECLARE_SOA_COLUMN(V0VtxZErr, v0vtxzerr, float);
DECLARE_SOA_COLUMN(V0VtxXKF, v0vtxxkf, float);
DECLARE_SOA_COLUMN(V0VtxYKF, v0vtxykf, float);
DECLARE_SOA_COLUMN(V0VtxZKF, v0vtxzkf, float);
DECLARE_SOA_COLUMN(V0VtxXErrKF, v0vtxxerrkf, float);
DECLARE_SOA_COLUMN(V0VtxYErrKF, v0vtxyerrkf, float);
DECLARE_SOA_COLUMN(V0VtxZErrKF, v0vtxzerrkf, float);
DECLARE_SOA_COLUMN(V0VtxXrec, v0vtxxrec, float);
DECLARE_SOA_COLUMN(V0VtxYrec, v0vtxyrec, float);
DECLARE_SOA_COLUMN(V0VtxZrec, v0vtxzrec, float);
DECLARE_SOA_COLUMN(V0VtxXrecErr, v0vtxxrecerr, float);
DECLARE_SOA_COLUMN(V0VtxYrecErr, v0vtxyrecerr, float);
DECLARE_SOA_COLUMN(V0VtxZrecErr, v0vtxzrecerr, float);
DECLARE_SOA_COLUMN(V0VtxXrecKF, v0vtxxreckf, float);
DECLARE_SOA_COLUMN(V0VtxYrecKF, v0vtxyreckf, float);
DECLARE_SOA_COLUMN(V0VtxZrecKF, v0vtxzreckf, float);
DECLARE_SOA_COLUMN(V0VtxXrecErrKF, v0vtxxrecerrkf, float);
DECLARE_SOA_COLUMN(V0VtxYrecErrKF, v0vtxyrecerrkf, float);
DECLARE_SOA_COLUMN(V0VtxZrecErrKF, v0vtxzrecerrkf, float);
DECLARE_SOA_COLUMN(V0VtxXMC, v0vtxxmc, float);
DECLARE_SOA_COLUMN(V0VtxYMC, v0vtxymc, float);
DECLARE_SOA_COLUMN(V0VtxZMC, v0vtxzmc, float);
// production vertex
DECLARE_SOA_COLUMN(V0ProdVtxXMC, v0prodvtxxmc, float);
DECLARE_SOA_COLUMN(V0ProdVtxYMC, v0prodvtxymc, float);
DECLARE_SOA_COLUMN(V0ProdVtxZMC, v0prodvtxzmc, float);

// ********* Daughters **********
// DCA
DECLARE_SOA_COLUMN(DCAProtonToPV, dcaprotontopv, float);
DECLARE_SOA_COLUMN(DCAProtonToPVKF, dcaprotontopvkf, float);
DECLARE_SOA_COLUMN(DCAPionToPV, dcapiontopv, float);
DECLARE_SOA_COLUMN(DCAPionToPVKF, dcapiontopvkf, float);
DECLARE_SOA_COLUMN(DCABachToPV, dcabachtopv, float);
DECLARE_SOA_COLUMN(DCABachToPVKF, dcabachtopvkf, float);
// momentum
DECLARE_SOA_COLUMN(PxProtonIURec, pxprotoniurec, float);
DECLARE_SOA_COLUMN(PyProtonIURec, pyprotoniurec, float);
DECLARE_SOA_COLUMN(PzProtonIURec, pzprotoniurec, float);
DECLARE_SOA_COLUMN(PxProtonIURecErr, pxprotoniurecerr, float);
DECLARE_SOA_COLUMN(PyProtonIURecErr, pyprotoniurecerr, float);
DECLARE_SOA_COLUMN(PzProtonIURecErr, pzprotoniurecerr, float);
DECLARE_SOA_COLUMN(PxPionIURec, pxpioniurec, float);
DECLARE_SOA_COLUMN(PyPionIURec, pypioniurec, float);
DECLARE_SOA_COLUMN(PzPionIURec, pzpioniurec, float);
DECLARE_SOA_COLUMN(PxPionIURecErr, pxpioniurecerr, float);
DECLARE_SOA_COLUMN(PyPionIURecErr, pypioniurecerr, float);
DECLARE_SOA_COLUMN(PzPionIURecErr, pzpioniurecerr, float);
DECLARE_SOA_COLUMN(PxProtonRec, pxprotonrec, float);
DECLARE_SOA_COLUMN(PyProtonRec, pyprotonrec, float);
DECLARE_SOA_COLUMN(PzProtonRec, pzprotonrec, float);
DECLARE_SOA_COLUMN(PxProtonRecErr, pxprotonrecerr, float);
DECLARE_SOA_COLUMN(PyProtonRecErr, pyprotonrecerr, float);
DECLARE_SOA_COLUMN(PzProtonRecErr, pzprotonrecerr, float);
DECLARE_SOA_COLUMN(PxPionRec, pxpionrec, float);
DECLARE_SOA_COLUMN(PyPionRec, pypionrec, float);
DECLARE_SOA_COLUMN(PzPionRec, pzpionrec, float);
DECLARE_SOA_COLUMN(PxPionRecErr, pxpionrecerr, float);
DECLARE_SOA_COLUMN(PyPionRecErr, pypionrecerr, float);
DECLARE_SOA_COLUMN(PzPionRecErr, pzpionrecerr, float);
DECLARE_SOA_COLUMN(PxProtonMC, pxprotonmc, float);
DECLARE_SOA_COLUMN(PyProtonMC, pyprotonmc, float);
DECLARE_SOA_COLUMN(PzProtonMC, pzprotonmc, float);
DECLARE_SOA_COLUMN(PxPionMC, pxpionmc, float);
DECLARE_SOA_COLUMN(PyPionMC, pypionmc, float);
DECLARE_SOA_COLUMN(PzPionMC, pzpionmc, float);
// position
DECLARE_SOA_COLUMN(XProtonRec, xprotonrec, float);
DECLARE_SOA_COLUMN(YProtonRec, yprotonrec, float);
DECLARE_SOA_COLUMN(ZProtonRec, zprotonrec, float);
DECLARE_SOA_COLUMN(XProtonRecErr, xprotonrecerr, float);
DECLARE_SOA_COLUMN(YProtonRecErr, yprotonrecerr, float);
DECLARE_SOA_COLUMN(ZProtonRecErr, zprotonrecerr, float);
DECLARE_SOA_COLUMN(XPionRec, xpionrec, float);
DECLARE_SOA_COLUMN(YPionRec, ypionrec, float);
DECLARE_SOA_COLUMN(ZPionRec, zpionrec, float);
DECLARE_SOA_COLUMN(XPionRecErr, xpionrecerr, float);
DECLARE_SOA_COLUMN(YPionRecErr, ypionrecerr, float);
DECLARE_SOA_COLUMN(ZPionRecErr, zpionrecerr, float);
DECLARE_SOA_COLUMN(XProtonIURec, xprotoniurec, float);
DECLARE_SOA_COLUMN(YProtonIURec, yprotoniurec, float);
DECLARE_SOA_COLUMN(ZProtonIURec, zprotoniurec, float);
DECLARE_SOA_COLUMN(XPionIURec, xpioniurec, float);
DECLARE_SOA_COLUMN(YPionIURec, ypioniurec, float);
DECLARE_SOA_COLUMN(ZPionIURec, zpioniurec, float);
DECLARE_SOA_COLUMN(XBachIURec, xbachiurec, float);
DECLARE_SOA_COLUMN(YBachIURec, ybachiurec, float);
DECLARE_SOA_COLUMN(ZBachIURec, zbachiurec, float);

// eta
DECLARE_SOA_COLUMN(EtaProton, etaproton, float);
DECLARE_SOA_COLUMN(EtaPion, etapion, float);
// TPC clusters
DECLARE_SOA_COLUMN(TPCNClsProton, tpcnclsproton, float);
DECLARE_SOA_COLUMN(TPCNClsPion, tpcnclspion, float);

// ********* MC info **********
DECLARE_SOA_COLUMN(IsTrueCasc, istruecasc, int);
DECLARE_SOA_COLUMN(Source, source, int);

DECLARE_SOA_COLUMN(IsDCAFitter, isdcafitter, int);
DECLARE_SOA_COLUMN(IsKF, iskf, int);

} // namespace kfStrangenessStudy

DECLARE_SOA_TABLE(CascCand, "AOD", "CASCCAND",
                  kfStrangenessStudy::CollID,
                  kfStrangenessStudy::Pt,
                  kfStrangenessStudy::PtKF,
                  kfStrangenessStudy::MassXi,
                  kfStrangenessStudy::MassXiKF,
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
                  kfStrangenessStudy::V0Pt,
                  kfStrangenessStudy::V0PtKF,
                  kfStrangenessStudy::MassLambda,
                  kfStrangenessStudy::MassLambdaKF,
                  kfStrangenessStudy::V0Rad,
                  kfStrangenessStudy::V0RadKF,
                  kfStrangenessStudy::V0VtxX,
                  kfStrangenessStudy::V0VtxY,
                  kfStrangenessStudy::V0VtxZ,
                  kfStrangenessStudy::V0VtxXErr,
                  kfStrangenessStudy::V0VtxYErr,
                  kfStrangenessStudy::V0VtxZErr,
                  kfStrangenessStudy::V0VtxXKF,
                  kfStrangenessStudy::V0VtxYKF,
                  kfStrangenessStudy::V0VtxZKF,
                  kfStrangenessStudy::V0VtxXErrKF,
                  kfStrangenessStudy::V0VtxYErrKF,
                  kfStrangenessStudy::V0VtxZErrKF,
                  kfStrangenessStudy::DCAV0Dau,
                  kfStrangenessStudy::DCAV0DauKF,
                  kfStrangenessStudy::DCAProtonToPV,
                  kfStrangenessStudy::DCAProtonToPVKF,
                  kfStrangenessStudy::DCAPionToPV,
                  kfStrangenessStudy::DCAPionToPVKF,
                  kfStrangenessStudy::DCABachToPV,
                  kfStrangenessStudy::DCABachToPVKF,
                  kfStrangenessStudy::PAV0ToPV,
                  kfStrangenessStudy::PAV0ToPVKF,
                  kfStrangenessStudy::PxProtonIURec,
                  kfStrangenessStudy::PyProtonIURec,
                  kfStrangenessStudy::PzProtonIURec,
                  kfStrangenessStudy::PxProtonIURecErr,
                  kfStrangenessStudy::PyProtonIURecErr,
                  kfStrangenessStudy::PzProtonIURecErr,
                  kfStrangenessStudy::PxPionIURec,
                  kfStrangenessStudy::PyPionIURec,
                  kfStrangenessStudy::PzPionIURec,
                  kfStrangenessStudy::PxPionIURecErr,
                  kfStrangenessStudy::PyPionIURecErr,
                  kfStrangenessStudy::PzPionIURecErr,
                  kfStrangenessStudy::PxProtonRec,
                  kfStrangenessStudy::PyProtonRec,
                  kfStrangenessStudy::PzProtonRec,
                  kfStrangenessStudy::PxProtonRecErr,
                  kfStrangenessStudy::PyProtonRecErr,
                  kfStrangenessStudy::PzProtonRecErr,
                  kfStrangenessStudy::PxPionRec,
                  kfStrangenessStudy::PyPionRec,
                  kfStrangenessStudy::PzPionRec,
                  kfStrangenessStudy::PxPionRecErr,
                  kfStrangenessStudy::PyPionRecErr,
                  kfStrangenessStudy::PzPionRecErr,
                  kfStrangenessStudy::XProtonIURec,
                  kfStrangenessStudy::YProtonIURec,
                  kfStrangenessStudy::ZProtonIURec,
                  kfStrangenessStudy::XPionIURec,
                  kfStrangenessStudy::YPionIURec,
                  kfStrangenessStudy::ZPionIURec,
                  kfStrangenessStudy::XBachIURec,
                  kfStrangenessStudy::YBachIURec,
                  kfStrangenessStudy::ZBachIURec,
                  kfStrangenessStudy::EtaProton,
                  kfStrangenessStudy::EtaPion,
                  kfStrangenessStudy::TPCNClsProton,
                  kfStrangenessStudy::TPCNClsPion,
                  kfStrangenessStudy::IsDCAFitter,
                  kfStrangenessStudy::IsKF);

DECLARE_SOA_TABLE(CascCandMC, "AOD", "CASCCANDMC",
                  kfStrangenessStudy::CollID,
                  kfStrangenessStudy::PtRec,
                  kfStrangenessStudy::PtRecKF,
                  kfStrangenessStudy::PtGen,
                  kfStrangenessStudy::MassXi,
                  kfStrangenessStudy::MassXiKF,
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
                  kfStrangenessStudy::ProdVtxXMC,
                  kfStrangenessStudy::ProdVtxYMC,
                  kfStrangenessStudy::ProdVtxZMC,
                  kfStrangenessStudy::DCAXYCascToPV,
                  kfStrangenessStudy::DCAXYCascToPVKF,
                  kfStrangenessStudy::DCAZCascToPV,
                  kfStrangenessStudy::DCAZCascToPVKF,
                  kfStrangenessStudy::DCACascDau,
                  kfStrangenessStudy::DCACascDauKF,
                  kfStrangenessStudy::PACascToPV,
                  kfStrangenessStudy::PACascToPVKF,
                  kfStrangenessStudy::Charge,
                  kfStrangenessStudy::V0PtRec,
                  kfStrangenessStudy::V0PtRecKF,
                  kfStrangenessStudy::V0PtGen,
                  kfStrangenessStudy::MassLambda,
                  kfStrangenessStudy::MassLambdaKF,
                  kfStrangenessStudy::V0Rad,
                  kfStrangenessStudy::V0RadKF,
                  kfStrangenessStudy::V0VtxXrec,
                  kfStrangenessStudy::V0VtxYrec,
                  kfStrangenessStudy::V0VtxZrec,
                  kfStrangenessStudy::V0VtxXrecErr,
                  kfStrangenessStudy::V0VtxYrecErr,
                  kfStrangenessStudy::V0VtxZrecErr,
                  kfStrangenessStudy::V0VtxXrecKF,
                  kfStrangenessStudy::V0VtxYrecKF,
                  kfStrangenessStudy::V0VtxZrecKF,
                  kfStrangenessStudy::V0VtxXrecErrKF,
                  kfStrangenessStudy::V0VtxYrecErrKF,
                  kfStrangenessStudy::V0VtxZrecErrKF,
                  kfStrangenessStudy::V0VtxXMC,
                  kfStrangenessStudy::V0VtxYMC,
                  kfStrangenessStudy::V0VtxZMC,
                  kfStrangenessStudy::DCAV0Dau,
                  kfStrangenessStudy::DCAV0DauKF,
                  kfStrangenessStudy::DCAProtonToPV,
                  kfStrangenessStudy::DCAProtonToPVKF,
                  kfStrangenessStudy::DCAPionToPV,
                  kfStrangenessStudy::DCAPionToPVKF,
                  kfStrangenessStudy::DCABachToPV,
                  kfStrangenessStudy::DCABachToPVKF,
                  kfStrangenessStudy::PAV0ToPV,
                  kfStrangenessStudy::PAV0ToPVKF,
                  kfStrangenessStudy::PxProtonIURec,
                  kfStrangenessStudy::PyProtonIURec,
                  kfStrangenessStudy::PzProtonIURec,
                  kfStrangenessStudy::PxProtonIURecErr,
                  kfStrangenessStudy::PyProtonIURecErr,
                  kfStrangenessStudy::PzProtonIURecErr,
                  kfStrangenessStudy::PxPionIURec,
                  kfStrangenessStudy::PyPionIURec,
                  kfStrangenessStudy::PzPionIURec,
                  kfStrangenessStudy::PxPionIURecErr,
                  kfStrangenessStudy::PyPionIURecErr,
                  kfStrangenessStudy::PzPionIURecErr,
                  kfStrangenessStudy::PxProtonRec,
                  kfStrangenessStudy::PyProtonRec,
                  kfStrangenessStudy::PzProtonRec,
                  kfStrangenessStudy::PxProtonRecErr,
                  kfStrangenessStudy::PyProtonRecErr,
                  kfStrangenessStudy::PzProtonRecErr,
                  kfStrangenessStudy::PxPionRec,
                  kfStrangenessStudy::PyPionRec,
                  kfStrangenessStudy::PzPionRec,
                  kfStrangenessStudy::PxPionRecErr,
                  kfStrangenessStudy::PyPionRecErr,
                  kfStrangenessStudy::PzPionRecErr,
                  kfStrangenessStudy::PxProtonMC,
                  kfStrangenessStudy::PyProtonMC,
                  kfStrangenessStudy::PzProtonMC,
                  kfStrangenessStudy::PxPionMC,
                  kfStrangenessStudy::PyPionMC,
                  kfStrangenessStudy::PzPionMC,
                  kfStrangenessStudy::XProtonRec,
                  kfStrangenessStudy::YProtonRec,
                  kfStrangenessStudy::ZProtonRec,
                  kfStrangenessStudy::XProtonRecErr,
                  kfStrangenessStudy::YProtonRecErr,
                  kfStrangenessStudy::ZProtonRecErr,
                  kfStrangenessStudy::XPionRec,
                  kfStrangenessStudy::YPionRec,
                  kfStrangenessStudy::ZPionRec,
                  kfStrangenessStudy::XPionRecErr,
                  kfStrangenessStudy::YPionRecErr,
                  kfStrangenessStudy::ZPionRecErr,
                  kfStrangenessStudy::XProtonIURec,
                  kfStrangenessStudy::YProtonIURec,
                  kfStrangenessStudy::ZProtonIURec,
                  kfStrangenessStudy::XPionIURec,
                  kfStrangenessStudy::YPionIURec,
                  kfStrangenessStudy::ZPionIURec,
                  kfStrangenessStudy::XBachIURec,
                  kfStrangenessStudy::YBachIURec,
                  kfStrangenessStudy::ZBachIURec,
                  kfStrangenessStudy::EtaProton,
                  kfStrangenessStudy::EtaPion,
                  kfStrangenessStudy::TPCNClsProton,
                  kfStrangenessStudy::TPCNClsPion,
                  kfStrangenessStudy::IsDCAFitter,
                  kfStrangenessStudy::IsKF,
                  kfStrangenessStudy::IsTrueCasc,
                  kfStrangenessStudy::Source);
} // namespace o2::aod

#endif // PWGLF_DATAMODEL_KFSTRANGENESSSTUDY_H_
