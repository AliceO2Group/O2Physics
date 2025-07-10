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

/// \file treeCreatorXic0XicpToHadronic.cxx
/// \brief tree creator of Xic0 and Xicp candidates 

/// \author Jinhyun Park <jinhyun.park@cern.ch>, Pusan National University
/// \author Krista Smith <krista.lizbeth.smith@cern.ch>, Pusan National University

// Mandatory
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

// STL
#include <vector>

// AOB
#include "CommonConstants/PhysicsConstants.h"

// DataModel
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

namespace o2::aod
{

namespace full
{
DECLARE_SOA_COLUMN(CandidateSelFlag, candidateSelFlag, int); //! Selection flag of candidate from candidate selector output
// vertices
DECLARE_SOA_COLUMN(Chi2Sv, chi2Sv, float);
DECLARE_SOA_COLUMN(Chi2XiVtx, chi2XiVtx, float);
DECLARE_SOA_COLUMN(Chi2LamVtx, chi2LamVtx, float);
// properties of Xic0
DECLARE_SOA_COLUMN(E, e, float);
DECLARE_SOA_COLUMN(M, m, float);
DECLARE_SOA_COLUMN(P, p, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Ct, ct, float);
DECLARE_SOA_COLUMN(DecayLength, decayLength, float);
DECLARE_SOA_COLUMN(DecayLengthXY, decayLengthXY, float);
DECLARE_SOA_COLUMN(DecayLengthNormalised, decayLengthNormalised, float);
DECLARE_SOA_COLUMN(DecayLengthXYNormalised, decayLengthXYNormalised, float);
DECLARE_SOA_COLUMN(Cpa, cpa, float);
DECLARE_SOA_COLUMN(CpaXY, cpaXY, float);
DECLARE_SOA_COLUMN(Chi2TopoXic0ToPVBeforeConstraint, chi2TopoXic0ToPVBeforeConstraint, float);
DECLARE_SOA_COLUMN(Chi2TopoXic0ToPV, chi2TopoXic0ToPV, float);
DECLARE_SOA_COLUMN(Chi2TopoXiToXic0BeforeConstraint, chi2TopoXiToXic0BeforeConstraint, float);
DECLARE_SOA_COLUMN(Chi2TopoXiToXic0, chi2TopoXiToXic0, float);
// properties of Xic0 daughter tracks
DECLARE_SOA_COLUMN(PtXi, ptXi, float);													//! Transverse momentum of Xi(Prong0)
DECLARE_SOA_COLUMN(ImpactParameterXi, impactParameterXi, float);						//! Impact parameter of Xi(prong0)
DECLARE_SOA_COLUMN(ImpactParameterNormalisedXi, impactParameterNormalisedXi, float);	//! Normalised impact parameter of Xi(Prong0)
DECLARE_SOA_COLUMN(PPi, pPi, float);													
DECLARE_SOA_COLUMN(PtPi, ptPi, float);													//! Transverse momentum of Pi(Prong1)
DECLARE_SOA_COLUMN(ImpactParameterPi, impactParameterPi, float);						//! Impact parameter of Pi(Prong1)
DECLARE_SOA_COLUMN(ImpactParameterNormalisedPi, impactParameterNormalisedPi, float);	//! Normalised impact parameter of Pi(Prong1)
DECLARE_SOA_COLUMN(MaxNormalisedDeltaIP, maxNormalisedDeltaIP, float);
DECLARE_SOA_COLUMN(CpaXi, cpaXi, float);
DECLARE_SOA_COLUMN(CpaXYXi, cpaXYXi, float);
DECLARE_SOA_COLUMN(CpaLam, cpaLam, float);
DECLARE_SOA_COLUMN(CpaXYLam, cpaXYLam, float);
DECLARE_SOA_COLUMN(CpaLamToXi, cpaLamToXi, float);
DECLARE_SOA_COLUMN(CpaXYLamToXi, cpaXYLamToXi, float);
DECLARE_SOA_COLUMN(DcaPiXi, dcaPiXi, float);	//! DCA of Xi and charm bachelor pion
DECLARE_SOA_COLUMN(DcaXYPiXi, dcaXYPiXi, float);//! DCAXY of Xi and charm bachelor pion
DECLARE_SOA_COLUMN(InvMassXi, invMassXi, float);
DECLARE_SOA_COLUMN(InvMassLambda, invMassLambda, float);
DECLARE_SOA_COLUMN(PBachelorPi, pBachelorPi, float);
DECLARE_SOA_COLUMN(PPiFromLambda, pPiFromLambda, float);
DECLARE_SOA_COLUMN(PPrFromLambda, pPrFromLambda, float);
// residuals and pulls for Xic0
DECLARE_SOA_COLUMN(PtResidual, ptResidual, float);
DECLARE_SOA_COLUMN(PResidual, pResidual, float);
DECLARE_SOA_COLUMN(XPvResidual, xPvResidual, float);
DECLARE_SOA_COLUMN(YPvResidual, yPvResidual, float);
DECLARE_SOA_COLUMN(ZPvResidual, zPvResidual, float);
DECLARE_SOA_COLUMN(XPvPull, xPvPull, float);
DECLARE_SOA_COLUMN(YPvPull, yPvPull, float);
DECLARE_SOA_COLUMN(ZPvPull, zPvPull, float);
DECLARE_SOA_COLUMN(XSvResidual, xSvResidual, float);
DECLARE_SOA_COLUMN(YSvResidual, ySvResidual, float);
DECLARE_SOA_COLUMN(ZSvResidual, zSvResidual, float);
DECLARE_SOA_COLUMN(XSvPull, xSvPull, float);
DECLARE_SOA_COLUMN(YSvPull, ySvPull, float);
DECLARE_SOA_COLUMN(ZSvPull, zSvPull, float);

// Xicp specific columns
DECLARE_SOA_COLUMN(Sign, sign, float);
DECLARE_SOA_COLUMN(Chi2TopoXicPlusToPVBeforeConstraint, chi2TopoXicPlusToPVBeforeConstraint, float);
DECLARE_SOA_COLUMN(Chi2TopoXicPlusToPV, chi2TopoXicPlusToPV, float);
//DECLARE_SOA_COLUMN(Chi2TopoXiToXicPlusBeforeConstraint, chi2TopoXiToXicPlusBeforeConstraint, float);
//DECLARE_SOA_COLUMN(Chi2TopoXiToXicPlus, chi2TopoXiToXicPlus, float);
DECLARE_SOA_COLUMN(PPi0, pPi0, float);
DECLARE_SOA_COLUMN(PtPi0, ptPi0, float);                                               //! Transverse momentum of Pi0 (prong1) (GeV/c)
DECLARE_SOA_COLUMN(ImpactParameterPi0, impactParameterPi0, float);                     //! Impact parameter of Pi0 (prong1)
DECLARE_SOA_COLUMN(ImpactParameterNormalisedPi0, impactParameterNormalisedPi0, float); //! Normalised impact parameter of Pi0 (prong1)
DECLARE_SOA_COLUMN(PPi1, pPi1, float);
DECLARE_SOA_COLUMN(PtPi1, ptPi1, float);                                               //! Transverse momentum of Pi1 (prong2) (GeV/c)
DECLARE_SOA_COLUMN(ImpactParameterPi1, impactParameterPi1, float);                     //! Normalised impact parameter of Pi1 (prong2)
DECLARE_SOA_COLUMN(ImpactParameterNormalisedPi1, impactParameterNormalisedPi1, float); //! Normalised impact parameter of Pi1 (prong2)
DECLARE_SOA_COLUMN(DcaXYPi0Pi1, dcaXYPi0Pi1, float);
DECLARE_SOA_COLUMN(DcaXYPi0Xi, dcaXYPi0Xi, float);
DECLARE_SOA_COLUMN(DcaXYPi1Xi, dcaXYPi1Xi, float);
DECLARE_SOA_COLUMN(DcaPi0Pi1, dcaPi0Pi1, float);
DECLARE_SOA_COLUMN(DcaPi0Xi, dcaPi0Xi, float);
DECLARE_SOA_COLUMN(DcaPi1Xi, dcaPi1Xi, float);
DECLARE_SOA_COLUMN(InvMassXiPi0, invMassXiPi0, float);
DECLARE_SOA_COLUMN(InvMassXiPi1, invMassXiPi1, float);
}

DECLARE_SOA_TABLE(HfCandXic0Lites, "AOD", "HFCANDXIC0LITE",
				  hf_cand_xic0_xicp_to_hadronic::xic0::FlagMcMatchRec,
				  hf_cand_xic0_xicp_to_hadronic::xic0::DebugMcRec,
				  hf_cand_xic0_xicp_to_hadronic::xic0::OriginRec,
				  full::CandidateSelFlag,
				  full::Y,
				  full::Eta,
				  full::Phi,
				  full::P,
				  full::Pt,
				  full::PtXi,
				  full::PtPi,
				  full::M,
				  full::InvMassXi,
				  full::InvMassLambda,
				  full::Chi2Sv,
				  full::Ct,
				  full::DecayLength,
				  full::DecayLengthNormalised,
				  full::DecayLengthXY,
				  full::DecayLengthXYNormalised,
				  full::Cpa,
				  full::CpaXY,
				  full::CpaXi,
				  full::CpaXYXi,
				  full::CpaLam,
				  full::CpaXYLam,
				  full::ImpactParameterXi,
				  full::ImpactParameterNormalisedXi,
				  full::ImpactParameterPi,
				  full::ImpactParameterNormalisedPi,
				  full::MaxNormalisedDeltaIP);

DECLARE_SOA_TABLE(HfCandXic0LiteKfs, "AOD", "HFCANDXIC0LITKF",
				  hf_cand_xic0_xicp_to_hadronic::xic0::FlagMcMatchRec,
				  hf_cand_xic0_xicp_to_hadronic::xic0::DebugMcRec,
				  hf_cand_xic0_xicp_to_hadronic::xic0::OriginRec,
				  full::CandidateSelFlag,
				  full::Y,
				  full::Eta,
				  full::Phi,
				  full::P,
				  full::Pt,
				  full::PtXi,
				  full::PtPi,
				  full::M,
				  full::InvMassXi,
				  full::InvMassLambda,
				  full::Chi2Sv,
				  full::Ct,
				  full::DecayLength,
				  full::DecayLengthNormalised,
				  full::DecayLengthXY,
				  full::DecayLengthXYNormalised,
				  full::Cpa,
				  full::CpaXY,
				  full::CpaXi,
				  full::CpaXYXi,
				  full::CpaLam,
				  full::CpaXYLam,
				  full::ImpactParameterXi,
				  full::ImpactParameterNormalisedXi,
				  full::ImpactParameterPi,
				  full::ImpactParameterNormalisedPi,
				  full::MaxNormalisedDeltaIP,
				  // KF specific columns
				  full::Chi2XiVtx,
				  full::Chi2LamVtx,
				  full::Chi2TopoXic0ToPVBeforeConstraint,
				  full::Chi2TopoXic0ToPV,
				  full::Chi2TopoXiToXic0BeforeConstraint,
				  full::Chi2TopoXiToXic0,
				  full::DcaXYPiXi,
				  full::DcaPiXi);

DECLARE_SOA_TABLE(HfCandXic0Fulls, "AOD", "HFCANDXIC0FULL",
				  hf_cand_xic0_xicp_to_hadronic::xic0::FlagMcMatchRec,
				  hf_cand_xic0_xicp_to_hadronic::xic0::DebugMcRec,
				  hf_cand_xic0_xicp_to_hadronic::xic0::OriginRec,
				  full::CandidateSelFlag,
				  full::Y,
				  full::Eta,
				  full::Phi,
				  full::P,
				  full::Pt,
				  full::PtXi,
				  full::PtPi,
				  full::M,
				  full::InvMassXi,
				  full::InvMassLambda,
				  full::Chi2Sv,
				  full::Ct,
				  full::DecayLength,
				  full::DecayLengthNormalised,
				  full::DecayLengthXY,
				  full::DecayLengthXYNormalised,
				  full::Cpa,
				  full::CpaXY,
				  full::CpaXi,
				  full::CpaXYXi,
				  full::CpaLam,
				  full::CpaXYLam,
				  full::ImpactParameterXi,
				  full::ImpactParameterNormalisedXi,
				  full::ImpactParameterPi,
				  full::ImpactParameterNormalisedPi,
				  full::MaxNormalisedDeltaIP,
				  // additional columns only stored in the full candidate table
				  full::CpaLamToXi,
				  full::CpaXYLamToXi,
				  full::PPi,
				  full::PBachelorPi,
				  full::PPiFromLambda,
				  full::PPrFromLambda,
				  hf_cand_xic0_xicp_to_hadronic::xic0::DcaXiDaughters,
				  hf_cand_xic0_xicp_to_hadronic::xic0::DcaV0Daughters,
				  hf_cand_xic0_xicp_to_hadronic::xic0::DcaPosToPV,
				  hf_cand_xic0_xicp_to_hadronic::xic0::DcaNegToPV,
				  hf_cand_xic0_xicp_to_hadronic::xic0::DcaBachelorToPV,
				  hf_cand_xic0_xicp_to_hadronic::xic0::DcaXYCascToPV,
				  hf_cand_xic0_xicp_to_hadronic::xic0::DcaZCascToPV,
				  hf_cand_xic0_xicp_to_hadronic::xic0::NSigTpcPiFromXic0,
				  hf_cand_xic0_xicp_to_hadronic::xic0::NSigTpcBachelorPi,
				  hf_cand_xic0_xicp_to_hadronic::xic0::NSigTpcPiFromLambda,
				  hf_cand_xic0_xicp_to_hadronic::xic0::NSigTpcPrFromLambda,
				  hf_cand_xic0_xicp_to_hadronic::xic0::NSigTofPiFromXic0,
				  hf_cand_xic0_xicp_to_hadronic::xic0::NSigTofBachelorPi,
				  hf_cand_xic0_xicp_to_hadronic::xic0::NSigTofPiFromLambda,
				  hf_cand_xic0_xicp_to_hadronic::xic0::NSigTofPrFromLambda);

DECLARE_SOA_TABLE(HfCandXic0FullKfs, "AOD", "HFCANDXIC0FULKF",
				  hf_cand_xic0_xicp_to_hadronic::xic0::FlagMcMatchRec,
				  hf_cand_xic0_xicp_to_hadronic::xic0::DebugMcRec,
				  hf_cand_xic0_xicp_to_hadronic::xic0::OriginRec,
				  full::CandidateSelFlag,
				  full::Y,
				  full::Eta,
				  full::Phi,
				  full::P,
				  full::Pt,
				  full::PtXi,
				  full::PtPi,
				  full::M,
				  full::InvMassXi,
				  full::InvMassLambda,
				  full::Chi2Sv,
				  full::Ct,
				  full::DecayLength,
				  full::DecayLengthNormalised,
				  full::DecayLengthXY,
				  full::DecayLengthXYNormalised,
				  full::Cpa,
				  full::CpaXY,
				  full::CpaXi,
				  full::CpaXYXi,
				  full::CpaLam,
				  full::CpaXYLam,
				  full::ImpactParameterXi,
				  full::ImpactParameterNormalisedXi,
				  full::ImpactParameterPi,
				  full::ImpactParameterNormalisedPi,
				  full::MaxNormalisedDeltaIP,
				  // additional columns only stored in the full candidate table
				  full::CpaLamToXi,
				  full::CpaXYLamToXi,
				  full::PPi,
				  full::PBachelorPi,
				  full::PPiFromLambda,
				  full::PPrFromLambda,
				  hf_cand_xic0_xicp_to_hadronic::xic0::DcaXiDaughters,
				  hf_cand_xic0_xicp_to_hadronic::xic0::DcaV0Daughters,
				  hf_cand_xic0_xicp_to_hadronic::xic0::DcaPosToPV,
				  hf_cand_xic0_xicp_to_hadronic::xic0::DcaNegToPV,
				  hf_cand_xic0_xicp_to_hadronic::xic0::DcaBachelorToPV,
				  hf_cand_xic0_xicp_to_hadronic::xic0::DcaXYCascToPV,
				  hf_cand_xic0_xicp_to_hadronic::xic0::DcaZCascToPV,
				  hf_cand_xic0_xicp_to_hadronic::xic0::NSigTpcPiFromXic0,
				  hf_cand_xic0_xicp_to_hadronic::xic0::NSigTpcBachelorPi,
				  hf_cand_xic0_xicp_to_hadronic::xic0::NSigTpcPiFromLambda,
				  hf_cand_xic0_xicp_to_hadronic::xic0::NSigTpcPrFromLambda,
				  hf_cand_xic0_xicp_to_hadronic::xic0::NSigTofPiFromXic0,
				  hf_cand_xic0_xicp_to_hadronic::xic0::NSigTofBachelorPi,
				  hf_cand_xic0_xicp_to_hadronic::xic0::NSigTofPiFromLambda,
				  hf_cand_xic0_xicp_to_hadronic::xic0::NSigTofPrFromLambda,
				  // KF specific columns
				  full::Chi2XiVtx,
				  full::Chi2LamVtx,
				  full::Chi2TopoXic0ToPVBeforeConstraint,
				  full::Chi2TopoXic0ToPV,
				  full::Chi2TopoXiToXic0BeforeConstraint,
				  full::Chi2TopoXiToXic0,
				  full::DcaXYPiXi,
				  full::DcaPiXi);

DECLARE_SOA_TABLE(HfCandXic0ToXiPiFullPs, "AOD", "HFXIC0TOXIPIFULLP",	
				  hf_cand_xic0_xicp_to_hadronic::xic0::FlagMcMatchGen,
				  hf_cand_xic0_xicp_to_hadronic::xic0::DebugMcGen,
				  hf_cand_xic0_xicp_to_hadronic::xic0::OriginGen,
				  full::Pt,
				  full::Eta,
				  full::Phi,
				  full::Y);

DECLARE_SOA_TABLE(HfCandXic0ToXiPiResiduals, "AOD", "HFXIC0TOXIPIRESID",
				  hf_cand_xic0_xicp_to_hadronic::xic0::OriginGen,
				  full::PResidual,
				  full::PtResidual,
				  full::XPvResidual,
				  full::YPvResidual,
				  full::ZPvResidual,
				  full::XPvPull,
				  full::YPvPull,
				  full::ZPvPull,
				  full::XSvResidual,
				  full::YSvResidual,
				  full::ZSvResidual,
				  full::XSvPull,
				  full::YSvPull,
				  full::ZSvPull);

DECLARE_SOA_TABLE(HfCandXicToXiPiPiLites, "AOD", "HFXICXI2PILITE",
                  hf_cand_xic0_xicp_to_hadronic::xicp::FlagMcMatchRec,
                  hf_cand_xic0_xicp_to_hadronic::xicp::OriginRec,
                  full::CandidateSelFlag,
                  full::Sign,
                  full::Y,
                  full::Eta,
                  full::Phi,
                  full::P,
                  full::Pt,
                  full::PtXi,
                  full::PtPi0,
                  full::PtPi1,
                  full::M,
                  full::InvMassXi,
                  full::InvMassLambda,
                  full::InvMassXiPi0,
                  full::InvMassXiPi1,
                  full::Chi2Sv,
                  full::Ct,
                  full::DecayLength,
                  full::DecayLengthNormalised,
                  full::DecayLengthXY,
                  full::DecayLengthXYNormalised,
                  full::Cpa,
                  full::CpaXY,
                  full::CpaXi,
                  full::CpaXYXi,
                  full::CpaLam,
                  full::CpaXYLam,
                  full::ImpactParameterXi,
                  full::ImpactParameterNormalisedXi,
                  full::ImpactParameterPi0,
                  full::ImpactParameterNormalisedPi0,
                  full::ImpactParameterPi1,
                  full::ImpactParameterNormalisedPi1,
                  full::MaxNormalisedDeltaIP);

DECLARE_SOA_TABLE(HfCandXicToXiPiPiLiteKfs, "AOD", "HFXICXI2PILITKF",
                  hf_cand_xic0_xicp_to_hadronic::xicp::FlagMcMatchRec,
                  hf_cand_xic0_xicp_to_hadronic::xicp::OriginRec,
                  full::CandidateSelFlag,
                  full::Sign,
                  full::Y,
                  full::Eta,
                  full::Phi,
                  full::P,
                  full::Pt,
                  full::PtXi,
                  full::PtPi0,
                  full::PtPi1,
                  full::M,
                  full::InvMassXi,
                  full::InvMassLambda,
                  full::InvMassXiPi0,
                  full::InvMassXiPi1,
                  full::Chi2Sv,
                  full::Ct,
                  full::DecayLength,
                  full::DecayLengthNormalised,
                  full::DecayLengthXY,
                  full::DecayLengthXYNormalised,
                  full::Cpa,
                  full::CpaXY,
                  full::CpaXi,
                  full::CpaXYXi,
                  full::CpaLam,
                  full::CpaXYLam,
                  full::ImpactParameterXi,
                  full::ImpactParameterNormalisedXi,
                  full::ImpactParameterPi0,
                  full::ImpactParameterNormalisedPi0,
                  full::ImpactParameterPi1,
                  full::ImpactParameterNormalisedPi1,
                  full::MaxNormalisedDeltaIP,
                  // KF specific columns
                  full::Chi2XiVtx,
                  full::Chi2LamVtx,
                  full::Chi2TopoXicPlusToPVBeforeConstraint,
                  full::Chi2TopoXicPlusToPV,
                  //full::Chi2TopoXiToXicPlusBeforeConstraint,
                  //full::Chi2TopoXiToXicPlus,
                  full::DcaXYPi0Pi1,
                  full::DcaXYPi0Xi,
                  full::DcaXYPi1Xi,
                  full::DcaPi0Pi1,
                  full::DcaPi0Xi,
                  full::DcaPi1Xi);

DECLARE_SOA_TABLE(HfCandXicToXiPiPiFulls, "AOD", "HFXICXI2PIFULL",
                  hf_cand_xic0_xicp_to_hadronic::xicp::FlagMcMatchRec,
                  hf_cand_xic0_xicp_to_hadronic::xicp::OriginRec,
                  full::CandidateSelFlag,
                  full::Sign,
                  full::Y,
                  full::Eta,
                  full::Phi,
                  full::P,
                  full::Pt,
                  full::PtXi,
                  full::PtPi0,
                  full::PtPi1,
                  full::M,
                  full::InvMassXi,
                  full::InvMassLambda,
                  full::InvMassXiPi0,
                  full::InvMassXiPi1,
                  full::Chi2Sv,
                  full::Ct,
                  full::DecayLength,
                  full::DecayLengthNormalised,
                  full::DecayLengthXY,
                  full::DecayLengthXYNormalised,
                  full::Cpa,
                  full::CpaXY,
                  full::CpaXi,
                  full::CpaXYXi,
                  full::CpaLam,
                  full::CpaXYLam,
                  full::ImpactParameterXi,
                  full::ImpactParameterNormalisedXi,
                  full::ImpactParameterPi0,
                  full::ImpactParameterNormalisedPi0,
                  full::ImpactParameterPi1,
                  full::ImpactParameterNormalisedPi1,
                  full::MaxNormalisedDeltaIP,
                  // additional columns only stored in the full candidate table
                  full::CpaLamToXi,
                  full::CpaXYLamToXi,
                  full::PPi0,
                  full::PPi1,
                  full::PBachelorPi,
                  full::PPiFromLambda,
                  full::PPrFromLambda,
                  hf_cand_xic0_xicp_to_hadronic::xicp::DcaXiDaughters,
                  hf_cand_xic0_xicp_to_hadronic::xicp::DcaV0Daughters,
                  hf_cand_xic0_xicp_to_hadronic::xicp::DcaPosToPV,
                  hf_cand_xic0_xicp_to_hadronic::xicp::DcaNegToPV,
                  hf_cand_xic0_xicp_to_hadronic::xicp::DcaBachelorToPV,
                  hf_cand_xic0_xicp_to_hadronic::xicp::DcaXYCascToPV,
                  hf_cand_xic0_xicp_to_hadronic::xicp::DcaZCascToPV,
                  hf_cand_xic0_xicp_to_hadronic::xicp::NSigTpcPiFromXicPlus0,
                  hf_cand_xic0_xicp_to_hadronic::xicp::NSigTpcPiFromXicPlus1,
                  hf_cand_xic0_xicp_to_hadronic::xicp::NSigTpcBachelorPi,
                  hf_cand_xic0_xicp_to_hadronic::xicp::NSigTpcPiFromLambda,
                  hf_cand_xic0_xicp_to_hadronic::xicp::NSigTpcPrFromLambda,
                  hf_cand_xic0_xicp_to_hadronic::xicp::NSigTofPiFromXicPlus0,
                  hf_cand_xic0_xicp_to_hadronic::xicp::NSigTofPiFromXicPlus1,
                  hf_cand_xic0_xicp_to_hadronic::xicp::NSigTofBachelorPi,
                  hf_cand_xic0_xicp_to_hadronic::xicp::NSigTofPiFromLambda,
                  hf_cand_xic0_xicp_to_hadronic::xicp::NSigTofPrFromLambda);

DECLARE_SOA_TABLE(HfCandXicToXiPiPiFullKfs, "AOD", "HFXICXI2PIFULKF",
                  hf_cand_xic0_xicp_to_hadronic::xicp::FlagMcMatchRec,
                  hf_cand_xic0_xicp_to_hadronic::xicp::OriginRec,
                  full::CandidateSelFlag,
                  full::Sign,
                  full::Y,
                  full::Eta,
                  full::Phi,
                  full::P,
                  full::Pt,
                  full::PtXi,
                  full::PtPi0,
                  full::PtPi1,
                  full::M,
                  full::InvMassXi,
                  full::InvMassLambda,
                  full::InvMassXiPi0,
                  full::InvMassXiPi1,
                  full::Chi2Sv,
                  full::Ct,
                  full::DecayLength,
                  full::DecayLengthNormalised,
                  full::DecayLengthXY,
                  full::DecayLengthXYNormalised,
                  full::Cpa,
                  full::CpaXY,
                  full::CpaXi,
                  full::CpaXYXi,
                  full::CpaLam,
                  full::CpaXYLam,
                  full::ImpactParameterXi,
                  full::ImpactParameterNormalisedXi,
                  full::ImpactParameterPi0,
                  full::ImpactParameterNormalisedPi0,
                  full::ImpactParameterPi1,
                  full::ImpactParameterNormalisedPi1,
                  full::MaxNormalisedDeltaIP,
                  // additional columns only stored in the full candidate table
                  full::CpaLamToXi,
                  full::CpaXYLamToXi,
                  full::PPi0,
                  full::PPi1,
                  full::PBachelorPi,
                  full::PPiFromLambda,
                  full::PPrFromLambda,
                  hf_cand_xic0_xicp_to_hadronic::xicp::DcaXiDaughters,
                  hf_cand_xic0_xicp_to_hadronic::xicp::DcaV0Daughters,
                  hf_cand_xic0_xicp_to_hadronic::xicp::DcaPosToPV,
                  hf_cand_xic0_xicp_to_hadronic::xicp::DcaNegToPV,
                  hf_cand_xic0_xicp_to_hadronic::xicp::DcaBachelorToPV,
                  hf_cand_xic0_xicp_to_hadronic::xicp::DcaXYCascToPV,
                  hf_cand_xic0_xicp_to_hadronic::xicp::DcaZCascToPV,
                  hf_cand_xic0_xicp_to_hadronic::xicp::NSigTpcPiFromXicPlus0,
                  hf_cand_xic0_xicp_to_hadronic::xicp::NSigTpcPiFromXicPlus1,
                  hf_cand_xic0_xicp_to_hadronic::xicp::NSigTpcBachelorPi,
                  hf_cand_xic0_xicp_to_hadronic::xicp::NSigTpcPiFromLambda,
                  hf_cand_xic0_xicp_to_hadronic::xicp::NSigTpcPrFromLambda,
                  hf_cand_xic0_xicp_to_hadronic::xicp::NSigTofPiFromXicPlus0,
                  hf_cand_xic0_xicp_to_hadronic::xicp::NSigTofPiFromXicPlus1,
                  hf_cand_xic0_xicp_to_hadronic::xicp::NSigTofBachelorPi,
                  hf_cand_xic0_xicp_to_hadronic::xicp::NSigTofPiFromLambda,
                  hf_cand_xic0_xicp_to_hadronic::xicp::NSigTofPrFromLambda,
                  // KF-specific columns
                  full::Chi2XiVtx,
                  full::Chi2LamVtx,
                  full::Chi2TopoXicPlusToPVBeforeConstraint,
                  full::Chi2TopoXicPlusToPV,
                  //full::Chi2TopoXiToXicPlusBeforeConstraint,
                  //full::Chi2TopoXiToXicPlus,
                  full::DcaXYPi0Pi1,
                  full::DcaXYPi0Xi,
                  full::DcaXYPi1Xi,
                  full::DcaPi0Pi1,
                  full::DcaPi0Xi,
                  full::DcaPi1Xi);

DECLARE_SOA_TABLE(HfCandXicToXiPiPiFullPs, "AOD", "HFXICXI2PIFULLP",
                  hf_cand_xic0_xicp_to_hadronic::xicp::FlagMcMatchGen,
                  hf_cand_xic0_xicp_to_hadronic::xicp::OriginGen,
                  hf_cand::PdgBhadMotherPart,
                  full::Pt,
                  full::Eta,
                  full::Phi,
                  full::Y);

DECLARE_SOA_TABLE(HfCandXicToXiPiPiResiduals, "AOD", "HFXICXI2PIRESID",
                  hf_cand_xic0_xicp_to_hadronic::xicp::OriginGen,
                  full::PResidual,
                  full::PtResidual,
                  full::XPvResidual,
                  full::YPvResidual,
                  full::ZPvResidual,
                  full::XPvPull,
                  full::YPvPull,
                  full::ZPvPull,
                  full::XSvResidual,
                  full::YSvResidual,
                  full::ZSvResidual,
                  full::XSvPull,
                  full::YSvPull,
                  full::ZSvPull);
}// o2::aod
 
struct HfTreeCreatorXic0XicpToHadronic {

	struct : ProducesGroup {
		Produces<aod::HfCandXic0Lites> rowCandXic0Lites;
		Produces<aod::HfCandXic0LiteKfs> rowCandXic0LiteKfs;
		Produces<aod::HfCandXic0Fulls> rowCandXic0Fulls;
		Produces<aod::HfCandXic0FullKfs> rowCandXic0FullKfs;
		Produces<aod::HfCandXic0ToXiPiFullPs> rowCandXic0FullParticles;
		Produces<aod::HfCandXic0ToXiPiResiduals> rowCandResiduals;
			        
		Produces<o2::aod::HfCandXicToXiPiPiLites> rowCandXicpLites;
		Produces<o2::aod::HfCandXicToXiPiPiLiteKfs> rowCandXicpLiteKfs;
		Produces<o2::aod::HfCandXicToXiPiPiFulls> rowCandXicpFulls;
		Produces<o2::aod::HfCandXicToXiPiPiFullKfs> rowCandXicpFullKfs;
		Produces<o2::aod::HfCandXicToXiPiPiFullPs> rowCandXicpFullParticles;
		Produces<o2::aod::HfCandXicToXiPiPiResiduals> rowCandXicpResiduals;
	} cursors;
	
	struct : ConfigurableGroup {
		Configurable<int> selectionFlagXic0{"selectionFlagXic0", 1, "Selection flag for Xic0"};
		Configurable<int> selectionFlagXicp{"selectionFlagXicp", 1, "Selection flag for Xicp"};
		Configurable<bool> fillCandidateLite{"fillCandidateLite", true, "Switch to fill lite table with candidate properties"};
		Configurable<bool> fillGenParticleTable{"fillGenParticleTable", false, "Switch to fill table with MC truth for generated particles"};
		Configurable<bool> fillOnlySignal{"fillOnlySignal", false, "Flag to fill derived tables with signal for ML trainings"};
		Configurable<bool> fillOnlyBackground{"fillOnlyBackground", false, "Flag to fill derived tables with background for ML trainings"};
		Configurable<float> downSampleBkgFactor{"downSampleBkgFactor", 1., "Fraction of background candidates to keep for ML trainings"};// ?
		Configurable<float> ptMaxForDownSample{"ptMaxForDownSample", 10., "Maximim pt for the application of the downsampling factor"}; // ?
		// Parameter for MC matching in residual process function
		Configurable<bool> matchDecayedPions{"matchDecayedPions", true, "Match also candidates with daughter pion tracks that decay with kinked topology"};
		Configurable<bool> matchInteractionsWithMaterials{"matchInteractionsWithMaterials", true, "Match also candidates with daughter tracks that interact with material"};
	} configs;

	// Table alises for Xic0
	using SelectedXic0Candidates = soa::Filtered<soa::Join<aod::HfCandXic0, aod::HfSelXic0ToXiPi>>;
	using SelectedXic0CandidatesKf = soa::Filtered<soa::Join<aod::HfCandXic0, aod::HfCandXic0KF, aod::HfSelXic0ToXiPi>>;
	using SelectedXic0CandidatesMc = soa::Filtered<soa::Join<aod::HfCandXic0, aod::HfCandXic0McRec, aod::HfSelXic0ToXiPi>>;
	using SelectedXic0CandidatesKfMc = soa::Filtered<soa::Join<aod::HfCandXic0, aod::HfCandXic0KF, aod::HfCandXic0McRec, aod::HfSelXic0ToXiPi>>;

	Filter filterSelectXic0Candidates = aod::hf_sel_xic0_xicp_to_hadronic::xic0::isSelXic0ToXiPi >= configs.selectionFlagXic0;

	Partition<SelectedXic0CandidatesMc> recXic0Sig = nabs(aod::hf_cand_xic0_xicp_to_hadronic::xic0::flagMcMatchRec) != int8_t(0);
	Partition<SelectedXic0CandidatesMc> recXic0Bkg = nabs(aod::hf_cand_xic0_xicp_to_hadronic::xic0::flagMcMatchRec) == int8_t(0);
	Partition<SelectedXic0CandidatesKfMc> recXic0SigKf = nabs(aod::hf_cand_xic0_xicp_to_hadronic::xic0::flagMcMatchRec) != int8_t(0);
	Partition<SelectedXic0CandidatesKfMc> recXic0BkgKf = nabs(aod::hf_cand_xic0_xicp_to_hadronic::xic0::flagMcMatchRec) == int8_t(0);
	
	// Table alises for Xicp
	using SelectedXicpCandidates = soa::Filtered<soa::Join<aod::HfCandXic, aod::HfSelXicToXiPiPi>>;
	using SelectedXicpCandidatesKf = soa::Filtered<soa::Join<aod::HfCandXic, aod::HfCandXicKF, aod::HfSelXicToXiPiPi>>;
	using SelectedXicpCandidatesMc = soa::Filtered<soa::Join<aod::HfCandXic, aod::HfCandXicMcRec, aod::HfSelXicToXiPiPi>>;
	using SelectedXicpCandidatesKfMc = soa::Filtered<soa::Join<aod::HfCandXic, aod::HfCandXicKF, aod::HfCandXicMcRec, aod::HfSelXicToXiPiPi>>;
	using MatchedGenXicToXiPiPi = soa::Filtered<soa::Join<aod::McParticles, aod::HfCandXicMcGen>>;

	Filter filterSelecXicptCandidates = aod::hf_sel_xic0_xicp_to_hadronic::xicp::isSelXicpToXiPiPi >= configs.selectionFlagXicp;
	Filter filterGenXicToXiPiPi = (nabs(aod::hf_cand_xic0_xicp_to_hadronic::xicp::flagMcMatchGen) == static_cast<int8_t>(BIT(aod::hf_cand_xic0_xicp_to_hadronic::DecayTypeXicp::XicToXiPiPi)) || 
								   nabs(aod::hf_cand_xic0_xicp_to_hadronic::xicp::flagMcMatchGen) == static_cast<int8_t>(BIT(aod::hf_cand_xic0_xicp_to_hadronic::DecayTypeXicp::XicToXiResPiToXiPiPi)));

	Partition<SelectedXicpCandidatesMc> recXicpSig = nabs(aod::hf_cand_xic0_xicp_to_hadronic::xicp::flagMcMatchRec) != int8_t(0);
	Partition<SelectedXicpCandidatesMc> recXicpBkg = nabs(aod::hf_cand_xic0_xicp_to_hadronic::xicp::flagMcMatchRec) == int8_t(0);
	Partition<SelectedXicpCandidatesKfMc> recXicpSigKf = nabs(aod::hf_cand_xic0_xicp_to_hadronic::xicp::flagMcMatchRec) != int8_t(0);
	Partition<SelectedXicpCandidatesKfMc> recXicpBkgKf = nabs(aod::hf_cand_xic0_xicp_to_hadronic::xicp::flagMcMatchRec) == int8_t(0);

	void init(InitContext const&)
	{
	}

	template <bool doMc, bool doKf, typename T>
	void fillXic0CandidateTable(const T& candidate)
	{
		int8_t flagMc{0}, debugMc{0}, originMc{0};
		if constexpr (doMc) {
			flagMc = candidate.flagMcMatchRec();
			debugMc = candidate.debugMcRec();
			originMc = candidate.originRec();
		}
		if constexpr (!doKf) {
			if (configs.fillCandidateLite) {
				cursors.rowCandXic0Lites( flagMc,
										  debugMc,
										  originMc,
										  candidate.isSelXic0ToXiPi(),
										  candidate.y(o2::constants::physics::MassXiC0),
										  candidate.eta(),
										  candidate.phi(),
										  candidate.p(),
										  candidate.pt(),
										  candidate.ptProng0(),
										  candidate.ptProng1(),
										  candidate.invMassXic0(),
										  candidate.invMassXi(),
										  candidate.invMassLambda(),
										  candidate.chi2PCA(),
										  candidate.ct(o2::constants::physics::MassXiC0),
										  candidate.decayLength(),
										  candidate.decayLengthNormalised(),
										  candidate.decayLengthXY(),
										  candidate.decayLengthXYNormalised(),
										  candidate.cpa(),
										  candidate.cpaXY(),
										  candidate.cosPaXi(),
										  candidate.cosPaXYXi(),
										  candidate.cosPaLambda(),
										  candidate.cosPaXYLambda(),
										  candidate.impactParameter0(),
										  candidate.impactParameterNormalised0(),
										  candidate.impactParameter1(),
										  candidate.impactParameterNormalised1(),
										  candidate.maxNormalisedDeltaIP() );
			} else { // Fill full table
				cursors.rowCandXic0Fulls( flagMc,
										  debugMc,
										  originMc,
										  candidate.isSelXic0ToXiPi(),
										  candidate.y(o2::constants::physics::MassXiC0),
										  candidate.eta(),
										  candidate.phi(),
										  candidate.p(),
										  candidate.pt(),
										  candidate.ptProng0(),
										  candidate.ptProng1(),
										  candidate.invMassXic0(),
										  candidate.invMassXi(),
										  candidate.invMassLambda(),
										  candidate.chi2PCA(),
										  candidate.ct(o2::constants::physics::MassXiC0),
										  candidate.decayLength(),
										  candidate.decayLengthNormalised(),
										  candidate.decayLengthXY(),
										  candidate.decayLengthXYNormalised(),
										  candidate.cpa(),
										  candidate.cpaXY(),
										  candidate.cosPaXi(),
										  candidate.cosPaXYXi(),
										  candidate.cosPaLambda(),
										  candidate.cosPaXYLambda(),
										  candidate.impactParameter0(),
										  candidate.impactParameterNormalised0(),
										  candidate.impactParameter1(),
										  candidate.impactParameterNormalised1(),
										  candidate.maxNormalisedDeltaIP(),
										  // Additional colums only stored in the full candidate table
										  candidate.cosPaLambdaToXi(),
										  candidate.cosPaXYLambdaToXi(),
										  candidate.pProng1(), // Pi <- Xic0
										  candidate.pBachelorPi(), // Pi <- Xi-
										  candidate.pPiFromLambda(), // Pi <- V0 <- Xi-
										  candidate.pPrFromLambda(), // Pr <- V0 <- Xi-
										  candidate.dcaXiDaughters(),
										  candidate.dcaV0Daughters(),
										  candidate.dcaPosToPV(),
										  candidate.dcaNegToPV(),
										  candidate.dcaBachelorToPV(),
										  candidate.dcaXYCascToPV(),
										  candidate.dcaZCascToPV(),
										  candidate.nSigTpcPiFromXic0(),
										  candidate.nSigTpcBachelorPi(),
										  candidate.nSigTpcPiFromLambda(),
										  candidate.nSigTpcPrFromLambda(),
										  candidate.nSigTofPiFromXic0(),
										  candidate.nSigTofBachelorPi(),
										  candidate.nSigTofPiFromLambda(),
										  candidate.nSigTofPrFromLambda() );
			}
		} else {// Fill with KF informations
			if (configs.fillCandidateLite) {
				cursors.rowCandXic0LiteKfs( flagMc,
										  	debugMc,
										  	originMc,
										  	candidate.isSelXic0ToXiPi(),
										  	candidate.y(o2::constants::physics::MassXiC0),
										  	candidate.eta(),
										  	candidate.phi(),
										  	candidate.p(),
										  	candidate.pt(),
											candidate.ptProng0(),
											candidate.ptProng1(),
											candidate.invMassXic0(),
											candidate.invMassXi(),
											candidate.invMassLambda(),
											candidate.chi2PCA(),
											candidate.ct(o2::constants::physics::MassXiC0),
											candidate.decayLength(),
											candidate.decayLengthNormalised(),
											candidate.decayLengthXY(),
											candidate.decayLengthXYNormalised(),
											candidate.cpa(),
											candidate.cpaXY(),
											candidate.cosPaXi(),
											candidate.cosPaXYXi(),
											candidate.cosPaLambda(),
											candidate.cosPaXYLambda(),
											candidate.impactParameter0(),
											candidate.impactParameterNormalised0(),
											candidate.impactParameter1(),
											candidate.impactParameterNormalised1(),
											candidate.maxNormalisedDeltaIP(),
											// KF-specific columns
											candidate.kfCascadeChi2(),
											candidate.kfV0Chi2(),
											candidate.chi2TopoXic0ToPVBeforeConstraint(),
											candidate.chi2TopoXic0ToPV(),
											candidate.chi2TopoXiToXic0BeforeConstraint(),
											candidate.chi2TopoXiToXic0(),
											candidate.dcaXYPiXi(),
											candidate.dcaPiXi() );

			} else { // Fill full table with KF
				cursors.rowCandXic0FullKfs( flagMc,
										  	debugMc,
										  	originMc,
										  	candidate.isSelXic0ToXiPi(),
										  	candidate.y(o2::constants::physics::MassXiC0),
										  	candidate.eta(),
										  	candidate.phi(),
										  	candidate.p(),
										  	candidate.pt(),
											candidate.ptProng0(),
											candidate.ptProng1(),
											candidate.invMassXic0(),
											candidate.invMassXi(),
											candidate.invMassLambda(),
											candidate.chi2PCA(),
											candidate.ct(o2::constants::physics::MassXiC0),
											candidate.decayLength(),
											candidate.decayLengthNormalised(),
											candidate.decayLengthXY(),
											candidate.decayLengthXYNormalised(),
											candidate.cpa(),
											candidate.cpaXY(),
											candidate.cosPaXi(),
											candidate.cosPaXYXi(),
											candidate.cosPaLambda(),
											candidate.cosPaXYLambda(),
											candidate.impactParameter0(),
											candidate.impactParameterNormalised0(),
											candidate.impactParameter1(),
											candidate.impactParameterNormalised1(),
											candidate.maxNormalisedDeltaIP(),
										  	// Additional colums only stored in the full candidate table
										  	candidate.cosPaLambdaToXi(),
										  	candidate.cosPaXYLambdaToXi(),
										  	candidate.pProng1(), // Pi <- Xic0
										  	candidate.pBachelorPi(), // Pi <- Xi-
										  	candidate.pPiFromLambda(), // Pi <- V0 <- Xi-
										  	candidate.pPrFromLambda(), // Pr <- V0 <- Xi-
										  	candidate.dcaXiDaughters(),
										  	candidate.dcaV0Daughters(),
										  	candidate.dcaPosToPV(),
										  	candidate.dcaNegToPV(),
										  	candidate.dcaBachelorToPV(),
										  	candidate.dcaXYCascToPV(),
										  	candidate.dcaZCascToPV(),
										  	candidate.nSigTpcPiFromXic0(),
										  	candidate.nSigTpcBachelorPi(),
										  	candidate.nSigTpcPiFromLambda(),
										  	candidate.nSigTpcPrFromLambda(),
										  	candidate.nSigTofPiFromXic0(),
										  	candidate.nSigTofBachelorPi(),
										  	candidate.nSigTofPiFromLambda(),
										  	candidate.nSigTofPrFromLambda(), 
											// KF-specific columns
											candidate.kfCascadeChi2(),
											candidate.kfV0Chi2(),
											candidate.chi2TopoXic0ToPVBeforeConstraint(),
											candidate.chi2TopoXic0ToPV(),
											candidate.chi2TopoXiToXic0BeforeConstraint(),
											candidate.chi2TopoXiToXic0(),
											candidate.dcaXYPiXi(),
											candidate.dcaPiXi() );
			}
		}
	} // end fillXic0CandidateTable

	template <bool doMc, bool doKf, typename T2>
	void fillXicpCandidateTable(const T2& candidate)
	{
		int8_t flagMc = 0;
		int8_t originMc = 0;
		if constexpr (doMc) {
			flagMc = candidate.flagMcMatchRec();
			originMc = candidate.originRec();
		}
		if constexpr (!doKf) {
			if (configs.fillCandidateLite) {
				cursors.rowCandXicpLites( flagMc,
										  originMc,
										  candidate.isSelXicToXiPiPi(),
										  candidate.sign(),
										  candidate.y(o2::constants::physics::MassXiCPlus),
										  candidate.eta(),
										  candidate.phi(),
										  candidate.p(),
										  candidate.pt(),
										  candidate.ptProng0(),
										  candidate.ptProng1(),
										  candidate.ptProng2(),
										  candidate.invMassXicPlus(),
										  candidate.invMassXi(),
										  candidate.invMassLambda(),
										  candidate.invMassXiPi0(),
										  candidate.invMassXiPi1(),
										  candidate.chi2PCA(),
										  candidate.ct(o2::constants::physics::MassXiCPlus),
										  candidate.decayLength(),
										  candidate.decayLengthNormalised(),
										  candidate.decayLengthXY(),
										  candidate.decayLengthXYNormalised(),
										  candidate.cpa(), 
										  candidate.cpaXY(),
										  candidate.cpaXi(), //candidate.cosPaXi(),
										  candidate.cpaXYXi(), //candidate.cosPaXYXi(),
										  candidate.cpaLambda(), //candidate.cosPaLambda(),
										  candidate.cpaXYLambda(), //candidate.cosPaXYLambda(),
										  candidate.impactParameter0(),
										  candidate.impactParameterNormalised0(),
										  candidate.impactParameter1(),
										  candidate.impactParameterNormalised1(),
										  candidate.impactParameter2(),
										  candidate.impactParameterNormalised2(),
										  candidate.maxNormalisedDeltaIP() );
      
			} else {
				cursors.rowCandXicpFulls( flagMc,
										  originMc,
										  candidate.isSelXicToXiPiPi(),
										  candidate.sign(),
										  candidate.y(o2::constants::physics::MassXiCPlus),
										  candidate.eta(),
										  candidate.phi(), 
										  candidate.p(),
										  candidate.pt(),
										  candidate.ptProng0(),
										  candidate.ptProng1(),
										  candidate.ptProng2(),
										  candidate.invMassXicPlus(),
										  candidate.invMassXi(),
										  candidate.invMassLambda(), 
										  candidate.invMassXiPi0(),
										  candidate.invMassXiPi1(),
										  candidate.chi2PCA(),
										  candidate.ct(o2::constants::physics::MassXiCPlus),
										  candidate.decayLength(),
										  candidate.decayLengthNormalised(),
										  candidate.decayLengthXY(),
										  candidate.decayLengthXYNormalised(),
										  candidate.cpa(), 
										  candidate.cpaXY(),
										  candidate.cpaXi(), //candidate.cosPaXi(),
										  candidate.cpaXYXi(), //candidate.cosPaXYXi(),
										  candidate.cpaLambda(), //candidate.cosPaLambda(), 
										  candidate.cpaXYLambda(), //candidate.cosPaXYLambda(),
										  candidate.impactParameter0(), 
										  candidate.impactParameterNormalised0(),
										  candidate.impactParameter1(),
										  candidate.impactParameterNormalised1(),
										  candidate.impactParameter2(), 
										  candidate.impactParameterNormalised2(),
										  candidate.maxNormalisedDeltaIP(),
										  // additional columns only stored in the full candidate table
										  candidate.cpaLambdaToXi(), //candidate.cosPaLambdaToXi(),
										  candidate.cpaXYLambdaToXi(), //candidate.cosPaXYLambdaToXi(), 
										  candidate.pProng1(),
										  candidate.pProng2(),
										  candidate.pBachelorPi(),
										  candidate.pPiFromLambda(),
										  candidate.pPrFromLambda(), 
										  candidate.dcaXiDaughters(),
										  candidate.dcaV0Daughters(),
										  candidate.dcaPosToPV(),
										  candidate.dcaNegToPV(),
										  candidate.dcaBachelorToPV(),
										  candidate.dcaXYCascToPV(),
										  candidate.dcaZCascToPV(),
										  candidate.nSigTpcPiFromXicPlus0(),
										  candidate.nSigTpcPiFromXicPlus1(),
										  candidate.nSigTpcBachelorPi(),
										  candidate.nSigTpcPiFromLambda(),
										  candidate.nSigTpcPrFromLambda(),
										  candidate.nSigTofPiFromXicPlus0(),
										  candidate.nSigTofPiFromXicPlus1(),
										  candidate.nSigTofBachelorPi(),
										  candidate.nSigTofPiFromLambda(),
										  candidate.nSigTofPrFromLambda() );
			}
    
		} else {
			if (configs.fillCandidateLite) {
				cursors.rowCandXicpLiteKfs( flagMc,
											originMc,
											candidate.isSelXicToXiPiPi(),
											candidate.sign(),
											candidate.y(o2::constants::physics::MassXiCPlus),
											candidate.eta(),
											candidate.phi(),
											candidate.p(),
											candidate.pt(),
											candidate.ptProng0(),
											candidate.ptProng1(),
											candidate.ptProng2(),
											candidate.invMassXicPlus(), 
											candidate.invMassXi(),
											candidate.invMassLambda(),
											candidate.invMassXiPi0(),
											candidate.invMassXiPi1(),
											candidate.chi2PCA(),
											candidate.ct(o2::constants::physics::MassXiCPlus),
											candidate.kfDecayLength(),
											candidate.kfDecayLengthNormalised(),
											candidate.kfDecayLengthXY(),
											candidate.kfDecayLengthXYNormalised(),
											candidate.cpa(),
											candidate.cpaXY(),
											candidate.cpaXi(), //candidate.cosPaXi(),
											candidate.cpaXYXi(), //candidate.cosPaXYXi(),
											candidate.cpaLambda(), //candidate.cosPaLambda(),
											candidate.cpaXYLambda(), //candidate.cosPaXYLambda(),
											candidate.impactParameter0(),
											candidate.impactParameterNormalised0(),
											candidate.impactParameter1(),
											candidate.impactParameterNormalised1(),
											candidate.impactParameter2(),
											candidate.impactParameterNormalised2(),
											candidate.maxNormalisedDeltaIP(),
											// KF-specific columns
											candidate.kfCascadeChi2(),
											candidate.kfV0Chi2(),
											candidate.chi2TopoXicPlusToPVBefConst(), //candidate.chi2TopoXicPlusToPVBeforeConstraint(),
											candidate.chi2TopoXicPlusToPV(),
											//candidate.chi2TopoXiToXicPlusBeforeConstraint(),
											//candidate.chi2TopoXiToXicPlus(),
											candidate.dcaXYPi0Pi1(),
											candidate.dcaXYPi0Xi(),
											candidate.dcaXYPi1Xi(),
											candidate.dcaPi0Pi1(),
											candidate.dcaPi0Xi(),
											candidate.dcaPi1Xi());
      
			} else {
				cursors.rowCandXicpFullKfs( flagMc,
											originMc,
											candidate.isSelXicToXiPiPi(),
											candidate.sign(),
											candidate.y(o2::constants::physics::MassXiCPlus),
											candidate.eta(),
											candidate.phi(),
											candidate.p(),
											candidate.pt(),
											candidate.ptProng0(),
											candidate.ptProng1(),
											candidate.ptProng2(),
											candidate.invMassXicPlus(),
											candidate.invMassXi(),
											candidate.invMassLambda(),
											candidate.invMassXiPi0(),
											candidate.invMassXiPi1(),
											candidate.chi2PCA(),
											candidate.ct(o2::constants::physics::MassXiCPlus),
											candidate.kfDecayLength(),
											candidate.kfDecayLengthNormalised(),
											candidate.kfDecayLengthXY(),
											candidate.kfDecayLengthXYNormalised(),
											candidate.cpa(),
											candidate.cpaXY(),
											candidate.cpaXi(), //candidate.cosPaXi(),
											candidate.cpaXYXi(), //candidate.cosPaXYXi(),
											candidate.cpaLambda(), //candidate.cosPaLambda(),
											candidate.cpaXYLambda(), //candidate.cosPaXYLambda(),
											candidate.impactParameter0(),
											candidate.impactParameterNormalised0(),
											candidate.impactParameter1(),
											candidate.impactParameterNormalised1(),
											candidate.impactParameter2(),
											candidate.impactParameterNormalised2(),
											candidate.maxNormalisedDeltaIP(),
											// additional columns only stored in the full candidate table
											candidate.cpaLambdaToXi(), //candidate.cosPaLambdaToXi(),
											candidate.cpaXYLambdaToXi(), //candidate.cosPaXYLambdaToXi(),
											candidate.pProng1(),
											candidate.pProng2(),
											candidate.pBachelorPi(),
											candidate.pPiFromLambda(),
											candidate.pPrFromLambda(),
											candidate.dcaXiDaughters(),
											candidate.dcaV0Daughters(),
											candidate.dcaPosToPV(), 
											candidate.dcaNegToPV(),
											candidate.dcaBachelorToPV(),
											candidate.dcaXYCascToPV(),
											candidate.dcaZCascToPV(),
											candidate.nSigTpcPiFromXicPlus0(),
											candidate.nSigTpcPiFromXicPlus1(),
											candidate.nSigTpcBachelorPi(),
											candidate.nSigTpcPiFromLambda(),
											candidate.nSigTpcPrFromLambda(),
											candidate.nSigTofPiFromXicPlus0(),
											candidate.nSigTofPiFromXicPlus1(),
											candidate.nSigTofBachelorPi(),
											candidate.nSigTofPiFromLambda(),
											candidate.nSigTofPrFromLambda(),
											// KF-specific columns
											candidate.kfCascadeChi2(),
											candidate.kfV0Chi2(),
											candidate.chi2TopoXicPlusToPVBefConst(), //candidate.chi2TopoXicPlusToPVBeforeConstraint(),
											candidate.chi2TopoXicPlusToPV(),
											//candidate.chi2TopoXiToXicPlusBeforeConstraint(),
											//candidate.chi2TopoXiToXicPlus(),
											candidate.dcaXYPi0Pi1(),
											candidate.dcaXYPi0Xi(),
											candidate.dcaXYPi1Xi(),
											candidate.dcaPi0Pi1(),
											candidate.dcaPi0Xi(),
											candidate.dcaPi1Xi());
			}
		}
	} // end of xicp candidate fill function


	////////////////////////////////////////////////////
	////											////
	////			Make TTree with Data			////
	////											////
	////////////////////////////////////////////////////


	void processDataXic0(SelectedXic0Candidates const& candidates) // -> To be filled
	{
		if (configs.fillCandidateLite) {
			cursors.rowCandXic0Lites.reserve(candidates.size());
		} else {
			cursors.rowCandXic0Fulls.reserve(candidates.size());
		}

		for (const auto& candidate : candidates) {
			if (configs.fillOnlyBackground && configs.downSampleBkgFactor <1.) {
				float pseudoRndm = candidate.ptProng1()*1000. - static_cast<int64_t>(candidate.ptProng1()*1000);
				if (pseudoRndm >= configs.downSampleBkgFactor && candidate.pt() < configs.ptMaxForDownSample) {
					continue;
				}
			}
			fillXic0CandidateTable<false, false>(candidate);
		}
	}
	PROCESS_SWITCH(HfTreeCreatorXic0XicpToHadronic, processDataXic0, "Process Data with DCAFitter reconstruction", true);
	

	void processDataXic0Kf(SelectedXic0CandidatesKf const& candidates)
	{
		// Reserve memory 	
		if (configs.fillCandidateLite) {
			cursors.rowCandXic0LiteKfs.reserve(candidates.size());
		} else {
			cursors.rowCandXic0FullKfs.reserve(candidates.size());
		}

		// Fill candidate properties
		for (auto const& candidate : candidates) {
			if (configs.fillOnlyBackground && configs.downSampleBkgFactor < 1.) {
				float pseudoRndm = candidate.ptProng1()*1000. - static_cast<int64_t>(candidate.ptProng1()*1000.);
				if (pseudoRndm >= configs.downSampleBkgFactor && candidate.pt() < configs.ptMaxForDownSample) {
					continue;
				}
			}
			fillXic0CandidateTable<false, true>(candidate);
		}
	}
	PROCESS_SWITCH(HfTreeCreatorXic0XicpToHadronic, processDataXic0Kf, "Process Data with KFParticle reconstruction", false);

	void processDataXicp(SelectedXicpCandidates const& candidates)
	{
		// Filling candidate properties
		if (configs.fillCandidateLite) {
			cursors.rowCandXicpLites.reserve(candidates.size());
		} else {
			cursors.rowCandXicpFulls.reserve(candidates.size());
		}
		for (const auto& candidate : candidates) {
			if (configs.fillOnlyBackground && configs.downSampleBkgFactor < 1.) {
				float pseudoRndm = candidate.ptProng1() * 1000. - static_cast<int64_t>(candidate.ptProng1() * 1000);
				if (pseudoRndm >= configs.downSampleBkgFactor && candidate.pt() < configs.ptMaxForDownSample) {
					continue;
				}
			}
			fillXicpCandidateTable<false, false>(candidate);
		}
	}
	PROCESS_SWITCH(HfTreeCreatorXic0XicpToHadronic, processDataXicp, "Process Xicp data with DCAFitter reconstruction", false);

	void processDataXicpKf(SelectedXicpCandidatesKf const& candidates)
	{
		// Filling candidate properties
		if (configs.fillCandidateLite) {
			cursors.rowCandXicpLites.reserve(candidates.size());
		} else {
			cursors.rowCandXicpFulls.reserve(candidates.size());
		}
		for (const auto& candidate : candidates) {
			if (configs.fillOnlyBackground && configs.downSampleBkgFactor < 1.) {
				float pseudoRndm = candidate.ptProng1() * 1000. - static_cast<int64_t>(candidate.ptProng1() * 1000);
				if (pseudoRndm >= configs.downSampleBkgFactor && candidate.pt() < configs.ptMaxForDownSample) {
					continue;
				}
			}
			fillXicpCandidateTable<false, true>(candidate);
		}
	}
	PROCESS_SWITCH(HfTreeCreatorXic0XicpToHadronic, processDataXicpKf, "Process Xicp data with KFParticle reconstruction", false);

	////////////////////////////////////////////////////
	////											////
	////			Make TTree with MC				////
	////											////
	////////////////////////////////////////////////////


	void processMcXic0(SelectedXic0CandidatesMc const& candidates,
				   soa::Join<aod::McParticles, aod::HfCandXic0McGen> const& particles) // -> To be filled
	{
		//! Fill only signal for ML Training
		//!! In the reference code, the cursors for DCAFitter was used instead....
		if (configs.fillOnlySignal) {

			if (configs.fillCandidateLite) {
				cursors.rowCandXic0Lites.reserve(recXic0Sig.size());
			} else {
				cursors.rowCandXic0Fulls.reserve(recXic0Sig.size());
			}

			for (const auto& candidate : recXic0Sig) {
				fillXic0CandidateTable<true, false>(candidate);
			}

		//! Fill only Bkg for ML Training
		} else if (configs.fillOnlyBackground) {

			if (configs.fillCandidateLite) {
				cursors.rowCandXic0Lites.reserve(recXic0Bkg.size());
			} else {
				cursors.rowCandXic0Fulls.reserve(recXic0Bkg.size());
			}

			for (const auto& candidate : recXic0Bkg) {
				float pseudoRndm = candidate.ptProng1()*1000. - static_cast<int64_t>(candidate.ptProng1()*1000);
				if (candidate.pt() < configs.ptMaxForDownSample && pseudoRndm >= configs.downSampleBkgFactor) {
					continue;
				}
				fillXic0CandidateTable<true, false>(candidate);
			}

		//! Fill whole reconstructed(?) candidates
		} else {

			if (configs.fillCandidateLite) {
				cursors.rowCandXic0Lites.reserve(candidates.size());
			} else {
				cursors.rowCandXic0Fulls.reserve(candidates.size());
			}

			for (auto const& candidate : candidates) {
				fillXic0CandidateTable<true, false>(candidate);
			}

		}
	
		if (configs.fillGenParticleTable) {
			cursors.rowCandXic0FullParticles.reserve(particles.size());
			for (auto const& particle : particles) {
				if (TESTBIT(std::abs(particle.flagMcMatchGen()), aod::hf_cand_xic0_xicp_to_hadronic::Xic0ToXiPi)) {

					cursors.rowCandXic0FullParticles( particle.flagMcMatchGen(),
													  particle.debugMcGen(),
													  particle.originGen(),
													  particle.pt(),
													  particle.eta(),
													  particle.phi(),
													  RecoDecay::y(particle.pVector(), o2::constants::physics::MassXiC0) );
				}
			}// Gen particle loop
		}
	}
	PROCESS_SWITCH(HfTreeCreatorXic0XicpToHadronic, processMcXic0, "Process MC with DCAFitter reconstruction", false);


	void processMcXic0Kf(SelectedXic0CandidatesKfMc const& candidates,
					 soa::Join<aod::McParticles, aod::HfCandXic0McGen> const& particles)
	{
		//! Fill only signal for ML Training
		//!! In the reference code, the cursors for DCAFitter was used instead....
		if (configs.fillOnlySignal) {

			if (configs.fillCandidateLite) {
				cursors.rowCandXic0LiteKfs.reserve(recXic0SigKf.size());
			} else {
				cursors.rowCandXic0FullKfs.reserve(recXic0SigKf.size());
			}

			for (const auto& candidate : recXic0SigKf) {
				fillXic0CandidateTable<true, true>(candidate);
			}

		//! Fill only Bkg for ML Training
		} else if (configs.fillOnlyBackground) {

			if (configs.fillCandidateLite) {
				cursors.rowCandXic0LiteKfs.reserve(recXic0BkgKf.size());
			} else {
				cursors.rowCandXic0FullKfs.reserve(recXic0BkgKf.size());
			}

			for (const auto& candidate : recXic0BkgKf) {
				float pseudoRndm = candidate.ptProng1()*1000. - static_cast<int64_t>(candidate.ptProng1()*1000);
				if (candidate.pt() < configs.ptMaxForDownSample && pseudoRndm >= configs.downSampleBkgFactor) {
					continue;
				}
				fillXic0CandidateTable<true, true>(candidate);
			}

		//! Fill whole reconstructed(?) candidates
		} else {

			if (configs.fillCandidateLite) {
				cursors.rowCandXic0LiteKfs.reserve(candidates.size());
			} else {
				cursors.rowCandXic0FullKfs.reserve(candidates.size());
			}

			for (auto const& candidate : candidates) {
				fillXic0CandidateTable<true, true>(candidate);
			}

		}
	
		if (configs.fillGenParticleTable) {
			cursors.rowCandXic0FullParticles.reserve(particles.size());
			for (auto const& particle : particles) {
				if (TESTBIT(std::abs(particle.flagMcMatchGen()), aod::hf_cand_xic0_xicp_to_hadronic::Xic0ToXiPi)) {

					cursors.rowCandXic0FullParticles( particle.flagMcMatchGen(),
													  particle.debugMcGen(),
													  particle.originGen(),
													  particle.pt(),
													  particle.eta(),
													  particle.phi(),
													  RecoDecay::y(particle.pVector(), o2::constants::physics::MassXiC0) );
				}
			}// Gen particle loop
		}
	}
	PROCESS_SWITCH(HfTreeCreatorXic0XicpToHadronic, processMcXic0Kf, "Process MC with KFParticle reconstruction", false);

	void processXic0Residuals(SelectedXic0CandidatesMc const&,
						  aod::TracksWMc const&  tracks,
						  aod::McParticles const& particles,
						  aod::McCollisions const&)
	{
		cursors.rowCandResiduals.reserve(recXic0Sig.size());

		recXic0Sig->bindExternalIndices(&tracks); // -> What does this feature do exactly?
											  
		std::vector<int> arrDaughIdx;
		int idxRecXic0;
		int8_t sign;
		int8_t origin;
		std::array<float, 3> pvResiduals;
		std::array<float, 3> svResiduals;
		std::array<float, 3> pvPulls;
		std::array<float, 3> svPulls;
		
		for (auto const& candidate : recXic0Sig) {

			// Initialization
			arrDaughIdx.clear();
			idxRecXic0 = -1;
			origin = 0;
			pvResiduals = {-9999.9};
			svResiduals = {-9999.9};
			pvPulls = {-9999.9};
			svPulls = {-9999.9};

			auto arrDaughters = std::array{ candidate.pi_as<aod::TracksWMc>(),				// pi <- Xic0
											candidate.bachelor_as<aod::TracksWMc>(),		// pi <- cascade
											candidate.posTrack_as<aod::TracksWMc>(),		// p <- V0
											candidate.negTrack_as<aod::TracksWMc>() };		// pi <- V0
			// Get Xic0 as MC particle
			idxRecXic0 = RecoDecay::getMatchedMCRec(particles, arrDaughters, Pdg::kXiC0, std::array{+kPiPlus, +kPiMinus, +kProton, +kPiMinus}, true, &sign, 3);
			if (idxRecXic0 == -1) {
				continue;
			}
			auto particleXic0Gen = particles.rawIteratorAt(idxRecXic0);
			origin = RecoDecay::getCharmHadronOrigin(particles, particleXic0Gen, true); // -> What does this do?

			// Get Mc collision
			auto mcCollision = particleXic0Gen.mcCollision_as<aod::McCollisions>();

			// Get Xic0 daughters as MC Particle
			RecoDecay::getDaughters(particleXic0Gen, &arrDaughIdx, std::array{+kXiMinus, +kPiPlus}, 1); // -> MaxDepth set to 1...correct?
			auto daughXic0 = particles.rawIteratorAt(arrDaughIdx[0]);

			// Calculate residuals and pulls
			float pResidual = candidate.p() - particleXic0Gen.p();
			float ptResidual = candidate.pt() - particleXic0Gen.pt();
			pvResiduals[0] = candidate.posX() - mcCollision.posX();
			pvResiduals[1] = candidate.posY() - mcCollision.posY();
			pvResiduals[2] = candidate.posZ() - mcCollision.posZ();
			svResiduals[0] = candidate.xSecondaryVertex() - daughXic0.vx();
			svResiduals[1] = candidate.ySecondaryVertex() - daughXic0.vy();
			svResiduals[2] = candidate.zSecondaryVertex() - daughXic0.vz();
			try {
				pvPulls[0] = pvResiduals[0] / candidate.xPvErr();
				pvPulls[1] = pvResiduals[1] / candidate.yPvErr();
				pvPulls[2] = pvResiduals[2] / candidate.zPvErr();
				svPulls[0] = svResiduals[0] / candidate.xSvErr();
				svPulls[1] = svResiduals[1] / candidate.ySvErr();
				svPulls[2] = svResiduals[2] / candidate.zSvErr();
			} catch (const std::runtime_error& error) {
				LOG(info) << "Run time error found: " << error.what() << ". Set values of vertex pulls to -9999.9";
			}

			cursors.rowCandResiduals( origin,
									  pResidual,
									  ptResidual,
									  pvResiduals[0],
									  pvResiduals[1],
									  pvResiduals[2],
									  pvPulls[0],
									  pvPulls[1],
									  pvPulls[2],
									  svResiduals[0],
									  svResiduals[1],
									  svResiduals[2],
									  svPulls[0],
									  svPulls[1],
									  svPulls[2] );
																		
		}
	}
	PROCESS_SWITCH(HfTreeCreatorXic0XicpToHadronic, processXic0Residuals, "Process residuals and pulls for both DCAFitter and KFParticle reconstruction", false);

	void processMcXicp(SelectedXicpCandidatesMc const& candidates,
                 MatchedGenXicToXiPiPi const& particles)
	{
		// Filling candidate properties
		if (configs.fillOnlySignal) {
			if (configs.fillCandidateLite) {
				cursors.rowCandXicpLites.reserve(recXicpSig.size());
			} else {
				cursors.rowCandXicpFulls.reserve(recXicpSig.size());
			}
			for (const auto& candidate : recXicpSig) {
				fillXicpCandidateTable<true, false>(candidate);
			}
		} else if (configs.fillOnlyBackground) {
			if (configs.fillCandidateLite) {
				cursors.rowCandXicpLites.reserve(recXicpBkg.size());
			} else {
				cursors.rowCandXicpFulls.reserve(recXicpBkg.size());
			}
			for (const auto& candidate : recXicpBkg) {
				float pseudoRndm = candidate.ptProng1() * 1000. - static_cast<int64_t>(candidate.ptProng1() * 1000);
				if (candidate.pt() < configs.ptMaxForDownSample && pseudoRndm >= configs.downSampleBkgFactor) {
					continue;
				}
				fillXicpCandidateTable<true, false>(candidate);
			}
		} else {
			if (configs.fillCandidateLite) {
				cursors.rowCandXicpLites.reserve(candidates.size());
			} else {
				cursors.rowCandXicpFulls.reserve(candidates.size());
			}
			for (const auto& candidate : candidates) {
				fillXicpCandidateTable<true, false>(candidate);
			}
		}

		if (configs.fillGenParticleTable) {
			cursors.rowCandXicpFullParticles.reserve(particles.size());
			for (const auto& particle : particles) {
				cursors.rowCandXicpFullParticles(
						particle.flagMcMatchGen(),
						particle.originGen(),
						particle.pdgBhadMotherPart(),
						particle.pt(),
						particle.eta(),
						particle.phi(),
						RecoDecay::y(particle.pVector(), 
						o2::constants::physics::MassXiCPlus));
			}
		}
	}
	PROCESS_SWITCH(HfTreeCreatorXic0XicpToHadronic, processMcXicp, "Process Xicp MC with DCAFitter reconstruction", false);

	void processMcXicpKf(SelectedXicpCandidatesKfMc const& candidates,
                   MatchedGenXicToXiPiPi const& particles)
	{
		// Filling candidate properties
		if (configs.fillOnlySignal) {
			if (configs.fillCandidateLite) {
				cursors.rowCandXicpLites.reserve(recXicpSigKf.size());
			} 
			else {
				cursors.rowCandXicpFulls.reserve(recXicpSigKf.size());
			}
			for (const auto& candidate : recXicpSigKf) {
				fillXicpCandidateTable<true, true>(candidate);
			}
		} else if (configs.fillOnlyBackground) {
			if (configs.fillCandidateLite) {
				cursors.rowCandXicpLites.reserve(recXicpBkgKf.size());
			} else {
				cursors.rowCandXicpFulls.reserve(recXicpBkgKf.size());
			}
			for (const auto& candidate : recXicpBkgKf) {
				float pseudoRndm = candidate.ptProng1() * 1000. - static_cast<int64_t>(candidate.ptProng1() * 1000);
				if (candidate.pt() < configs.ptMaxForDownSample && pseudoRndm >= configs.downSampleBkgFactor) {
					continue;
				}
				fillXicpCandidateTable<true, true>(candidate);
			}
		} else {
			if (configs.fillCandidateLite) {
				cursors.rowCandXicpLites.reserve(candidates.size());
			} else {
				cursors.rowCandXicpFulls.reserve(candidates.size());
			}
			for (const auto& candidate : candidates) {
				fillXicpCandidateTable<true, true>(candidate);
			}
		}

		if (configs.fillGenParticleTable) {
			cursors.rowCandXicpFullParticles.reserve(particles.size());
			for (const auto& particle : particles) {
				cursors.rowCandXicpFullParticles(
						particle.flagMcMatchGen(),
						particle.originGen(),
						particle.pdgBhadMotherPart(),
						particle.pt(),
						particle.eta(),
						particle.phi(),
						RecoDecay::y(particle.pVector(), 
						o2::constants::physics::MassXiCPlus));
			}
		}
	}
	PROCESS_SWITCH(HfTreeCreatorXic0XicpToHadronic, processMcXicpKf, "Process Xicp MC with KF Particle reconstruction", false);

	void processResidualsXicp(SelectedXicpCandidatesMc const&,
                        aod::TracksWMc const& tracks,
                        aod::McParticles const& particles,
                        aod::McCollisions const&)
	{
		cursors.rowCandXicpResiduals.reserve(recXicpSig.size());

		recXicpSig->bindExternalIndices(&tracks);

		std::vector<int> arrDaughIndex;
		int indexRecXicPlus;
		int origin;
		std::array<float, 3> pvResiduals;
		std::array<float, 3> svResiduals;
		std::array<float, 3> pvPulls;
		std::array<float, 3> svPulls;

		for (const auto& candidate : recXicpSig) {
			arrDaughIndex.clear();
			indexRecXicPlus = -1;
			origin = RecoDecay::OriginType::None;
			pvResiduals = {-9999.9};
			svResiduals = {-9999.9};
			pvPulls = {-9999.9};
			svPulls = {-9999.9};

			auto arrayDaughters = std::array{candidate.pi0_as<aod::TracksWMc>(),       // pi <- Xic
                                       		 candidate.pi1_as<aod::TracksWMc>(),       // pi <- Xic
											 candidate.bachelor_as<aod::TracksWMc>(),  // pi <- cascade
											 candidate.posTrack_as<aod::TracksWMc>(),  // p <- lambda
											 candidate.negTrack_as<aod::TracksWMc>()}; // pi <- lambda

			// get XicPlus as MC particle
			if (configs.matchDecayedPions && configs.matchInteractionsWithMaterials) {
				indexRecXicPlus = RecoDecay::getMatchedMCRec<false, true, false, true, true>(particles, arrayDaughters, Pdg::kXiCPlus, std::array{+kPiPlus, +kPiPlus, +kPiMinus, +kProton, +kPiMinus}, true, nullptr, 4);
			} else if (configs.matchDecayedPions && !configs.matchInteractionsWithMaterials) {
				indexRecXicPlus = RecoDecay::getMatchedMCRec<false, true, false, true, false>(particles, arrayDaughters, Pdg::kXiCPlus, std::array{+kPiPlus, +kPiPlus, +kPiMinus, +kProton, +kPiMinus}, true, nullptr, 4);
			} else if (!configs.matchDecayedPions && configs.matchInteractionsWithMaterials) {
				indexRecXicPlus = RecoDecay::getMatchedMCRec<false, true, false, false, true>(particles, arrayDaughters, Pdg::kXiCPlus, std::array{+kPiPlus, +kPiPlus, +kPiMinus, +kProton, +kPiMinus}, true, nullptr, 4);
			} else {
				indexRecXicPlus = RecoDecay::getMatchedMCRec<false, true, false, false, false>(particles, arrayDaughters, Pdg::kXiCPlus, std::array{+kPiPlus, +kPiPlus, +kPiMinus, +kProton, +kPiMinus}, true, nullptr, 4);
			}
			if (indexRecXicPlus == -1) {
				continue;
			}
			auto particleXicPlusGen = particles.rawIteratorAt(indexRecXicPlus);
			origin = RecoDecay::getCharmHadronOrigin(particles, particleXicPlusGen, false);

			// get MC collision
			auto mcCollision = particleXicPlusGen.mcCollision_as<aod::McCollisions>();

			// get XicPlus daughters as MC particle
			RecoDecay::getDaughters(particleXicPlusGen, &arrDaughIndex, std::array{+kXiMinus, +kPiPlus, +kPiPlus}, 2);
			auto daugh0XicPlus = particles.rawIteratorAt(arrDaughIndex[0]);

			// calculate residuals and pulls
			float pResidual = candidate.p() - particleXicPlusGen.p();
			float ptResidual = candidate.pt() - particleXicPlusGen.pt();
			pvResiduals[0] = candidate.posX() - mcCollision.posX();
			pvResiduals[1] = candidate.posY() - mcCollision.posY();
			pvResiduals[2] = candidate.posZ() - mcCollision.posZ();
			svResiduals[0] = candidate.xSecondaryVertex() - daugh0XicPlus.vx();
			svResiduals[1] = candidate.ySecondaryVertex() - daugh0XicPlus.vy();
			svResiduals[2] = candidate.zSecondaryVertex() - daugh0XicPlus.vz();
			try {
				pvPulls[0] = pvResiduals[0] / candidate.xPvErr();
				pvPulls[1] = pvResiduals[1] / candidate.yPvErr();
				pvPulls[2] = pvResiduals[2] / candidate.zPvErr();
				svPulls[0] = svResiduals[0] / candidate.xSvErr();
				svPulls[1] = svResiduals[1] / candidate.ySvErr();
				svPulls[2] = svResiduals[2] / candidate.zSvErr();
			} catch (const std::runtime_error& error) {
				LOG(info) << "Run time error found: " << error.what() << ". Set values of vertex pulls to -9999.9.";
			}

			// fill table
			cursors.rowCandXicpResiduals( origin,
										  pResidual,
										  ptResidual,
										  pvResiduals[0],
										  pvResiduals[1],
										  pvResiduals[2],
										  pvPulls[0],
										  pvPulls[1],
										  pvPulls[2],
										  svResiduals[0],
										  svResiduals[1],
										  svResiduals[2],
										  svPulls[0],
										  svPulls[1],
										  svPulls[2]);
		} // loop over reconstructed signal
	}
	PROCESS_SWITCH(HfTreeCreatorXic0XicpToHadronic, processResidualsXicp, "Process Xicp Residuals and pulls for both DCAFitter and KFParticle reconstruction", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
	return WorkflowSpec{adaptAnalysisTask<HfTreeCreatorXic0XicpToHadronic>(cfgc)};
}
