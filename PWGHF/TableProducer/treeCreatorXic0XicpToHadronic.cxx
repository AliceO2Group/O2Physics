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
// residuals and pulls -> What does this mean?
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
}

DECLARE_SOA_TABLE(HfCandXic0Lites, "AOD", "HFCANDXIC0LITE",
				  hf_cand_xic0_xicp_to_hadronic::FlagMcMatchRec,
				  hf_cand_xic0_xicp_to_hadronic::DebugMcRec,
				  hf_cand_xic0_xicp_to_hadronic::OriginRec,
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
				  hf_cand_xic0_xicp_to_hadronic::FlagMcMatchRec,
				  hf_cand_xic0_xicp_to_hadronic::DebugMcRec,
				  hf_cand_xic0_xicp_to_hadronic::OriginRec,
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
				  hf_cand_xic0_xicp_to_hadronic::FlagMcMatchRec,
				  hf_cand_xic0_xicp_to_hadronic::DebugMcRec,
				  hf_cand_xic0_xicp_to_hadronic::OriginRec,
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
				  hf_cand_xic0_xicp_to_hadronic::DcaXiDaughters,
				  hf_cand_xic0_xicp_to_hadronic::DcaV0Daughters,
				  hf_cand_xic0_xicp_to_hadronic::DcaPosToPV,
				  hf_cand_xic0_xicp_to_hadronic::DcaNegToPV,
				  hf_cand_xic0_xicp_to_hadronic::DcaBachelorToPV,
				  hf_cand_xic0_xicp_to_hadronic::DcaXYCascToPV,
				  hf_cand_xic0_xicp_to_hadronic::DcaZCascToPV,
				  hf_cand_xic0_xicp_to_hadronic::NSigTpcPiFromXic0,
				  hf_cand_xic0_xicp_to_hadronic::NSigTpcBachelorPi,
				  hf_cand_xic0_xicp_to_hadronic::NSigTpcPiFromLambda,
				  hf_cand_xic0_xicp_to_hadronic::NSigTpcPrFromLambda,
				  hf_cand_xic0_xicp_to_hadronic::NSigTofPiFromXic0,
				  hf_cand_xic0_xicp_to_hadronic::NSigTofBachelorPi,
				  hf_cand_xic0_xicp_to_hadronic::NSigTofPiFromLambda,
				  hf_cand_xic0_xicp_to_hadronic::NSigTofPrFromLambda);

DECLARE_SOA_TABLE(HfCandXic0FullKfs, "AOD", "HFCANDXIC0FULKF",
				  hf_cand_xic0_xicp_to_hadronic::FlagMcMatchRec,
				  hf_cand_xic0_xicp_to_hadronic::DebugMcRec,
				  hf_cand_xic0_xicp_to_hadronic::OriginRec,
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
				  hf_cand_xic0_xicp_to_hadronic::DcaXiDaughters,
				  hf_cand_xic0_xicp_to_hadronic::DcaV0Daughters,
				  hf_cand_xic0_xicp_to_hadronic::DcaPosToPV,
				  hf_cand_xic0_xicp_to_hadronic::DcaNegToPV,
				  hf_cand_xic0_xicp_to_hadronic::DcaBachelorToPV,
				  hf_cand_xic0_xicp_to_hadronic::DcaXYCascToPV,
				  hf_cand_xic0_xicp_to_hadronic::DcaZCascToPV,
				  hf_cand_xic0_xicp_to_hadronic::NSigTpcPiFromXic0,
				  hf_cand_xic0_xicp_to_hadronic::NSigTpcBachelorPi,
				  hf_cand_xic0_xicp_to_hadronic::NSigTpcPiFromLambda,
				  hf_cand_xic0_xicp_to_hadronic::NSigTpcPrFromLambda,
				  hf_cand_xic0_xicp_to_hadronic::NSigTofPiFromXic0,
				  hf_cand_xic0_xicp_to_hadronic::NSigTofBachelorPi,
				  hf_cand_xic0_xicp_to_hadronic::NSigTofPiFromLambda,
				  hf_cand_xic0_xicp_to_hadronic::NSigTofPrFromLambda,
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
				  hf_cand_xic0_xicp_to_hadronic::FlagMcMatchGen,
				  hf_cand_xic0_xicp_to_hadronic::DebugMcGen,
				  hf_cand_xic0_xicp_to_hadronic::OriginGen,
				  full::Pt,
				  full::Eta,
				  full::Phi,
				  full::Y);

DECLARE_SOA_TABLE(HfCandXic0ToXiPiResiduals, "AOD", "HFXIC0TOXIPIRESID",
				  hf_cand_xic0_xicp_to_hadronic::OriginGen,
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
	} cursors;
	
	struct : ConfigurableGroup {
		Configurable<int> selectionFlagXic0{"selectionFlagXic0", 1, "Selection flag for Xic0"};
		Configurable<bool> fillCandidateLite{"fillCandidateLite", true, "Switch to fill lite table with candidate properties"};
		Configurable<bool> fillGenParticleTable{"fillGenParticleTable", false, "Switch to fill table with MC truth for generated particles"};
		Configurable<bool> fillOnlySignal{"fillOnlySignal", false, "Flag to fill derived tables with signal for ML trainings"};
		Configurable<bool> fillOnlyBackground{"fillOnlyBackground", false, "Flag to fill derived tables with background for ML trainings"};
		Configurable<float> downSampleBkgFactor{"downSampleBkgFactor", 1., "Fraction of background candidates to keep for ML trainings"};// ?
		Configurable<float> ptMaxForDownSample{"ptMaxForDownSample", 10., "Maximim pt for the application of the downsampling factor"}; // ?
	} configs;

	using SelectedCandidates = soa::Filtered<soa::Join<aod::HfCandXic0, aod::HfSelXic0ToXiPi>>;
	using SelectedCandidatesKf = soa::Filtered<soa::Join<aod::HfCandXic0, aod::HfCandXic0KF, aod::HfSelXic0ToXiPi>>;
	using SelectedCandidatesMc = soa::Filtered<soa::Join<aod::HfCandXic0, aod::HfCandXic0McRec, aod::HfSelXic0ToXiPi>>;
	using SelectedCandidatesKfMc = soa::Filtered<soa::Join<aod::HfCandXic0, aod::HfCandXic0KF, aod::HfCandXic0McRec, aod::HfSelXic0ToXiPi>>;

	Filter filterSelectCandidates = aod::hf_sel_xic0_xicp_to_hadronic::isSelXic0ToXiPi >= configs.selectionFlagXic0;

	Partition<SelectedCandidatesMc> recSig = nabs(aod::hf_cand_xic0_xicp_to_hadronic::flagMcMatchRec) != int8_t(0);
	Partition<SelectedCandidatesMc> recBkg= nabs(aod::hf_cand_xic0_xicp_to_hadronic::flagMcMatchRec) == int8_t(0);
	Partition<SelectedCandidatesKfMc> recSigKf = nabs(aod::hf_cand_xic0_xicp_to_hadronic::flagMcMatchRec) != int8_t(0);
	Partition<SelectedCandidatesKfMc> recBkgKf = nabs(aod::hf_cand_xic0_xicp_to_hadronic::flagMcMatchRec) == int8_t(0);

	void init(InitContext const&)
	{
	}

	template <bool doMc, bool doKf, typename T>
	void fillCandidateTable(const T& candidate)
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
	}


	////////////////////////////////////////////////////
	////											////
	////			Make TTree with Data			////
	////											////
	////////////////////////////////////////////////////


	void processData(SelectedCandidates const& candidates) // -> To be filled
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
			fillCandidateTable<false, false>(candidate);
		}
	}
	PROCESS_SWITCH(HfTreeCreatorXic0XicpToHadronic, processData, "Process Data with DCAFitter reconstruction", false);
	

	void processDataKf(SelectedCandidatesKf const& candidates)
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
			fillCandidateTable<false, true>(candidate);
		}
	}
	PROCESS_SWITCH(HfTreeCreatorXic0XicpToHadronic, processDataKf, "Process Data with KFParticle reconstruction", true);

	////////////////////////////////////////////////////
	////											////
	////			Make TTree with MC				////
	////											////
	////////////////////////////////////////////////////


	void processMc(SelectedCandidatesMc const& candidates,
				   soa::Join<aod::McParticles, aod::HfCandXic0McGen> const& particles) // -> To be filled
	{
		//! Fill only signal for ML Training
		//!! In the reference code, the cursors for DCAFitter was used instead....
		if (configs.fillOnlySignal) {

			if (configs.fillCandidateLite) {
				cursors.rowCandXic0Lites.reserve(recSig.size());
			} else {
				cursors.rowCandXic0Fulls.reserve(recSig.size());
			}

			for (const auto& candidate : recSig) {
				fillCandidateTable<true, false>(candidate);
			}

		//! Fill only Bkg for ML Training
		} else if (configs.fillOnlyBackground) {

			if (configs.fillCandidateLite) {
				cursors.rowCandXic0Lites.reserve(recBkg.size());
			} else {
				cursors.rowCandXic0Fulls.reserve(recBkg.size());
			}

			for (const auto& candidate : recBkg) {
				float pseudoRndm = candidate.ptProng1()*1000. - static_cast<int64_t>(candidate.ptProng1()*1000);
				if (candidate.pt() < configs.ptMaxForDownSample && pseudoRndm >= configs.downSampleBkgFactor) {
					continue;
				}
				fillCandidateTable<true, false>(candidate);
			}

		//! Fill whole reconstructed(?) candidates
		} else {

			if (configs.fillCandidateLite) {
				cursors.rowCandXic0Lites.reserve(candidates.size());
			} else {
				cursors.rowCandXic0Fulls.reserve(candidates.size());
			}

			for (auto const& candidate : candidates) {
				fillCandidateTable<true, false>(candidate);
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
	PROCESS_SWITCH(HfTreeCreatorXic0XicpToHadronic, processMc, "Process MC with DCAFitter reconstruction", false);


	void processMcKf(SelectedCandidatesKfMc const& candidates,
					 soa::Join<aod::McParticles, aod::HfCandXic0McGen> const& particles)
	{
		//! Fill only signal for ML Training
		//!! In the reference code, the cursors for DCAFitter was used instead....
		if (configs.fillOnlySignal) {

			if (configs.fillCandidateLite) {
				cursors.rowCandXic0LiteKfs.reserve(recSigKf.size());
			} else {
				cursors.rowCandXic0FullKfs.reserve(recSigKf.size());
			}

			for (const auto& candidate : recSigKf) {
				fillCandidateTable<true, true>(candidate);
			}

		//! Fill only Bkg for ML Training
		} else if (configs.fillOnlyBackground) {

			if (configs.fillCandidateLite) {
				cursors.rowCandXic0LiteKfs.reserve(recBkgKf.size());
			} else {
				cursors.rowCandXic0FullKfs.reserve(recBkgKf.size());
			}

			for (const auto& candidate : recBkgKf) {
				float pseudoRndm = candidate.ptProng1()*1000. - static_cast<int64_t>(candidate.ptProng1()*1000);
				if (candidate.pt() < configs.ptMaxForDownSample && pseudoRndm >= configs.downSampleBkgFactor) {
					continue;
				}
				fillCandidateTable<true, true>(candidate);
			}

		//! Fill whole reconstructed(?) candidates
		} else {

			if (configs.fillCandidateLite) {
				cursors.rowCandXic0LiteKfs.reserve(candidates.size());
			} else {
				cursors.rowCandXic0FullKfs.reserve(candidates.size());
			}

			for (auto const& candidate : candidates) {
				fillCandidateTable<true, true>(candidate);
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
	PROCESS_SWITCH(HfTreeCreatorXic0XicpToHadronic, processMcKf, "Process MC with DCAFitter reconstruction", false);

	void processResiduals(SelectedCandidatesMc const&,
						  aod::TracksWMc const&  tracks,
						  aod::McParticles const& particles,
						  aod::McCollisions const&)
	{
		cursors.rowCandResiduals.reserve(recSig.size());

		recSig->bindExternalIndices(&tracks); // -> What does this feature do exactly?
											  
		std::vector<int> arrDaughIdx;
		int idxRecXic0;
		int8_t sign;
		int8_t origin;
		std::array<float, 3> pvResiduals;
		std::array<float, 3> svResiduals;
		std::array<float, 3> pvPulls;
		std::array<float, 3> svPulls;
		
		for (auto const& candidate : recSig) {

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
	PROCESS_SWITCH(HfTreeCreatorXic0XicpToHadronic, processResiduals, "Process residuals and pulls for both DCAFitter and KFParticle reconstruction", false);

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
	return WorkflowSpec{adaptAnalysisTask<HfTreeCreatorXic0XicpToHadronic>(cfgc)};
}
