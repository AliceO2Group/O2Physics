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

/// \file treeCreatorXicToXiPiPi.cxx
/// \brief Writer of Ξc± → Ξ∓ π± π± candidates in the form of flat tables to be stored in TTrees.
///
/// \author Phil Lennart Stahlhut <phil.lennart.stahlhut@cern.ch>, Heidelberg University
/// \author Carolina Reetz <c.reetz@cern.ch>, Heidelberg University

#include <vector>

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

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
DECLARE_SOA_COLUMN(CandidateSelFlag, candidateSelFlag, int); //! Selection flag of candidate (output of candidateSelector)
// vertices
DECLARE_SOA_COLUMN(XPv, xPv, float);
DECLARE_SOA_COLUMN(YPv, yPv, float);
DECLARE_SOA_COLUMN(ZPv, zPv, float);
DECLARE_SOA_COLUMN(XPvErr, xPvErr, float);
DECLARE_SOA_COLUMN(YPvErr, yPvErr, float);
DECLARE_SOA_COLUMN(ZPvErr, zPvErr, float);
DECLARE_SOA_COLUMN(XPvGen, xPvGen, float);
DECLARE_SOA_COLUMN(YPvGen, yPvGen, float);
DECLARE_SOA_COLUMN(ZPvGen, zPvGen, float);
DECLARE_SOA_COLUMN(XSv, xSv, float);
DECLARE_SOA_COLUMN(YSv, ySv, float);
DECLARE_SOA_COLUMN(ZSv, zSv, float);
DECLARE_SOA_COLUMN(XSvErr, xSvErr, float);
DECLARE_SOA_COLUMN(YSvErr, ySvErr, float);
DECLARE_SOA_COLUMN(ZSvErr, zSvErr, float);
DECLARE_SOA_COLUMN(Chi2Sv, chi2Sv, float);
DECLARE_SOA_COLUMN(XSvGen, xSvGen, float);
DECLARE_SOA_COLUMN(YSvGen, ySvGen, float);
DECLARE_SOA_COLUMN(ZSvGen, zSvGen, float);
DECLARE_SOA_COLUMN(XDecVtxXi, xDecVtxXi, float);
DECLARE_SOA_COLUMN(YDecVtxXi, yDecVtxXi, float);
DECLARE_SOA_COLUMN(ZDecVtxXi, zDecVtxXi, float);
DECLARE_SOA_COLUMN(Chi2XiVtx, chi2XiVtx, float);
DECLARE_SOA_COLUMN(XDecVtxLam, xDecVtxLam, float);
DECLARE_SOA_COLUMN(YDecVtxLam, yDecVtxLam, float);
DECLARE_SOA_COLUMN(ZDecVtxLam, zDecVtxLam, float);
DECLARE_SOA_COLUMN(Chi2LamVtx, chi2LamVtx, float);
// properties of XicPlus
DECLARE_SOA_COLUMN(Sign, sign, float);
DECLARE_SOA_COLUMN(E, e, float);                                             //! Energy of candidate (GeV)
DECLARE_SOA_COLUMN(M, m, float);                                             //! Invariant mass of candidate (GeV/c2)
DECLARE_SOA_COLUMN(P, p, float);                                             //! Momentum of candidate (GeV/c)
DECLARE_SOA_COLUMN(Pt, pt, float);                                           //! Transverse momentum of candidate (GeV/c)
DECLARE_SOA_COLUMN(Y, y, float);                                             //! Rapidity of candidate
DECLARE_SOA_COLUMN(Eta, eta, float);                                         //! Pseudorapidity of candidate
DECLARE_SOA_COLUMN(Phi, phi, float);                                         //! Azimuth angle of candidate
DECLARE_SOA_COLUMN(Ct, ct, float);                                           //! Proper lifetime time ctau of candidate (cm)
DECLARE_SOA_COLUMN(DecayLength, decayLength, float);                         //! Decay length of candidate (cm)
DECLARE_SOA_COLUMN(DecayLengthXY, decayLengthXY, float);                     //! Transverse decay length of candidate (cm)
DECLARE_SOA_COLUMN(DecayLengthNormalised, decayLengthNormalised, float);     //! Normalised decay length of candidate
DECLARE_SOA_COLUMN(DecayLengthXYNormalised, decayLengthXYNormalised, float); //! Normalised transverse decay length of candidate
DECLARE_SOA_COLUMN(Cpa, cpa, float);                                         //! Cosine pointing angle of candidate
DECLARE_SOA_COLUMN(CpaXY, cpaXY, float);                                     //! Cosine pointing angle of candidate in transverse plane
DECLARE_SOA_COLUMN(Chi2TopoXicPlusToPVBeforeConstraint, chi2TopoXicPlusToPVBeforeConstraint, float);
DECLARE_SOA_COLUMN(Chi2TopoXicPlusToPV, chi2TopoXicPlusToPV, float);
DECLARE_SOA_COLUMN(Chi2TopoXiToXicPlusBeforeConstraint, chi2TopoXiToXicPlusBeforeConstraint, float);
DECLARE_SOA_COLUMN(Chi2TopoXiToXicPlus, chi2TopoXiToXicPlus, float);
// properties of daughter tracks
DECLARE_SOA_COLUMN(PtXi, ptXi, float);                                                 //! Transverse momentum of Xi (prong0) (GeV/c)
DECLARE_SOA_COLUMN(ImpactParameterXi, impactParameterXi, float);                       //! Impact parameter of Xi (prong0)
DECLARE_SOA_COLUMN(ImpactParameterNormalisedXi, impactParameterNormalisedXi, float);   //! Normalised impact parameter of Xi (prong0)
DECLARE_SOA_COLUMN(PtPi0, ptPi0, float);                                               //! Transverse momentum of Pi0 (prong1) (GeV/c)
DECLARE_SOA_COLUMN(ImpactParameterPi0, impactParameterPi0, float);                     //! Impact parameter of Pi0 (prong1)
DECLARE_SOA_COLUMN(ImpactParameterNormalisedPi0, impactParameterNormalisedPi0, float); //! Normalised impact parameter of Pi0 (prong1)
DECLARE_SOA_COLUMN(PtPi1, ptPi1, float);                                               //! Transverse momentum of Pi1 (prong2) (GeV/c)
DECLARE_SOA_COLUMN(ImpactParameterPi1, impactParameterPi1, float);                     //! Normalised impact parameter of Pi1 (prong2)
DECLARE_SOA_COLUMN(ImpactParameterNormalisedPi1, impactParameterNormalisedPi1, float); //! Normalised impact parameter of Pi1 (prong2)
DECLARE_SOA_COLUMN(MaxNormalisedDeltaIP, maxNormalisedDeltaIP, float);                 //! Maximum normalized difference between measured and expected impact parameter of candidate prongs
DECLARE_SOA_COLUMN(CpaXi, cpaXi, float);
DECLARE_SOA_COLUMN(CpaXYXi, cpaXYXi, float);
DECLARE_SOA_COLUMN(CpaLam, cpaLam, float);
DECLARE_SOA_COLUMN(CpaXYLam, cpaXYLam, float);
DECLARE_SOA_COLUMN(DcaXYPi0Pi1, dcaXYPi0Pi1, float);
DECLARE_SOA_COLUMN(DcaXYPi0Xi, dcaXYPi0Xi, float);
DECLARE_SOA_COLUMN(DcaXYPi1Xi, dcaXYPi1Xi, float);
DECLARE_SOA_COLUMN(DcaPi0Pi1, dcaPi0Pi1, float);
DECLARE_SOA_COLUMN(DcaPi0Xi, dcaPi0Xi, float);
DECLARE_SOA_COLUMN(DcaPi1Xi, dcaPi1Xi, float);
DECLARE_SOA_COLUMN(DcaXiDaughters, dcaXiDaughters, float);
DECLARE_SOA_COLUMN(InvMassXiPi0, invMassXiPi0, float);
DECLARE_SOA_COLUMN(InvMassXiPi1, invMassXiPi1, float);
// residuals and pulls
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
} // namespace full

DECLARE_SOA_TABLE(HfCandXicToXiPiPiLites, "AOD", "HFXICXI2PILITE",
                  hf_cand_xic_to_xi_pi_pi::FlagMcMatchRec,
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
                  hf_cand_xic_to_xi_pi_pi::FlagMcMatchRec,
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
                  full::Chi2TopoXiToXicPlusBeforeConstraint,
                  full::Chi2TopoXiToXicPlus,
                  full::DcaXYPi0Pi1,
                  full::DcaXYPi0Xi,
                  full::DcaXYPi1Xi,
                  full::DcaPi0Pi1,
                  full::DcaPi0Xi,
                  full::DcaPi1Xi,
                  full::DcaXiDaughters);

DECLARE_SOA_TABLE(HfCandXicToXiPiPiFulls, "AOD", "HFXICXI2PIFULL",
                  hf_cand_xic_to_xi_pi_pi::FlagMcMatchRec,
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
                  full::E,
                  full::XPv,
                  full::YPv,
                  full::ZPv,
                  full::XPvErr,
                  full::YPvErr,
                  full::ZPvErr,
                  full::XSv,
                  full::YSv,
                  full::ZSv,
                  full::XSvErr,
                  full::YSvErr,
                  full::ZSvErr,
                  full::XDecVtxXi,
                  full::YDecVtxXi,
                  full::ZDecVtxXi,
                  full::XDecVtxLam,
                  full::YDecVtxLam,
                  full::ZDecVtxLam);

DECLARE_SOA_TABLE(HfCandXicToXiPiPiFullKfs, "AOD", "HFXICXI2PIFULKF",
                  hf_cand_xic_to_xi_pi_pi::FlagMcMatchRec,
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
                  full::E,
                  full::XPv,
                  full::YPv,
                  full::ZPv,
                  full::XPvErr,
                  full::YPvErr,
                  full::ZPvErr,
                  full::XSv,
                  full::YSv,
                  full::ZSv,
                  full::XSvErr,
                  full::YSvErr,
                  full::ZSvErr,
                  full::XDecVtxXi,
                  full::YDecVtxXi,
                  full::ZDecVtxXi,
                  full::XDecVtxLam,
                  full::YDecVtxLam,
                  full::ZDecVtxLam,
                  // KF-specific columns
                  full::Chi2XiVtx,
                  full::Chi2LamVtx,
                  full::Chi2TopoXicPlusToPVBeforeConstraint,
                  full::Chi2TopoXicPlusToPV,
                  full::Chi2TopoXiToXicPlusBeforeConstraint,
                  full::Chi2TopoXiToXicPlus,
                  full::DcaXYPi0Pi1,
                  full::DcaXYPi0Xi,
                  full::DcaXYPi1Xi,
                  full::DcaPi0Pi1,
                  full::DcaPi0Xi,
                  full::DcaPi1Xi,
                  full::DcaXiDaughters);

DECLARE_SOA_TABLE(HfCandXicToXiPiPiFullPs, "AOD", "HFXICXI2PIFULLP",
                  full::Pt,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  full::XPvGen,
                  full::YPvGen,
                  full::ZPvGen,
                  full::XSvGen,
                  full::YSvGen,
                  full::ZSvGen,
                  hf_cand_xic_to_xi_pi_pi::FlagMcMatchGen);

DECLARE_SOA_TABLE(HfCandXicToXiPiPiResiduals, "AOD", "HFXICXI2PIRESID",
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
} // namespace o2::aod

/// Writes the full information in an output TTree
struct HfTreeCreatorXicToXiPiPi {
  Produces<o2::aod::HfCandXicToXiPiPiLites> rowCandidateLite;
  Produces<o2::aod::HfCandXicToXiPiPiLiteKfs> rowCandidateLiteKf;
  Produces<o2::aod::HfCandXicToXiPiPiFulls> rowCandidateFull;
  Produces<o2::aod::HfCandXicToXiPiPiFullKfs> rowCandidateFullKf;
  Produces<o2::aod::HfCandXicToXiPiPiFullPs> rowCandidateFullParticles;
  Produces<o2::aod::HfCandXicToXiPiPiResiduals> rowCandidateResiduals;

  Configurable<int> selectionFlagXic{"selectionFlagXic", 1, "Selection Flag for Xic"};
  Configurable<bool> fillCandidateLiteTable{"fillCandidateLiteTable", false, "Switch to fill lite table with candidate properties"};
  Configurable<bool> fillGenParticleTable{"fillGenParticleTable", false, "Switch to fill table with MC truth for generated particles"};
  // parameters for production of training samples
  Configurable<bool> fillOnlySignal{"fillOnlySignal", false, "Flag to fill derived tables with signal for ML trainings"};
  Configurable<bool> fillOnlyBackground{"fillOnlyBackground", false, "Flag to fill derived tables with background for ML trainings"};
  Configurable<float> downSampleBkgFactor{"downSampleBkgFactor", 1., "Fraction of background candidates to keep for ML trainings"};
  Configurable<float> ptMaxForDownSample{"ptMaxForDownSample", 10., "Maximum pt for the application of the downsampling factor"};

  using SelectedCandidates = soa::Filtered<soa::Join<aod::HfCandXic, aod::HfSelXicToXiPiPi>>;
  using SelectedCandidatesKf = soa::Filtered<soa::Join<aod::HfCandXic, aod::HfCandXicKF, aod::HfSelXicToXiPiPi>>;
  using SelectedCandidatesMc = soa::Filtered<soa::Join<aod::HfCandXic, aod::HfCandXicMcRec, aod::HfSelXicToXiPiPi>>;
  using SelectedCandidatesKfMc = soa::Filtered<soa::Join<aod::HfCandXic, aod::HfCandXicKF, aod::HfCandXicMcRec, aod::HfSelXicToXiPiPi>>;

  Filter filterSelectCandidates = aod::hf_sel_candidate_xic::isSelXicToXiPiPi >= selectionFlagXic;

  Partition<SelectedCandidatesMc> recSig = nabs(aod::hf_cand_xic_to_xi_pi_pi::flagMcMatchRec) != int8_t(0);
  Partition<SelectedCandidatesMc> recBg = nabs(aod::hf_cand_xic_to_xi_pi_pi::flagMcMatchRec) == int8_t(0);
  Partition<SelectedCandidatesKfMc> recSigKf = nabs(aod::hf_cand_xic_to_xi_pi_pi::flagMcMatchRec) != int8_t(0);
  Partition<SelectedCandidatesKfMc> recBgKf = nabs(aod::hf_cand_xic_to_xi_pi_pi::flagMcMatchRec) == int8_t(0);

  void init(InitContext const&)
  {
  }

  template <bool doMc, bool doKf, typename T>
  void fillCandidateTable(const T& candidate)
  {
    int8_t flagMc = 0;
    if constexpr (doMc) {
      flagMc = candidate.flagMcMatchRec();
    }
    if constexpr (!doKf) {
      if (fillCandidateLiteTable) {
        rowCandidateLite(
          flagMc,
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
          candidate.invMassXic(),
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
          candidate.cosPaXi(),
          candidate.cosPaXYXi(),
          candidate.cosPaLambda(),
          candidate.cosPaXYLambda(),
          candidate.impactParameter0(),
          candidate.impactParameterNormalised0(),
          candidate.impactParameter1(),
          candidate.impactParameterNormalised1(),
          candidate.impactParameter2(),
          candidate.impactParameterNormalised2(),
          candidate.maxNormalisedDeltaIP());
      } else {
        rowCandidateFull(
          flagMc,
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
          candidate.invMassXic(),
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
          candidate.cosPaXi(),
          candidate.cosPaXYXi(),
          candidate.cosPaLambda(),
          candidate.cosPaXYLambda(),
          candidate.impactParameter0(),
          candidate.impactParameterNormalised0(),
          candidate.impactParameter1(),
          candidate.impactParameterNormalised1(),
          candidate.impactParameter2(),
          candidate.impactParameterNormalised2(),
          candidate.maxNormalisedDeltaIP(),
          // additional columns only stored in the full candidate table
          candidate.e(o2::constants::physics::MassXiCPlus),
          candidate.posX(),
          candidate.posY(),
          candidate.posZ(),
          candidate.xPvErr(),
          candidate.yPvErr(),
          candidate.zPvErr(),
          candidate.xSecondaryVertex(),
          candidate.ySecondaryVertex(),
          candidate.zSecondaryVertex(),
          candidate.xSvErr(),
          candidate.ySvErr(),
          candidate.zSvErr(),
          candidate.xDecayVtxXi(),
          candidate.yDecayVtxXi(),
          candidate.zDecayVtxXi(),
          candidate.xDecayVtxLambda(),
          candidate.yDecayVtxLambda(),
          candidate.zDecayVtxLambda());
      }
    } else {
      if (fillCandidateLiteTable) {
        rowCandidateLiteKf(
          flagMc,
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
          candidate.invMassXic(),
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
          candidate.cosPaXi(),
          candidate.cosPaXYXi(),
          candidate.cosPaLambda(),
          candidate.cosPaXYLambda(),
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
          candidate.chi2TopoXicPlusToPVBeforeConstraint(),
          candidate.chi2TopoXicPlusToPV(),
          candidate.chi2TopoXiToXicPlusBeforeConstraint(),
          candidate.chi2TopoXiToXicPlus(),
          candidate.dcaXYPi0Pi1(),
          candidate.dcaXYPi0Xi(),
          candidate.dcaXYPi1Xi(),
          candidate.dcaPi0Pi1(),
          candidate.dcaPi0Xi(),
          candidate.dcaPi1Xi(),
          candidate.dcacascdaughters());
      } else {
        rowCandidateFullKf(
          flagMc,
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
          candidate.invMassXic(),
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
          candidate.cosPaXi(),
          candidate.cosPaXYXi(),
          candidate.cosPaLambda(),
          candidate.cosPaXYLambda(),
          candidate.impactParameter0(),
          candidate.impactParameterNormalised0(),
          candidate.impactParameter1(),
          candidate.impactParameterNormalised1(),
          candidate.impactParameter2(),
          candidate.impactParameterNormalised2(),
          candidate.maxNormalisedDeltaIP(),
          // additional columns only stored in the full candidate table
          candidate.e(o2::constants::physics::MassXiCPlus),
          candidate.posX(),
          candidate.posY(),
          candidate.posZ(),
          candidate.xPvErr(),
          candidate.yPvErr(),
          candidate.zPvErr(),
          candidate.xSecondaryVertex(),
          candidate.ySecondaryVertex(),
          candidate.zSecondaryVertex(),
          candidate.xSvErr(),
          candidate.ySvErr(),
          candidate.zSvErr(),
          candidate.xDecayVtxXi(),
          candidate.yDecayVtxXi(),
          candidate.zDecayVtxXi(),
          candidate.xDecayVtxLambda(),
          candidate.yDecayVtxLambda(),
          candidate.zDecayVtxLambda(),
          // KF-specific columns
          candidate.kfCascadeChi2(),
          candidate.kfV0Chi2(),
          candidate.chi2TopoXicPlusToPVBeforeConstraint(),
          candidate.chi2TopoXicPlusToPV(),
          candidate.chi2TopoXiToXicPlusBeforeConstraint(),
          candidate.chi2TopoXiToXicPlus(),
          candidate.dcaXYPi0Pi1(),
          candidate.dcaXYPi0Xi(),
          candidate.dcaXYPi1Xi(),
          candidate.dcaPi0Pi1(),
          candidate.dcaPi0Xi(),
          candidate.dcaPi1Xi(),
          candidate.dcacascdaughters());
      }
    }
  }

  void processData(SelectedCandidates const& candidates)
  {
    // Filling candidate properties
    if (fillCandidateLiteTable) {
      rowCandidateLite.reserve(candidates.size());
    } else {
      rowCandidateFull.reserve(candidates.size());
    }
    for (const auto& candidate : candidates) {
      if (fillOnlyBackground && downSampleBkgFactor < 1.) {
        float pseudoRndm = candidate.ptProng1() * 1000. - static_cast<int64_t>(candidate.ptProng1() * 1000);
        if (pseudoRndm >= downSampleBkgFactor && candidate.pt() < ptMaxForDownSample) {
          continue;
        }
      }
      fillCandidateTable<false, false>(candidate);
    }
  }
  PROCESS_SWITCH(HfTreeCreatorXicToXiPiPi, processData, "Process data with DCAFitter reconstruction", true);

  void processDataKf(SelectedCandidatesKf const& candidates)
  {
    // Filling candidate properties
    if (fillCandidateLiteTable) {
      rowCandidateLite.reserve(candidates.size());
    } else {
      rowCandidateFull.reserve(candidates.size());
    }
    for (const auto& candidate : candidates) {
      if (fillOnlyBackground && downSampleBkgFactor < 1.) {
        float pseudoRndm = candidate.ptProng1() * 1000. - static_cast<int64_t>(candidate.ptProng1() * 1000);
        if (pseudoRndm >= downSampleBkgFactor && candidate.pt() < ptMaxForDownSample) {
          continue;
        }
      }
      fillCandidateTable<false, true>(candidate);
    }
  }
  PROCESS_SWITCH(HfTreeCreatorXicToXiPiPi, processDataKf, "Process data with KFParticle reconstruction", false);

  void processMc(SelectedCandidatesMc const& candidates,
                 soa::Join<aod::McParticles, aod::HfCandXicMcGen> const& particles)
  {
    std::vector<int> arrDaughIndex;

    // Filling candidate properties
    if (fillOnlySignal) {
      if (fillCandidateLiteTable) {
        rowCandidateLite.reserve(recSig.size());
      } else {
        rowCandidateFull.reserve(recSig.size());
      }
      for (const auto& candidate : recSig) {
        fillCandidateTable<true, false>(candidate);
      }
    } else if (fillOnlyBackground) {
      if (fillCandidateLiteTable) {
        rowCandidateLite.reserve(recBg.size());
      } else {
        rowCandidateFull.reserve(recBg.size());
      }
      for (const auto& candidate : recBg) {
        float pseudoRndm = candidate.ptProng1() * 1000. - static_cast<int64_t>(candidate.ptProng1() * 1000);
        if (candidate.pt() < ptMaxForDownSample && pseudoRndm >= downSampleBkgFactor) {
          continue;
        }
        fillCandidateTable<true, false>(candidate);
      }
    } else {
      if (fillCandidateLiteTable) {
        rowCandidateLite.reserve(candidates.size());
      } else {
        rowCandidateFull.reserve(candidates.size());
      }
      for (const auto& candidate : candidates) {
        fillCandidateTable<true, false>(candidate);
      }
    }

    if (fillGenParticleTable) {
      rowCandidateFullParticles.reserve(particles.size());

      for (const auto& particle : particles) {
        if (TESTBIT(std::abs(particle.flagMcMatchGen()), aod::hf_cand_xic_to_xi_pi_pi::DecayType::XicToXiPiPi) || TESTBIT(std::abs(particle.flagMcMatchGen()), aod::hf_cand_xic_to_xi_pi_pi::DecayType::XicToXiResPiToXiPiPi)) {
          arrDaughIndex.clear();
          RecoDecay::getDaughters(particle, &arrDaughIndex, std::array{+kXiMinus, +kPiPlus, +kPiPlus}, 2);
          auto XicDaugh0 = particles.rawIteratorAt(arrDaughIndex[0]);

          rowCandidateFullParticles(
            particle.pt(),
            particle.eta(),
            particle.phi(),
            RecoDecay::y(particle.pVector(), o2::constants::physics::MassXiCPlus),
            particle.vx(),
            particle.vy(),
            particle.vz(),
            XicDaugh0.vx(),
            XicDaugh0.vx(),
            XicDaugh0.vz(),
            particle.flagMcMatchGen());
        }
      } // loop over generated particles
    }
  }
  PROCESS_SWITCH(HfTreeCreatorXicToXiPiPi, processMc, "Process MC with DCAFitter reconstruction", false);

  void processMcKf(SelectedCandidatesKfMc const& candidates,
                   soa::Join<aod::McParticles, aod::HfCandXicMcGen> const& particles)
  {
    std::vector<int> arrDaughIndex;

    // Filling candidate properties
    if (fillOnlySignal) {
      if (fillCandidateLiteTable) {
        rowCandidateLite.reserve(recSigKf.size());
      } else {
        rowCandidateFull.reserve(recSigKf.size());
      }
      for (const auto& candidate : recSigKf) {
        fillCandidateTable<true, true>(candidate);
      }
    } else if (fillOnlyBackground) {
      if (fillCandidateLiteTable) {
        rowCandidateLite.reserve(recBgKf.size());
      } else {
        rowCandidateFull.reserve(recBgKf.size());
      }
      for (const auto& candidate : recBgKf) {
        float pseudoRndm = candidate.ptProng1() * 1000. - static_cast<int64_t>(candidate.ptProng1() * 1000);
        if (candidate.pt() < ptMaxForDownSample && pseudoRndm >= downSampleBkgFactor) {
          continue;
        }
        fillCandidateTable<true, true>(candidate);
      }
    } else {
      if (fillCandidateLiteTable) {
        rowCandidateLite.reserve(candidates.size());
      } else {
        rowCandidateFull.reserve(candidates.size());
      }
      for (const auto& candidate : candidates) {
        fillCandidateTable<true, true>(candidate);
      }
    }

    if (fillGenParticleTable) {
      rowCandidateFullParticles.reserve(particles.size());
      for (const auto& particle : particles) {
        if (TESTBIT(std::abs(particle.flagMcMatchGen()), aod::hf_cand_xic_to_xi_pi_pi::DecayType::XicToXiPiPi) || TESTBIT(std::abs(particle.flagMcMatchGen()), aod::hf_cand_xic_to_xi_pi_pi::DecayType::XicToXiResPiToXiPiPi)) {
          arrDaughIndex.clear();
          RecoDecay::getDaughters(particle, &arrDaughIndex, std::array{+kXiMinus, +kPiPlus, +kPiPlus}, 2);
          auto XicDaugh0 = particles.rawIteratorAt(arrDaughIndex[0]);

          rowCandidateFullParticles(
            particle.pt(),
            particle.eta(),
            particle.phi(),
            RecoDecay::y(particle.pVector(), o2::constants::physics::MassXiCPlus),
            particle.vx(),
            particle.vy(),
            particle.vz(),
            XicDaugh0.vx(),
            XicDaugh0.vx(),
            XicDaugh0.vz(),
            particle.flagMcMatchGen());
        }
      } // loop over generated particles
    }
  }
  PROCESS_SWITCH(HfTreeCreatorXicToXiPiPi, processMcKf, "Process MC with KF Particle reconstruction", false);

  void processResiduals(SelectedCandidatesMc const&,
                        aod::TracksWMc const& tracks,
                        aod::McParticles const& particles)
  {
    rowCandidateResiduals.reserve(recSig.size());

    recSig->bindExternalIndices(&tracks);

    std::vector<int> arrDaughIndex;
    int indexRecXic;
    int8_t sign;
    std::array<float, 3> pvResiduals;
    std::array<float, 3> svResiduals;
    std::array<float, 3> pvPulls;
    std::array<float, 3> svPulls;

    for (const auto& candidate : recSig) {
      arrDaughIndex.clear();
      indexRecXic = -1;
      sign = 0;
      pvResiduals = {-9999.9};
      svResiduals = {-9999.9};
      pvPulls = {-9999.9};
      svPulls = {-9999.9};

      auto arrayDaughters = std::array{candidate.pi0_as<aod::TracksWMc>(),       // pi <- Xic
                                       candidate.pi1_as<aod::TracksWMc>(),       // pi <- Xic
                                       candidate.bachelor_as<aod::TracksWMc>(),  // pi <- cascade
                                       candidate.posTrack_as<aod::TracksWMc>(),  // p <- lambda
                                       candidate.negTrack_as<aod::TracksWMc>()}; // pi <- lambda

      // get Xic and daughters as MC particle
      indexRecXic = RecoDecay::getMatchedMCRec(particles, arrayDaughters, Pdg::kXiCPlus, std::array{+kPiPlus, +kPiPlus, +kPiMinus, +kProton, +kPiMinus}, true, &sign, 4);
      if (indexRecXic == -1) {
        continue;
      }
      auto XicGen = particles.rawIteratorAt(indexRecXic);
      RecoDecay::getDaughters(XicGen, &arrDaughIndex, std::array{+kXiMinus, +kPiPlus, +kPiPlus}, 2);
      auto XicDaugh0 = particles.rawIteratorAt(arrDaughIndex[0]);

      // calculate residuals and pulls
      float pResidual = candidate.p() - XicGen.p();
      float ptResidual = candidate.pt() - XicGen.pt();
      pvResiduals[0] = candidate.posX() - XicGen.vx();
      pvResiduals[1] = candidate.posY() - XicGen.vy();
      pvResiduals[2] = candidate.posZ() - XicGen.vz();
      svResiduals[0] = candidate.xSecondaryVertex() - XicDaugh0.vx();
      svResiduals[1] = candidate.ySecondaryVertex() - XicDaugh0.vy();
      svResiduals[2] = candidate.zSecondaryVertex() - XicDaugh0.vz();
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
      rowCandidateResiduals(
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
  PROCESS_SWITCH(HfTreeCreatorXicToXiPiPi, processResiduals, "Process Residuals and pulls for both DCAFitter and KFParticle reconstruction", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTreeCreatorXicToXiPiPi>(cfgc)};
}
