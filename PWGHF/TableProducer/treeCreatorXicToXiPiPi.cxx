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

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

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
DECLARE_SOA_COLUMN(Chi2XicPlusTopoToPV, chi2XicPlusTopoToPV, float);
DECLARE_SOA_COLUMN(Chi2XicPlusTopoXiToXicPlus, chi2XicPlusTopoXiToXicPlus, float);
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
                  full::CandidateSelFlag,
                  full::XPv,
                  full::YPv,
                  full::ZPv,
                  full::XSv,
                  full::YSv,
                  full::ZSv,
                  full::Chi2Sv,
                  full::Sign,
                  full::E,
                  full::M,
                  full::P,
                  full::Pt,
                  full::Y,
                  full::Eta,
                  full::Phi,
                  full::Ct,
                  full::DecayLength,
                  full::DecayLengthXY,
                  full::Cpa,
                  full::CpaXY,
                  full::PtXi,
                  full::PtPi0,
                  full::PtPi1,
                  full::ImpactParameterXi,
                  full::ImpactParameterPi0,
                  full::ImpactParameterPi1,
                  full::CpaXi,
                  full::CpaXYXi,
                  full::CpaLam,
                  full::CpaXYLam,
                  hf_cand_xic_to_xi_pi_pi::FlagMcMatchRec);

DECLARE_SOA_TABLE(HfCandXicToXiPiPiLiteKfs, "AOD", "HFXICXI2PILITKF",
                  full::CandidateSelFlag,
                  full::XPv,
                  full::YPv,
                  full::ZPv,
                  full::XSv,
                  full::YSv,
                  full::ZSv,
                  full::Chi2Sv,
                  full::Sign,
                  full::E,
                  full::M,
                  full::P,
                  full::Pt,
                  full::Y,
                  full::Eta,
                  full::Phi,
                  full::Ct,
                  full::DecayLength,
                  full::DecayLengthXY,
                  full::Cpa,
                  full::CpaXY,
                  full::PtXi,
                  full::PtPi0,
                  full::PtPi1,
                  full::ImpactParameterXi,
                  full::ImpactParameterPi0,
                  full::ImpactParameterPi1,
                  full::CpaXi,
                  full::CpaXYXi,
                  full::CpaLam,
                  full::CpaXYLam,
                  full::Chi2XiVtx,
                  full::Chi2LamVtx,
                  full::Chi2XicPlusTopoToPV,
                  full::Chi2XicPlusTopoXiToXicPlus,
                  full::DcaXYPi0Pi1,
                  full::DcaXYPi0Xi,
                  full::DcaXYPi1Xi,
                  full::DcaPi0Pi1,
                  full::DcaPi0Xi,
                  full::DcaPi1Xi,
                  full::DcaXiDaughters,
                  hf_cand_xic_to_xi_pi_pi::FlagMcMatchRec);

DECLARE_SOA_TABLE(HfCandXicToXiPiPiFulls, "AOD", "HFXICXI2PIFULL",
                  full::CandidateSelFlag,
                  full::XPv,
                  full::YPv,
                  full::ZPv,
                  full::XPvErr,
                  full::YPvErr,
                  full::ZPvErr,
                  full::XSv,
                  full::YSv,
                  full::ZSv,
                  full::Chi2Sv,
                  full::XSvErr,
                  full::YSvErr,
                  full::ZSvErr,
                  full::XDecVtxXi,
                  full::YDecVtxXi,
                  full::ZDecVtxXi,
                  full::XDecVtxLam,
                  full::YDecVtxLam,
                  full::ZDecVtxLam,
                  full::Sign,
                  full::E,
                  full::M,
                  full::P,
                  full::Pt,
                  full::Y,
                  full::Eta,
                  full::Phi,
                  full::Ct,
                  full::DecayLength,
                  full::DecayLengthNormalised,
                  full::DecayLengthXY,
                  full::DecayLengthXYNormalised,
                  full::Cpa,
                  full::CpaXY,
                  full::PtXi,
                  full::PtPi0,
                  full::PtPi1,
                  full::ImpactParameterXi,
                  full::ImpactParameterNormalisedXi,
                  full::ImpactParameterPi0,
                  full::ImpactParameterNormalisedPi0,
                  full::ImpactParameterPi1,
                  full::ImpactParameterNormalisedPi1,
                  full::MaxNormalisedDeltaIP,
                  full::CpaXi,
                  full::CpaXYXi,
                  full::CpaLam,
                  full::CpaXYLam,
                  full::InvMassXiPi0,
                  full::InvMassXiPi1,
                  hf_cand_xic_to_xi_pi_pi::FlagMcMatchRec);

DECLARE_SOA_TABLE(HfCandXicToXiPiPiFullKfs, "AOD", "HFXICXI2PIFULKF",
                  full::CandidateSelFlag,
                  full::XPv,
                  full::YPv,
                  full::ZPv,
                  full::XPvErr,
                  full::YPvErr,
                  full::ZPvErr,
                  full::XSv,
                  full::YSv,
                  full::ZSv,
                  full::Chi2Sv,
                  full::XSvErr,
                  full::YSvErr,
                  full::ZSvErr,
                  full::XDecVtxXi,
                  full::YDecVtxXi,
                  full::ZDecVtxXi,
                  full::XDecVtxLam,
                  full::YDecVtxLam,
                  full::ZDecVtxLam,
                  full::Sign,
                  full::E,
                  full::M,
                  full::P,
                  full::Pt,
                  full::Y,
                  full::Eta,
                  full::Phi,
                  full::Ct,
                  full::DecayLength,
                  full::DecayLengthNormalised,
                  full::DecayLengthXY,
                  full::DecayLengthXYNormalised,
                  full::Cpa,
                  full::CpaXY,
                  full::PtXi,
                  full::PtPi0,
                  full::PtPi1,
                  full::ImpactParameterXi,
                  full::ImpactParameterNormalisedXi,
                  full::ImpactParameterPi0,
                  full::ImpactParameterNormalisedPi0,
                  full::ImpactParameterPi1,
                  full::ImpactParameterNormalisedPi1,
                  full::MaxNormalisedDeltaIP,
                  full::CpaXi,
                  full::CpaXYXi,
                  full::CpaLam,
                  full::CpaXYLam,
                  full::InvMassXiPi0,
                  full::InvMassXiPi1,
                  full::Chi2XiVtx,
                  full::Chi2LamVtx,
                  full::Chi2XicPlusTopoToPV,
                  full::Chi2XicPlusTopoXiToXicPlus,
                  full::DcaXYPi0Pi1,
                  full::DcaXYPi0Xi,
                  full::DcaXYPi1Xi,
                  full::DcaPi0Pi1,
                  full::DcaPi0Xi,
                  full::DcaPi1Xi,
                  full::DcaXiDaughters,
                  hf_cand_xic_to_xi_pi_pi::FlagMcMatchRec);

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

  Configurable<int> selectionFlagXic{"selectionXic", 1, "Selection Flag for Xic"};
  Configurable<bool> fillCandidateLiteTable{"fillCandidateLiteTable", false, "Switch to fill lite table with candidate properties"};
  Configurable<bool> fillResidualTable{"fillResidualTable", false, "Switch to fill table with residuals for MC candidates"};
  // parameters for production of training samples
  Configurable<bool> fillOnlySignal{"fillOnlySignal", false, "Flag to fill derived tables with signal for ML trainings"};
  Configurable<bool> fillOnlyBackground{"fillOnlyBackground", false, "Flag to fill derived tables with background for ML trainings"};
  Configurable<float> downSampleBkgFactor{"downSampleBkgFactor", 1., "Fraction of background candidates to keep for ML trainings"};
  Configurable<float> ptMaxForDownSample{"ptMaxForDownSample", 10., "Maximum pt for the application of the downsampling factor"};

  using SelectedCandidates = soa::Filtered<soa::Join<aod::HfCandXic, aod::HfSelXicToXiPiPi>>;
  using SelectedCandidatesKf = soa::Filtered<soa::Join<aod::HfCandXic, aod::HfCandXicKF, aod::HfSelXicToXiPiPi>>;
  using SelectedCandidatesMc = soa::Filtered<soa::Join<aod::HfCandXic, aod::HfCandXicMcRec, aod::HfSelXicToXiPiPi>>;
  using SelectedCandidatesKfMc = soa::Filtered<soa::Join<aod::HfCandXic, aod::HfCandXicKF, aod::HfCandXicMcRec, aod::HfSelXicToXiPiPi>>;
  using TracksWPid = soa::Join<aod::Tracks, aod::TracksPidPi>;

  Filter filterSelectCandidates = aod::hf_sel_candidate_xic::isSelXicToXiPiPi >= selectionFlagXic;

  Preslice<aod::McParticles> genParticlesPerCollision = aod::mcparticle::mcCollisionId;

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
          candidate.isSelXicToXiPiPi(),
          candidate.posX(),
          candidate.posY(),
          candidate.posZ(),
          candidate.xSecondaryVertex(),
          candidate.ySecondaryVertex(),
          candidate.zSecondaryVertex(),
          candidate.chi2PCA(),
          candidate.sign(),
          candidate.e(o2::constants::physics::MassXiCPlus),
          candidate.invMassXic(),
          candidate.p(),
          candidate.pt(),
          candidate.y(o2::constants::physics::MassXiCPlus),
          candidate.eta(),
          candidate.phi(),
          candidate.ct(o2::constants::physics::MassXiCPlus),
          candidate.decayLength(),
          candidate.decayLengthXY(),
          candidate.cpa(),
          candidate.cpaXY(),
          candidate.ptProng0(),
          candidate.ptProng1(),
          candidate.ptProng2(),
          candidate.impactParameter0(),
          candidate.impactParameter1(),
          candidate.impactParameter2(),
          candidate.cosPaXi(),
          candidate.cosPaXYXi(),
          candidate.cosPaLambda(),
          candidate.cosPaXYLambda(),
          flagMc);
      } else {
        rowCandidateFull(
          candidate.isSelXicToXiPiPi(),
          candidate.posX(),
          candidate.posY(),
          candidate.posZ(),
          candidate.xPvErr(),
          candidate.yPvErr(),
          candidate.zPvErr(),
          candidate.xSecondaryVertex(),
          candidate.ySecondaryVertex(),
          candidate.zSecondaryVertex(),
          candidate.chi2PCA(),
          candidate.xSvErr(),
          candidate.ySvErr(),
          candidate.zSvErr(),
          candidate.xDecayVtxXi(),
          candidate.yDecayVtxXi(),
          candidate.zDecayVtxXi(),
          candidate.xDecayVtxLambda(),
          candidate.yDecayVtxLambda(),
          candidate.zDecayVtxLambda(),
          candidate.sign(),
          candidate.e(o2::constants::physics::MassXiCPlus),
          candidate.invMassXic(),
          candidate.p(),
          candidate.pt(),
          candidate.y(o2::constants::physics::MassXiCPlus),
          candidate.eta(),
          candidate.phi(),
          candidate.ct(o2::constants::physics::MassXiCPlus),
          candidate.decayLength(),
          candidate.decayLengthNormalised(),
          candidate.decayLengthXY(),
          candidate.decayLengthXYNormalised(),
          candidate.cpa(),
          candidate.cpaXY(),
          candidate.ptProng0(),
          candidate.ptProng1(),
          candidate.ptProng2(),
          candidate.impactParameter0(),
          candidate.impactParameterNormalised0(),
          candidate.impactParameter1(),
          candidate.impactParameterNormalised1(),
          candidate.impactParameter2(),
          candidate.impactParameterNormalised2(),
          candidate.maxNormalisedDeltaIP(),
          candidate.cosPaXi(),
          candidate.cosPaXYXi(),
          candidate.cosPaLambda(),
          candidate.cosPaXYLambda(),
          candidate.invMassXiPi0(),
          candidate.invMassXiPi1(),
          flagMc);
      }
    } else {
      if (fillCandidateLiteTable) {
        rowCandidateLiteKf(
          candidate.isSelXicToXiPiPi(),
          candidate.posX(),
          candidate.posY(),
          candidate.posZ(),
          candidate.xSecondaryVertex(),
          candidate.ySecondaryVertex(),
          candidate.zSecondaryVertex(),
          candidate.chi2PCA(),
          candidate.sign(),
          candidate.e(o2::constants::physics::MassXiCPlus),
          candidate.invMassXic(),
          candidate.p(),
          candidate.pt(),
          candidate.y(o2::constants::physics::MassXiCPlus),
          candidate.eta(),
          candidate.phi(),
          candidate.ct(o2::constants::physics::MassXiCPlus),
          candidate.decayLength(),
          candidate.decayLengthXY(),
          candidate.cpa(),
          candidate.cpaXY(),
          candidate.ptProng0(),
          candidate.ptProng1(),
          candidate.ptProng2(),
          candidate.impactParameter0(),
          candidate.impactParameter1(),
          candidate.impactParameter2(),
          candidate.cosPaXi(),
          candidate.cosPaXYXi(),
          candidate.cosPaLambda(),
          candidate.cosPaXYLambda(),
          candidate.kfCascadeChi2(),
          candidate.kfV0Chi2(),
          candidate.chi2TopoXicPlusToPV(),
          candidate.chi2TopoXiToXicPlus(),
          candidate.dcaXYPi0Pi1(),
          candidate.dcaXYPi0Xi(),
          candidate.dcaXYPi1Xi(),
          candidate.dcaPi0Pi1(),
          candidate.dcaPi0Xi(),
          candidate.dcaPi1Xi(),
          candidate.dcacascdaughters(),
          flagMc);
      } else {
        rowCandidateFullKf(
          candidate.isSelXicToXiPiPi(),
          candidate.posX(),
          candidate.posY(),
          candidate.posZ(),
          candidate.xPvErr(),
          candidate.yPvErr(),
          candidate.zPvErr(),
          candidate.xSecondaryVertex(),
          candidate.ySecondaryVertex(),
          candidate.zSecondaryVertex(),
          candidate.chi2PCA(),
          candidate.xSvErr(),
          candidate.ySvErr(),
          candidate.zSvErr(),
          candidate.xDecayVtxXi(),
          candidate.yDecayVtxXi(),
          candidate.zDecayVtxXi(),
          candidate.xDecayVtxLambda(),
          candidate.yDecayVtxLambda(),
          candidate.zDecayVtxLambda(),
          candidate.sign(),
          candidate.e(o2::constants::physics::MassXiCPlus),
          candidate.invMassXic(),
          candidate.p(),
          candidate.pt(),
          candidate.y(o2::constants::physics::MassXiCPlus),
          candidate.eta(),
          candidate.phi(),
          candidate.ct(o2::constants::physics::MassXiCPlus),
          candidate.decayLength(),
          candidate.decayLengthNormalised(),
          candidate.decayLengthXY(),
          candidate.decayLengthXYNormalised(),
          candidate.cpa(),
          candidate.cpaXY(),
          candidate.ptProng0(),
          candidate.ptProng1(),
          candidate.ptProng2(),
          candidate.impactParameter0(),
          candidate.impactParameterNormalised0(),
          candidate.impactParameter1(),
          candidate.impactParameterNormalised1(),
          candidate.impactParameter2(),
          candidate.impactParameterNormalised2(),
          candidate.maxNormalisedDeltaIP(),
          candidate.cosPaXi(),
          candidate.cosPaXYXi(),
          candidate.cosPaLambda(),
          candidate.cosPaXYLambda(),
          candidate.invMassXiPi0(),
          candidate.invMassXiPi1(),
          candidate.kfCascadeChi2(),
          candidate.kfV0Chi2(),
          candidate.chi2TopoXicPlusToPV(),
          candidate.chi2TopoXiToXicPlus(),
          candidate.dcaXYPi0Pi1(),
          candidate.dcaXYPi0Xi(),
          candidate.dcaXYPi1Xi(),
          candidate.dcaPi0Pi1(),
          candidate.dcaPi0Xi(),
          candidate.dcaPi1Xi(),
          candidate.dcacascdaughters(),
          flagMc);
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
        float pseudoRndm = candidate.ptProng1() * 1000. - (int64_t)(candidate.ptProng1() * 1000);
        if (pseudoRndm >= downSampleBkgFactor && candidate.pt() < ptMaxForDownSample) {
          continue;
        }
      }
      fillCandidateTable<false, false>(candidate);
    }
  }
  PROCESS_SWITCH(HfTreeCreatorXicToXiPiPi, processData, "Process data", true);

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
        float pseudoRndm = candidate.ptProng1() * 1000. - (int64_t)(candidate.ptProng1() * 1000);
        if (pseudoRndm >= downSampleBkgFactor && candidate.pt() < ptMaxForDownSample) {
          continue;
        }
      }
      fillCandidateTable<false, true>(candidate);
    }
  }
  PROCESS_SWITCH(HfTreeCreatorXicToXiPiPi, processDataKf, "Process data with KF Particle reconstruction", false);

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
        float pseudoRndm = candidate.ptProng1() * 1000. - (int64_t)(candidate.ptProng1() * 1000);
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

    // Filling particle properties
    rowCandidateFullParticles.reserve(particles.size());
    for (const auto& particle : particles) {
      if (TESTBIT(std::abs(particle.flagMcMatchGen()), aod::hf_cand_xic_to_xi_pi_pi::DecayType::XicToXiPiPi) || TESTBIT(std::abs(particle.flagMcMatchGen()), aod::hf_cand_xic_to_xi_pi_pi::DecayType::XicToXiResPiToXiPiPi)) {
        arrDaughIndex.clear();
        RecoDecay::getDaughters(particle, &arrDaughIndex, std::array{+kXiMinus, +kPiPlus, +kPiPlus}, 2);
        auto xicDaugh0 = particles.rawIteratorAt(arrDaughIndex[0]);

        rowCandidateFullParticles(
          particle.pt(),
          particle.eta(),
          particle.phi(),
          RecoDecay::y(particle.pVector(), o2::constants::physics::MassXiCPlus),
          particle.vx(),
          particle.vy(),
          particle.vz(),
          xicDaugh0.vx(),
          xicDaugh0.vx(),
          xicDaugh0.vz(),
          particle.flagMcMatchGen());
      }
    }

    if (fillResidualTable) {
      rowCandidateResiduals.reserve(recSig.size());
      for (const auto& candidate : recSig) {
        auto thisCollId = candidate.collisionId();
        auto groupedParticles = particles.sliceBy(genParticlesPerCollision, thisCollId);

        for (const auto& particle : groupedParticles) {
          std::array<float, 3> pvResiduals;
          std::array<float, 3> svResiduals;
          std::array<float, 3> pvPulls = {-999.9};
          std::array<float, 3> svPulls = {-999.9};

          arrDaughIndex.clear();
          RecoDecay::getDaughters(particle, &arrDaughIndex, std::array{+kXiMinus, +kPiPlus, +kPiPlus}, 2);

          if ((candidate.cascadeId() == arrDaughIndex[0] || candidate.cascadeId() == arrDaughIndex[1] || candidate.cascadeId() == arrDaughIndex[2]) && (candidate.pi0Id() == arrDaughIndex[0] || candidate.pi0Id() == arrDaughIndex[1] || candidate.pi0Id() == arrDaughIndex[2]) && (candidate.pi1Id() == arrDaughIndex[0] || candidate.pi1Id() == arrDaughIndex[1] || candidate.pi1Id() == arrDaughIndex[2])) {
            auto xicDaugh0 = particles.rawIteratorAt(arrDaughIndex[0]);

            // residuals
            float pResidual = candidate.p() - particle.p();
            float ptResidual = candidate.pt() - particle.pt();
            pvResiduals[0] = candidate.posX() - particle.vx();
            pvResiduals[1] = candidate.posY() - particle.vy();
            pvResiduals[2] = candidate.posZ() - particle.vz();
            svResiduals[0] = candidate.xSecondaryVertex() - xicDaugh0.vx();
            svResiduals[1] = candidate.ySecondaryVertex() - xicDaugh0.vy();
            svResiduals[2] = candidate.zSecondaryVertex() - xicDaugh0.vz();
            // pulls
            try {
              pvPulls[0] = pvResiduals[0] / candidate.xPvErr();
              pvPulls[1] = pvResiduals[1] / candidate.yPvErr();
              pvPulls[2] = pvResiduals[2] / candidate.zPvErr();
              svPulls[0] = svResiduals[0] / candidate.xSvErr();
              svPulls[1] = svResiduals[1] / candidate.ySvErr();
              svPulls[2] = svResiduals[2] / candidate.zSvErr();
            } catch (const std::runtime_error& error) {
              LOG(info) << "Run time error found: " << error.what() << ". Set values of vertex pulls to -999.9.";
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
          }
        } // loop over generated particles
      } // loop over reconstructed candidates
    }
  }
  PROCESS_SWITCH(HfTreeCreatorXicToXiPiPi, processMc, "Process MC", false);

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
        float pseudoRndm = candidate.ptProng1() * 1000. - (int64_t)(candidate.ptProng1() * 1000);
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

    // Filling particle properties
    rowCandidateFullParticles.reserve(particles.size());
    for (const auto& particle : particles) {
      if (TESTBIT(std::abs(particle.flagMcMatchGen()), aod::hf_cand_xic_to_xi_pi_pi::DecayType::XicToXiPiPi) || TESTBIT(std::abs(particle.flagMcMatchGen()), aod::hf_cand_xic_to_xi_pi_pi::DecayType::XicToXiResPiToXiPiPi)) {
        arrDaughIndex.clear();
        RecoDecay::getDaughters(particle, &arrDaughIndex, std::array{+kXiMinus, +kPiPlus, +kPiPlus}, 2);
        auto xicDaugh0 = particles.rawIteratorAt(arrDaughIndex[0]);

        rowCandidateFullParticles(
          particle.pt(),
          particle.eta(),
          particle.phi(),
          RecoDecay::y(particle.pVector(), o2::constants::physics::MassXiCPlus),
          particle.vx(),
          particle.vy(),
          particle.vz(),
          xicDaugh0.vx(),
          xicDaugh0.vx(),
          xicDaugh0.vz(),
          particle.flagMcMatchGen());
      }
    }

    if (fillResidualTable) {
      rowCandidateResiduals.reserve(recSigKf.size());
      for (const auto& candidate : recSigKf) {
        auto thisCollId = candidate.collisionId();
        auto groupedParticles = particles.sliceBy(genParticlesPerCollision, thisCollId);

        for (const auto& particle : groupedParticles) {
          std::array<float, 3> pvResiduals;
          std::array<float, 3> svResiduals;
          std::array<float, 3> pvPulls = {-999.9};
          std::array<float, 3> svPulls = {-999.9};

          arrDaughIndex.clear();
          RecoDecay::getDaughters(particle, &arrDaughIndex, std::array{+kXiMinus, +kPiPlus, +kPiPlus}, 2);

          if ((candidate.cascadeId() == arrDaughIndex[0] || candidate.cascadeId() == arrDaughIndex[1] || candidate.cascadeId() == arrDaughIndex[2]) && (candidate.pi0Id() == arrDaughIndex[0] || candidate.pi0Id() == arrDaughIndex[1] || candidate.pi0Id() == arrDaughIndex[2]) && (candidate.pi1Id() == arrDaughIndex[0] || candidate.pi1Id() == arrDaughIndex[1] || candidate.pi1Id() == arrDaughIndex[2])) {
            auto xicDaugh0 = particles.rawIteratorAt(arrDaughIndex[0]);

            // residuals
            float pResidual = candidate.p() - particle.p();
            float ptResidual = candidate.pt() - particle.pt();
            pvResiduals[0] = candidate.posX() - particle.vx();
            pvResiduals[1] = candidate.posY() - particle.vy();
            pvResiduals[2] = candidate.posZ() - particle.vz();
            svResiduals[0] = candidate.xSecondaryVertex() - xicDaugh0.vx();
            svResiduals[1] = candidate.ySecondaryVertex() - xicDaugh0.vy();
            svResiduals[2] = candidate.zSecondaryVertex() - xicDaugh0.vz();
            // pulls
            try {
              pvPulls[0] = pvResiduals[0] / candidate.xPvErr();
              pvPulls[1] = pvResiduals[1] / candidate.yPvErr();
              pvPulls[2] = pvResiduals[2] / candidate.zPvErr();
              svPulls[0] = svResiduals[0] / candidate.xSvErr();
              svPulls[1] = svResiduals[1] / candidate.ySvErr();
              svPulls[2] = svResiduals[2] / candidate.zSvErr();
            } catch (const std::runtime_error& error) {
              LOG(info) << "Run time error found: " << error.what() << ". Set values of vertex pulls to -999.9.";
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
          }
        } // loop over generated particles
      } // loop over reconstructed candidates
    }
  }
  PROCESS_SWITCH(HfTreeCreatorXicToXiPiPi, processMcKf, "Process MC with KF Particle reconstruction", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTreeCreatorXicToXiPiPi>(cfgc)};
}
