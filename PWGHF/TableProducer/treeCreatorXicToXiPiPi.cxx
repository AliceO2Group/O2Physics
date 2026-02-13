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
/// \author Jaeyoon Cho <jaeyoon.cho@cern.ch>, Inha University

#include "PWGHF/Core/DecayChannelsLegacy.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <Rtypes.h>

#include <cstdint>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

namespace o2::aod
{
namespace full
{
DECLARE_SOA_COLUMN(ParticleFlag, particleFlag, int8_t);      //! hf_cand_xic_to_xi_pi_pi::Sign for data, hf_cand_xic_to_xi_pi_pi::FlagMcMatchRec for MC
DECLARE_SOA_COLUMN(CandidateSelFlag, candidateSelFlag, int); //! Selection flag of candidate (output of candidateSelector)
// vertices
DECLARE_SOA_COLUMN(Chi2SV, chi2SV, float);               //! Chi2 of candidate vertex
DECLARE_SOA_COLUMN(Chi2GeoXi, chi2GeoXi, float);         //! Chi2 of Xi vertex
DECLARE_SOA_COLUMN(Chi2GeoLambda, chi2GeoLambda, float); //! Chi2 of Lambda vertex
// properties of XicPlus
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
// properties of daughter tracks
DECLARE_SOA_COLUMN(PtXi, ptXi, float);                                                 //! Transverse momentum of Xi (prong0) (GeV/c)
DECLARE_SOA_COLUMN(ImpactParameterXi, impactParameterXi, float);                       //! Impact parameter of Xi (prong0)
DECLARE_SOA_COLUMN(ImpactParameterNormalisedXi, impactParameterNormalisedXi, float);   //! Normalised impact parameter of Xi (prong0)
DECLARE_SOA_COLUMN(PPi0, pPi0, float);                                                 //! Momentum of Pi0 (prong1) (GeV/c)
DECLARE_SOA_COLUMN(PtPi0, ptPi0, float);                                               //! Transverse momentum of Pi0 (prong1) (GeV/c)
DECLARE_SOA_COLUMN(ImpactParameterPi0, impactParameterPi0, float);                     //! Impact parameter of Pi0 (prong1)
DECLARE_SOA_COLUMN(ImpactParameterNormalisedPi0, impactParameterNormalisedPi0, float); //! Normalised impact parameter of Pi0 (prong1)
DECLARE_SOA_COLUMN(PPi1, pPi1, float);                                                 //! Momentum of Pi1 (prong2) (GeV/c)
DECLARE_SOA_COLUMN(PtPi1, ptPi1, float);                                               //! Transverse momentum of Pi1 (prong2) (GeV/c)
DECLARE_SOA_COLUMN(ImpactParameterPi1, impactParameterPi1, float);                     //! Normalised impact parameter of Pi1 (prong2)
DECLARE_SOA_COLUMN(ImpactParameterNormalisedPi1, impactParameterNormalisedPi1, float); //! Normalised impact parameter of Pi1 (prong2)
DECLARE_SOA_COLUMN(MaxNormalisedDeltaIP, maxNormalisedDeltaIP, float);                 //! Maximum normalized difference between measured and expected impact parameter of candidate prongs
} // namespace full

DECLARE_SOA_TABLE(HfCandXicToXiPiPiLites, "AOD", "HFXICXI2PILITE",
                  full::ParticleFlag,
                  hf_cand_xic_to_xi_pi_pi::OriginMcRec,
                  full::CandidateSelFlag,
                  full::Y,
                  full::Eta,
                  full::Phi,
                  full::P,
                  full::Pt,
                  full::PtXi,
                  full::PtPi0,
                  full::PtPi1,
                  full::M,
                  hf_cand_xic_to_xi_pi_pi::InvMassXi,
                  hf_cand_xic_to_xi_pi_pi::InvMassLambda,
                  hf_cand_xic_to_xi_pi_pi::InvMassXiPi0,
                  hf_cand_xic_to_xi_pi_pi::InvMassXiPi1,
                  full::Chi2SV,
                  full::Ct,
                  full::DecayLength,
                  full::DecayLengthNormalised,
                  full::DecayLengthXY,
                  full::DecayLengthXYNormalised,
                  full::Cpa,
                  full::CpaXY,
                  hf_cand_xic_to_xi_pi_pi::CpaXi,
                  hf_cand_xic_to_xi_pi_pi::CpaXYXi,
                  hf_cand_xic_to_xi_pi_pi::CpaLambda,
                  hf_cand_xic_to_xi_pi_pi::CpaXYLambda,
                  full::ImpactParameterXi,
                  full::ImpactParameterNormalisedXi,
                  full::ImpactParameterPi0,
                  full::ImpactParameterNormalisedPi0,
                  full::ImpactParameterPi1,
                  full::ImpactParameterNormalisedPi1,
                  full::MaxNormalisedDeltaIP);

DECLARE_SOA_TABLE(HfCandXicToXiPiPiLiteKfs, "AOD", "HFXICXI2PILITKF",
                  full::ParticleFlag,
                  hf_cand_xic_to_xi_pi_pi::OriginMcRec,
                  full::CandidateSelFlag,
                  full::Y,
                  full::Eta,
                  full::Phi,
                  full::Pt,
                  full::PtXi,
                  full::PtPi0,
                  full::PtPi1,
                  full::M,
                  full::Chi2SV,
                  full::Ct,
                  full::DecayLength,
                  full::DecayLengthNormalised,
                  full::DecayLengthXY,
                  full::DecayLengthXYNormalised,
                  full::Cpa,
                  full::CpaXY,
                  hf_cand_xic_to_xi_pi_pi::CpaXi,
                  hf_cand_xic_to_xi_pi_pi::CpaXYXi,
                  hf_cand_xic_to_xi_pi_pi::CpaLambda,
                  hf_cand_xic_to_xi_pi_pi::CpaXYLambda,
                  hf_cand_xic_to_xi_pi_pi::CpaLambdaToXi,
                  hf_cand_xic_to_xi_pi_pi::CpaXYLambdaToXi,
                  full::ImpactParameterXi,
                  full::ImpactParameterNormalisedXi,
                  full::ImpactParameterPi0,
                  full::ImpactParameterNormalisedPi0,
                  full::ImpactParameterPi1,
                  full::ImpactParameterNormalisedPi1,
                  full::MaxNormalisedDeltaIP,
                  hf_cand_xic_to_xi_pi_pi::DcaXiDaughters,
                  hf_cand_xic_to_xi_pi_pi::DcaV0Daughters,
                  hf_cand_xic_to_xi_pi_pi::DcaXYCascToPV,
                  hf_cand_xic_to_xi_pi_pi::DcaZCascToPV,
                  hf_cand_xic_to_xi_pi_pi::DcaBachelorToPV,
                  hf_cand_xic_to_xi_pi_pi::DcaPosToPV,
                  hf_cand_xic_to_xi_pi_pi::DcaNegToPV,
                  // PID information
                  hf_cand_xic_to_xi_pi_pi::NSigTpcPiFromXicPlus0,
                  hf_cand_xic_to_xi_pi_pi::NSigTpcPiFromXicPlus1,
                  hf_cand_xic_to_xi_pi_pi::NSigTpcBachelorPi,
                  hf_cand_xic_to_xi_pi_pi::NSigTpcPiFromLambda,
                  hf_cand_xic_to_xi_pi_pi::NSigTpcPrFromLambda,
                  hf_cand_xic_to_xi_pi_pi::NSigTofPiFromXicPlus0,
                  hf_cand_xic_to_xi_pi_pi::NSigTofPiFromXicPlus1,
                  hf_cand_xic_to_xi_pi_pi::NSigTofBachelorPi,
                  hf_cand_xic_to_xi_pi_pi::NSigTofPiFromLambda,
                  hf_cand_xic_to_xi_pi_pi::NSigTofPrFromLambda,
                  // KF specific columns
                  full::Chi2GeoXi,
                  full::Chi2GeoLambda,
                  hf_cand_xic_to_xi_pi_pi::Chi2TopoXicPlusToPVBefConst,
                  hf_cand_xic_to_xi_pi_pi::Chi2TopoXicPlusToPV,
                  hf_cand_xic_to_xi_pi_pi::Chi2PrimXi,
                  hf_cand_xic_to_xi_pi_pi::Chi2PrimPi0,
                  hf_cand_xic_to_xi_pi_pi::Chi2PrimPi1,
                  hf_cand_xic_to_xi_pi_pi::Chi2DevPi0Pi1,
                  hf_cand_xic_to_xi_pi_pi::Chi2DevPi0Xi,
                  hf_cand_xic_to_xi_pi_pi::Chi2DevPi1Xi,
                  hf_cand_xic_to_xi_pi_pi::DcaPi0Pi1,
                  hf_cand_xic_to_xi_pi_pi::DcaPi0Xi,
                  hf_cand_xic_to_xi_pi_pi::DcaPi1Xi,
                  hf_cand_xic_to_xi_pi_pi::DcaXYPi0Pi1,
                  hf_cand_xic_to_xi_pi_pi::DcaXYPi0Xi,
                  hf_cand_xic_to_xi_pi_pi::DcaXYPi1Xi);

DECLARE_SOA_TABLE(HfCandXicToXiPiPiFulls, "AOD", "HFXICXI2PIFULL",
                  full::ParticleFlag,
                  hf_cand_xic_to_xi_pi_pi::OriginMcRec,
                  full::CandidateSelFlag,
                  full::Y,
                  full::Eta,
                  full::Phi,
                  full::P,
                  full::Pt,
                  full::PtXi,
                  full::PtPi0,
                  full::PtPi1,
                  full::M,
                  hf_cand_xic_to_xi_pi_pi::InvMassXi,
                  hf_cand_xic_to_xi_pi_pi::InvMassLambda,
                  hf_cand_xic_to_xi_pi_pi::InvMassXiPi0,
                  hf_cand_xic_to_xi_pi_pi::InvMassXiPi1,
                  full::Chi2SV,
                  full::Ct,
                  full::DecayLength,
                  full::DecayLengthNormalised,
                  full::DecayLengthXY,
                  full::DecayLengthXYNormalised,
                  full::Cpa,
                  full::CpaXY,
                  hf_cand_xic_to_xi_pi_pi::CpaXi,
                  hf_cand_xic_to_xi_pi_pi::CpaXYXi,
                  hf_cand_xic_to_xi_pi_pi::CpaLambda,
                  hf_cand_xic_to_xi_pi_pi::CpaXYLambda,
                  full::ImpactParameterXi,
                  full::ImpactParameterNormalisedXi,
                  full::ImpactParameterPi0,
                  full::ImpactParameterNormalisedPi0,
                  full::ImpactParameterPi1,
                  full::ImpactParameterNormalisedPi1,
                  full::MaxNormalisedDeltaIP,
                  // additional columns only stored in the full candidate table
                  hf_cand_xic_to_xi_pi_pi::CpaLambdaToXi,
                  hf_cand_xic_to_xi_pi_pi::CpaXYLambdaToXi,
                  full::PPi0,
                  full::PPi1,
                  hf_cand_xic_to_xi_pi_pi::PBachelorPi,
                  hf_cand_xic_to_xi_pi_pi::PPiFromLambda,
                  hf_cand_xic_to_xi_pi_pi::PPrFromLambda,
                  hf_cand_xic_to_xi_pi_pi::DcaXiDaughters,
                  hf_cand_xic_to_xi_pi_pi::DcaV0Daughters,
                  hf_cand_xic_to_xi_pi_pi::DcaPosToPV,
                  hf_cand_xic_to_xi_pi_pi::DcaNegToPV,
                  hf_cand_xic_to_xi_pi_pi::DcaBachelorToPV,
                  hf_cand_xic_to_xi_pi_pi::DcaXYCascToPV,
                  hf_cand_xic_to_xi_pi_pi::DcaZCascToPV,
                  hf_cand_xic_to_xi_pi_pi::NSigTpcPiFromXicPlus0,
                  hf_cand_xic_to_xi_pi_pi::NSigTpcPiFromXicPlus1,
                  hf_cand_xic_to_xi_pi_pi::NSigTpcBachelorPi,
                  hf_cand_xic_to_xi_pi_pi::NSigTpcPiFromLambda,
                  hf_cand_xic_to_xi_pi_pi::NSigTpcPrFromLambda,
                  hf_cand_xic_to_xi_pi_pi::NSigTofPiFromXicPlus0,
                  hf_cand_xic_to_xi_pi_pi::NSigTofPiFromXicPlus1,
                  hf_cand_xic_to_xi_pi_pi::NSigTofBachelorPi,
                  hf_cand_xic_to_xi_pi_pi::NSigTofPiFromLambda,
                  hf_cand_xic_to_xi_pi_pi::NSigTofPrFromLambda);

DECLARE_SOA_TABLE(HfCandXicToXiPiPiFullKfs, "AOD", "HFXICXI2PIFULKF",
                  full::ParticleFlag,
                  hf_cand_xic_to_xi_pi_pi::OriginMcRec,
                  full::CandidateSelFlag,
                  full::Y,
                  full::Eta,
                  full::Phi,
                  full::Pt,
                  full::PtXi,
                  full::PtPi0,
                  full::PtPi1,
                  full::M,
                  full::Chi2SV,
                  full::Ct,
                  full::DecayLength,
                  full::DecayLengthNormalised,
                  full::DecayLengthXY,
                  full::DecayLengthXYNormalised,
                  full::Cpa,
                  full::CpaXY,
                  hf_cand_xic_to_xi_pi_pi::CpaXi,
                  hf_cand_xic_to_xi_pi_pi::CpaXYXi,
                  hf_cand_xic_to_xi_pi_pi::CpaLambda,
                  hf_cand_xic_to_xi_pi_pi::CpaXYLambda,
                  hf_cand_xic_to_xi_pi_pi::CpaLambdaToXi,
                  hf_cand_xic_to_xi_pi_pi::CpaXYLambdaToXi,
                  full::ImpactParameterXi,
                  full::ImpactParameterNormalisedXi,
                  full::ImpactParameterPi0,
                  full::ImpactParameterNormalisedPi0,
                  full::ImpactParameterPi1,
                  full::ImpactParameterNormalisedPi1,
                  full::MaxNormalisedDeltaIP,
                  hf_cand_xic_to_xi_pi_pi::DcaXiDaughters,
                  hf_cand_xic_to_xi_pi_pi::DcaV0Daughters,
                  hf_cand_xic_to_xi_pi_pi::DcaXYCascToPV,
                  hf_cand_xic_to_xi_pi_pi::DcaZCascToPV,
                  hf_cand_xic_to_xi_pi_pi::DcaBachelorToPV,
                  hf_cand_xic_to_xi_pi_pi::DcaPosToPV,
                  hf_cand_xic_to_xi_pi_pi::DcaNegToPV,
                  // additional columns only stored in the full candidate table
                  hf_cand_xic_to_xi_pi_pi::InvMassXi,
                  hf_cand_xic_to_xi_pi_pi::InvMassXiPi0,
                  hf_cand_xic_to_xi_pi_pi::InvMassXiPi1,
                  hf_cand_xic_to_xi_pi_pi::InvMassLambda,
                  full::P,
                  full::PPi0,
                  full::PPi1,
                  hf_cand_xic_to_xi_pi_pi::PBachelorPi,
                  hf_cand_xic_to_xi_pi_pi::PPiFromLambda,
                  hf_cand_xic_to_xi_pi_pi::PPrFromLambda,
                  // PID information
                  hf_cand_xic_to_xi_pi_pi::NSigTpcPiFromXicPlus0,
                  hf_cand_xic_to_xi_pi_pi::NSigTpcPiFromXicPlus1,
                  hf_cand_xic_to_xi_pi_pi::NSigTpcBachelorPi,
                  hf_cand_xic_to_xi_pi_pi::NSigTpcPiFromLambda,
                  hf_cand_xic_to_xi_pi_pi::NSigTpcPrFromLambda,
                  hf_cand_xic_to_xi_pi_pi::NSigTofPiFromXicPlus0,
                  hf_cand_xic_to_xi_pi_pi::NSigTofPiFromXicPlus1,
                  hf_cand_xic_to_xi_pi_pi::NSigTofBachelorPi,
                  hf_cand_xic_to_xi_pi_pi::NSigTofPiFromLambda,
                  hf_cand_xic_to_xi_pi_pi::NSigTofPrFromLambda,
                  // KF-specific columns
                  full::Chi2GeoXi,
                  full::Chi2GeoLambda,
                  hf_cand_xic_to_xi_pi_pi::Chi2TopoXicPlusToPVBefConst,
                  hf_cand_xic_to_xi_pi_pi::Chi2TopoXicPlusToPV,
                  hf_cand_xic_to_xi_pi_pi::Chi2PrimXi,
                  hf_cand_xic_to_xi_pi_pi::Chi2PrimPi0,
                  hf_cand_xic_to_xi_pi_pi::Chi2PrimPi1,
                  hf_cand_xic_to_xi_pi_pi::Chi2DevPi0Pi1,
                  hf_cand_xic_to_xi_pi_pi::Chi2DevPi0Xi,
                  hf_cand_xic_to_xi_pi_pi::Chi2DevPi1Xi,
                  hf_cand_xic_to_xi_pi_pi::DcaPi0Pi1,
                  hf_cand_xic_to_xi_pi_pi::DcaPi0Xi,
                  hf_cand_xic_to_xi_pi_pi::DcaPi1Xi,
                  hf_cand_xic_to_xi_pi_pi::DcaXYPi0Pi1,
                  hf_cand_xic_to_xi_pi_pi::DcaXYPi0Xi,
                  hf_cand_xic_to_xi_pi_pi::DcaXYPi1Xi);

DECLARE_SOA_TABLE(HfCandXicToXiPiPiFullPs, "AOD", "HFXICXI2PIFULLP",
                  hf_cand_xic_to_xi_pi_pi::FlagMcMatchGen,
                  hf_cand_xic_to_xi_pi_pi::OriginMcGen,
                  hf_cand_mc_flag::PdgBhadMotherPart,
                  full::Pt,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  hf_cand_xic_to_xi_pi_pi::DecayLengthMcGen);
} // namespace o2::aod

/// Writes the full information in an output TTree
struct HfTreeCreatorXicToXiPiPi {
  Produces<o2::aod::HfCandXicToXiPiPiLites> rowCandidateLite;
  Produces<o2::aod::HfCandXicToXiPiPiLiteKfs> rowCandidateLiteKf;
  Produces<o2::aod::HfCandXicToXiPiPiFulls> rowCandidateFull;
  Produces<o2::aod::HfCandXicToXiPiPiFullKfs> rowCandidateFullKf;
  Produces<o2::aod::HfCandXicToXiPiPiFullPs> rowCandidateFullParticles;

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
  using MatchedGenXicToXiPiPi = soa::Filtered<soa::Join<aod::McParticles, aod::HfCandXicMcGen>>;

  Filter filterSelectCandidates = aod::hf_sel_candidate_xic::isSelXicToXiPiPi >= selectionFlagXic;
  Filter filterGenXicToXiPiPi = (nabs(aod::hf_cand_xic_to_xi_pi_pi::flagMcMatchGen) == static_cast<int8_t>(BIT(aod::hf_cand_xic_to_xi_pi_pi::DecayType::XicToXiPiPi)) || nabs(aod::hf_cand_xic_to_xi_pi_pi::flagMcMatchGen) == static_cast<int8_t>(BIT(aod::hf_cand_xic_to_xi_pi_pi::DecayType::XicToXiResPiToXiPiPi)));

  Partition<SelectedCandidatesMc> recSig = nabs(aod::hf_cand_xic_to_xi_pi_pi::flagMcMatchRec) != int8_t(0);
  Partition<SelectedCandidatesMc> recBg = nabs(aod::hf_cand_xic_to_xi_pi_pi::flagMcMatchRec) == int8_t(0);
  Partition<SelectedCandidatesKfMc> recSigKf = nabs(aod::hf_cand_xic_to_xi_pi_pi::flagMcMatchRec) != int8_t(0);
  Partition<SelectedCandidatesKfMc> recBgKf = nabs(aod::hf_cand_xic_to_xi_pi_pi::flagMcMatchRec) == int8_t(0);

  void init(InitContext const&)
  {
  }

  template <bool DoMc, bool DoKf, typename T>
  void fillCandidateTable(const T& candidate)
  {
    int8_t particleFlag = candidate.sign();
    int8_t originMc = 0;
    if constexpr (DoMc) {
      particleFlag = candidate.flagMcMatchRec();
      originMc = candidate.originMcRec();
    }
    if constexpr (!DoKf) {
      if (fillCandidateLiteTable) {
        rowCandidateLite(
          particleFlag,
          originMc,
          candidate.isSelXicToXiPiPi(),
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
          candidate.cpaXi(),
          candidate.cpaXYXi(),
          candidate.cpaLambda(),
          candidate.cpaXYLambda(),
          candidate.impactParameter0(),
          candidate.impactParameterNormalised0(),
          candidate.impactParameter1(),
          candidate.impactParameterNormalised1(),
          candidate.impactParameter2(),
          candidate.impactParameterNormalised2(),
          candidate.maxNormalisedDeltaIP());
      } else {
        rowCandidateFull(
          particleFlag,
          originMc,
          candidate.isSelXicToXiPiPi(),
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
          candidate.cpaXi(),
          candidate.cpaXYXi(),
          candidate.cpaLambda(),
          candidate.cpaXYLambda(),
          candidate.impactParameter0(),
          candidate.impactParameterNormalised0(),
          candidate.impactParameter1(),
          candidate.impactParameterNormalised1(),
          candidate.impactParameter2(),
          candidate.impactParameterNormalised2(),
          candidate.maxNormalisedDeltaIP(),
          // additional columns only stored in the full candidate table
          candidate.cpaLambdaToXi(),
          candidate.cpaXYLambdaToXi(),
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
          candidate.nSigTofPrFromLambda());
      }
    } else {
      if (fillCandidateLiteTable) {
        rowCandidateLiteKf(
          particleFlag,
          originMc,
          candidate.isSelXicToXiPiPi(),
          candidate.y(o2::constants::physics::MassXiCPlus),
          candidate.eta(),
          candidate.phi(),
          candidate.pt(),
          candidate.ptProng0(),
          candidate.ptProng1(),
          candidate.ptProng2(),
          candidate.invMassXicPlus(),
          candidate.chi2PCA(),
          candidate.ct(o2::constants::physics::MassXiCPlus),
          candidate.kfDecayLength(),
          candidate.kfDecayLengthNormalised(),
          candidate.kfDecayLengthXY(),
          candidate.kfDecayLengthXYNormalised(),
          candidate.cpa(),
          candidate.cpaXY(),
          candidate.cpaXi(),
          candidate.cpaXYXi(),
          candidate.cpaLambda(),
          candidate.cpaXYLambda(),
          candidate.cpaLambdaToXi(),
          candidate.cpaXYLambdaToXi(),
          candidate.impactParameter0(),
          candidate.impactParameterNormalised0(),
          candidate.impactParameter1(),
          candidate.impactParameterNormalised1(),
          candidate.impactParameter2(),
          candidate.impactParameterNormalised2(),
          candidate.maxNormalisedDeltaIP(),
          candidate.dcaXiDaughters(),
          candidate.dcaV0Daughters(),
          candidate.dcaXYCascToPV(),
          candidate.dcaZCascToPV(),
          candidate.dcaBachelorToPV(),
          candidate.dcaPosToPV(),
          candidate.dcaNegToPV(),
          // PID information
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
          candidate.chi2TopoXicPlusToPVBefConst(),
          candidate.chi2TopoXicPlusToPV(),
          candidate.chi2PrimXi(),
          candidate.chi2PrimPi0(),
          candidate.chi2PrimPi1(),
          candidate.chi2DevPi0Pi1(),
          candidate.chi2DevPi0Xi(),
          candidate.chi2DevPi1Xi(),
          candidate.dcaPi0Pi1(),
          candidate.dcaPi0Xi(),
          candidate.dcaPi1Xi(),
          candidate.dcaXYPi0Pi1(),
          candidate.dcaXYPi0Xi(),
          candidate.dcaXYPi1Xi());
      } else {
        rowCandidateFullKf(
          particleFlag,
          originMc,
          candidate.isSelXicToXiPiPi(),
          candidate.y(o2::constants::physics::MassXiCPlus),
          candidate.eta(),
          candidate.phi(),
          candidate.pt(),
          candidate.ptProng0(),
          candidate.ptProng1(),
          candidate.ptProng2(),
          candidate.invMassXicPlus(),
          candidate.chi2PCA(),
          candidate.ct(o2::constants::physics::MassXiCPlus),
          candidate.kfDecayLength(),
          candidate.kfDecayLengthNormalised(),
          candidate.kfDecayLengthXY(),
          candidate.kfDecayLengthXYNormalised(),
          candidate.cpa(),
          candidate.cpaXY(),
          candidate.cpaXi(),
          candidate.cpaXYXi(),
          candidate.cpaLambda(),
          candidate.cpaXYLambda(),
          candidate.cpaLambdaToXi(),
          candidate.cpaXYLambdaToXi(),
          candidate.impactParameter0(),
          candidate.impactParameterNormalised0(),
          candidate.impactParameter1(),
          candidate.impactParameterNormalised1(),
          candidate.impactParameter2(),
          candidate.impactParameterNormalised2(),
          candidate.maxNormalisedDeltaIP(),
          candidate.dcaXiDaughters(),
          candidate.dcaV0Daughters(),
          candidate.dcaXYCascToPV(),
          candidate.dcaZCascToPV(),
          candidate.dcaBachelorToPV(),
          candidate.dcaPosToPV(),
          candidate.dcaNegToPV(),
          // additional columns only stored in the full candidate table
          candidate.invMassXi(),
          candidate.invMassXiPi0(),
          candidate.invMassXiPi1(),
          candidate.invMassLambda(),
          candidate.p(),
          candidate.pProng1(),
          candidate.pProng2(),
          candidate.pBachelorPi(),
          candidate.pPiFromLambda(),
          candidate.pPrFromLambda(),
          // PID information
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
          candidate.chi2TopoXicPlusToPVBefConst(),
          candidate.chi2TopoXicPlusToPV(),
          candidate.chi2PrimXi(),
          candidate.chi2PrimPi0(),
          candidate.chi2PrimPi1(),
          candidate.chi2DevPi0Pi1(),
          candidate.chi2DevPi0Xi(),
          candidate.chi2DevPi1Xi(),
          candidate.dcaPi0Pi1(),
          candidate.dcaPi0Xi(),
          candidate.dcaPi1Xi(),
          candidate.dcaXYPi0Pi1(),
          candidate.dcaXYPi0Xi(),
          candidate.dcaXYPi1Xi());
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
        float const pseudoRndm = candidate.ptProng1() * 1000. - static_cast<int64_t>(candidate.ptProng1() * 1000);
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
        float const pseudoRndm = candidate.ptProng1() * 1000. - static_cast<int64_t>(candidate.ptProng1() * 1000);
        if (pseudoRndm >= downSampleBkgFactor && candidate.pt() < ptMaxForDownSample) {
          continue;
        }
      }
      fillCandidateTable<false, true>(candidate);
    }
  }
  PROCESS_SWITCH(HfTreeCreatorXicToXiPiPi, processDataKf, "Process data with KFParticle reconstruction", false);

  void processMc(SelectedCandidatesMc const& candidates,
                 MatchedGenXicToXiPiPi const& particles)
  {
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
        float const pseudoRndm = candidate.ptProng1() * 1000. - static_cast<int64_t>(candidate.ptProng1() * 1000);
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
        rowCandidateFullParticles(
          particle.flagMcMatchGen(),
          particle.originMcGen(),
          particle.pdgBhadMotherPart(),
          particle.pt(),
          particle.eta(),
          particle.phi(),
          RecoDecay::y(particle.pVector(), o2::constants::physics::MassXiCPlus),
          particle.decayLengthMcGen());
      }
    }
  }
  PROCESS_SWITCH(HfTreeCreatorXicToXiPiPi, processMc, "Process MC with DCAFitter reconstruction", false);

  void processMcKf(SelectedCandidatesKfMc const& candidates,
                   MatchedGenXicToXiPiPi const& particles)
  {
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
        float const pseudoRndm = candidate.ptProng1() * 1000. - static_cast<int64_t>(candidate.ptProng1() * 1000);
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
        rowCandidateFullParticles(
          particle.flagMcMatchGen(),
          particle.originMcGen(),
          particle.pdgBhadMotherPart(),
          particle.pt(),
          particle.eta(),
          particle.phi(),
          RecoDecay::y(particle.pVector(), o2::constants::physics::MassXiCPlus),
          particle.decayLengthMcGen());
      }
    }
  }
  PROCESS_SWITCH(HfTreeCreatorXicToXiPiPi, processMcKf, "Process MC with KF Particle reconstruction", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTreeCreatorXicToXiPiPi>(cfgc)};
}
